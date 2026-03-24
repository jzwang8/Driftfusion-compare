function run_sweep_chunk_profiles(chunk_id)
%RUN_SWEEP_CHUNK Run a specific chunk of the parameter sweep
%   chunk_id: which chunk to run (0-based from SLURM array)
%
%   This version:
%   - Sweeps Phi_R relative to each N_D's Fermi level
%   - Sweeps emitter doping profile shape (uniform, gaussian, exponential)
%   - EF_vals represents PEAK doping; Phi_R is relative to surface EF
%   - Skips runs where identical parameters already exist (not by index)
%   - Safe to add new parameter values without re-running old combinations
%   - Computes EQE at V=0; reverse-bias spectral response when SR_bias_vals is non-empty
%   - Supports spectral (TMM) or constant reflection via config.reflection_model
%
%   Parameter definitions live in sweep_config_profiles.m (single source of truth).
%   Combination generation lives in generate_sweep_combinations.m.
%   Filename generation lives in generate_sweep_filename.m.

%% Setup
delete(gcp('nocreate'));

% Suppress cosmetic warning from Driftfusion's internal readtable calls
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');

% Add necessary paths
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Transfer matrix modeling'));
addpath(fullfile(thisDir, 'Doping'));

%% Load shared configuration (single source of truth)
config = sweep_config_profiles();

file_name       = config.file_name;
N_sublayers     = config.N_sublayers;
SR_bias_vals    = config.SR_bias_vals;
reflection      = config.reflection;
reflection_model = config.reflection_model;
light_intensity_JV = config.light_intensity_JV;
Vmin_JV         = config.Vmin_JV;
Vmax_JV         = config.Vmax_JV;
scan_rate_JV    = config.scan_rate_JV;
cycles_JV       = config.cycles_JV;
tpoints_JV      = config.tpoints_JV;
scan_range_1    = config.scan_range_1;
scan_range_2    = config.scan_range_2;
number_of_workers = config.number_of_workers;
outputDir       = config.outputDir;
runs_per_chunk  = config.runs_per_chunk;
profile_shapes  = config.profile_shapes;

%% Setup output directory
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

existingFiles = dir(fullfile(outputDir, 'SF*.mat'));
fprintf('Found %d existing result files\n', length(existingFiles));

%% Generate all parameter combinations (shared logic)
[combinations, combo_profile_names, ~, n_skipped_invalid] = ...
    generate_sweep_combinations(config);

n_baseline = sum(combinations(:, 1) == 0);
n_with_SF  = sum(combinations(:, 1) == 1);
total_combos = size(combinations, 1);

fprintf('\nParameter space:\n');
fprintf('  Profile shapes: %s\n', strjoin(profile_shapes, ', '));
fprintf('  SF_flag=0 (baseline): %d combinations\n', n_baseline);
fprintf('  SF_flag=1 (with SF): %d combinations\n', n_with_SF);
fprintf('  Total valid: %d combinations\n', total_combos);
fprintf('  Skipped (Phi_R > E_C): %d invalid combinations\n', n_skipped_invalid);
if isempty(SR_bias_vals)
    fprintf('  Spectral response bias points: none (EQE at V=0 only)\n');
else
    fprintf('  Spectral response bias points: %s V\n', mat2str(SR_bias_vals));
end
fprintf('  Reflection model: %s\n', reflection_model);

%% Determine which runs this chunk should do
start_idx = chunk_id * runs_per_chunk + 1;
end_idx = min((chunk_id + 1) * runs_per_chunk, total_combos);

if start_idx > total_combos
    fprintf('Chunk %d: No work to do (start_idx %d > total %d)\n', ...
            chunk_id, start_idx, total_combos);
    return;
end

fprintf('\nChunk %d: Assigned indices %d to %d (of %d total)\n', ...
        chunk_id, start_idx, end_idx, total_combos);

%% Initialize Driftfusion
initialise_df;

%% Pre-compute TMM reflection spectra (only if using spectral model)
if strcmp(reflection_model, 'spectral')
    fprintf('Pre-computing TMM reflection spectra...\n');
    [lambda_tmm_SF,   Reflection_tmm_SF]   = get_TMM_reflection(true);   % Tc/ZnPc/AlOx/Si
    [lambda_tmm_noSF, Reflection_tmm_noSF] = get_TMM_reflection(false);  % bare AlOx/Si
    fprintf('  SF stack:   %d-%d nm (%d pts)\n', lambda_tmm_SF(1), lambda_tmm_SF(end), numel(lambda_tmm_SF));
    fprintf('  bare Si:    %d-%d nm (%d pts)\n', lambda_tmm_noSF(1), lambda_tmm_noSF(end), numel(lambda_tmm_noSF));
else
    fprintf('Using constant reflection kappa = %.3f for EQE\n', reflection);
end

%% Create parallel pool once before the main loop
pool = gcp('nocreate');
if isempty(pool)
    parpool(number_of_workers);
end
% Propagate warning suppression to all workers
pctRunOnAll warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');

%% Main loop over assigned runs
n_completed_this_chunk = 0;
n_skipped_this_chunk = 0;

for iRun = start_idx:end_idx

    %% Progress tracking -- write current run index so partial failures
    %  can be diagnosed after SLURM walltime/OOM kills
    progress_file = fullfile(outputDir, sprintf('progress_chunk%d.txt', chunk_id));
    fid_prog = fopen(progress_file, 'w');
    fprintf(fid_prog, 'chunk_id: %d\n', chunk_id);
    fprintf(fid_prog, 'current_run: %d\n', iRun);
    fprintf(fid_prog, 'range: %d-%d\n', start_idx, end_idx);
    fprintf(fid_prog, 'completed: %d\n', n_completed_this_chunk);
    fprintf(fid_prog, 'skipped: %d\n', n_skipped_this_chunk);
    fprintf(fid_prog, 'timestamp: %s\n', datestr(now));
    fclose(fid_prog);

    % Extract parameters for this run
    params = combinations(iRun, :);
    SF_flag      = logical(params(1));
    SF_efficiency = params(2);
    junction_depth = params(3);
    sp_r         = params(4);
    N0_peak      = params(5);
    EF_peak      = params(6);
    Phi_R        = params(7);
    Phi_R_offset = params(8);
    profile_idx  = params(9);
    EF_surface   = params(10);
    profile_type = combo_profile_names{iRun};

    % Generate filename (shared function)
    filename = generate_sweep_filename(SF_flag, SF_efficiency, junction_depth, sp_r, ...
                                       N0_peak, Phi_R_offset, profile_type);
    filepath = fullfile(outputDir, filename);

    % Skip if file already exists
    if isfile(filepath)
        fprintf('  [%d/%d] Skipping (exists): %s\n', iRun, total_combos, filename);
        n_skipped_this_chunk = n_skipped_this_chunk + 1;
        continue;
    end

    fprintf('\n=== Run %d/%d (Chunk %d) ===\n', iRun, total_combos, chunk_id);
    fprintf('  SF=%d, eta=%.1f, jd=%d nm, sp=%.0e, N0=%.0e, profile=%s, Phi_R_offset=%.4f\n', ...
            SF_flag, SF_efficiency, junction_depth, sp_r, N0_peak, profile_type, Phi_R_offset);
    fprintf('  -> %s\n', filename);

    try
        %% Generate doping profile for emitter sublayers
        [EF_profile, mu_n_profile, mu_p_profile, taup_profile, N_D_profile, x_sublayers] = ...
            generate_emitter_EF_profile(N0_peak, junction_depth, profile_type, N_sublayers);

        %% Write temp CSV with expanded emitter rows, build par from it
        temp_csv = write_emitter_profile_csv(file_name, EF_profile, mu_n_profile, ...
                                             mu_p_profile, taup_profile, junction_depth);

        par = pc(temp_csv);
        delete(temp_csv);  % Clean up temp CSV after par is built

        par.AbsTol = 5e-6;
        par.RelTol = 5e-3;
        par.MaxStepFactor = 1;
        par.kappa = reflection;

        % Apply fixed parameters
        par.Tetracene_TF = SF_flag;
        par.eta = SF_efficiency;
        par.sp_r = sp_r;
        par.Phi_right = Phi_R;

        par = refresh_device(par);
        par.gx1 = generation(par, par.light_source1, par.laser_lambda1);

        %% Equilibrate (with fallback)
        try
            solEq = equilibrate(par);
        catch ME_eq
            fprintf('  Tight tolerances failed for equilibrate, trying looser...\n');
            par.AbsTol = 1e-5;
            par.RelTol = 1e-2;
            par = refresh_device(par);
            par.gx1 = generation(par, par.light_source1, par.laser_lambda1);
            solEq = equilibrate(par);
        end

        %% JV (with fallback)
        try
            solCV_JV = doCV(solEq.el, light_intensity_JV, ...
                            0, Vmax_JV, Vmin_JV, scan_rate_JV, cycles_JV, tpoints_JV);
        catch ME_cv
            fprintf('  Tight tolerances failed for doCV, trying looser...\n');
            par.AbsTol = 1e-5;
            par.RelTol = 1e-2;
            solEq.el.par.AbsTol = 1e-5;
            solEq.el.par.RelTol = 1e-2;
            solCV_JV = doCV(solEq.el, light_intensity_JV, ...
                            0, Vmax_JV, Vmin_JV, scan_rate_JV, cycles_JV, tpoints_JV);
        end

        xmesh_JV = solCV_JV.x;
        ppos_JV = getpointpos(par.d_midactive, xmesh_JV);

        J_CV_JV = dfana.calcJ(solCV_JV, "sub");
        Vapp_JV = dfana.calcVapp(solCV_JV);

        %% Extract forward sweep (0 -> Vmax segment) robustly
        %  doCV uses tri wave: V0 -> Vmax -> Vmin -> V0 per cycle.
        %  The forward sweep is the first monotonically increasing segment.
        [~, idx_maxV] = max(Vapp_JV);  % Robust: avoids float equality
        V_fwd = Vapp_JV(1:idx_maxV);
        J_fwd = J_CV_JV.tot(1:idx_maxV, ppos_JV);

        % Sanity check: forward sweep voltage should be monotonically increasing
        dV = diff(V_fwd);
        if any(dV < -1e-6)
            warning('Forward JV sweep is not monotonic (max backtrack = %.3e V)', min(dV));
        end

        power_JV = Vapp_JV .* transpose(J_CV_JV.tot(:, ppos_JV));
        MPP_JV = min(power_JV);
        Pin_JV = dfana.calcPin(solCV_JV);
        eff_JV = abs(100 * (MPP_JV / Pin_JV));

        % Voc: interpolate where J crosses zero on the forward sweep
        Voc_JV = interp1(J_fwd, V_fwd, 0, 'linear');
        if isnan(Voc_JV)
            warning('Voc extraction failed: J never crosses zero on forward sweep. J range = [%.3e, %.3e]', ...
                    min(J_fwd), max(J_fwd));
        end

        fprintf('  Voc=%.5fV, PCE=%.5f%%\n', Voc_JV, eff_JV);

        %% Build per-wavelength reflection for EQE
        Wavelength_scan = [scan_range_1, scan_range_2];
        number_of_wavelengths = numel(Wavelength_scan);
        n_bias = length(SR_bias_vals);

        if strcmp(reflection_model, 'spectral')
            % Select cached TMM spectrum for this run's stack
            if SF_flag
                lambda_tmm_raw     = lambda_tmm_SF;
                Reflection_tmm_raw = Reflection_tmm_SF;
            else
                lambda_tmm_raw     = lambda_tmm_noSF;
                Reflection_tmm_raw = Reflection_tmm_noSF;
            end
            kappa_scan = interp1(lambda_tmm_raw, Reflection_tmm_raw, ...
                                 Wavelength_scan, 'linear', reflection);
        else
            % Constant reflection at all wavelengths
            kappa_scan = reflection * ones(1, number_of_wavelengths);
        end

        EQE_array    = zeros(1, number_of_wavelengths);
        IQE_array    = zeros(1, number_of_wavelengths);
        Jsc_array    = zeros(1, number_of_wavelengths);
        SR_array     = zeros(number_of_wavelengths, n_bias);
        SR_EQE_array = zeros(number_of_wavelengths, n_bias);

        % Check pool health; recreate only if it died
        pool = gcp('nocreate');
        if isempty(pool)
            fprintf('  Pool not running, recreating...\n');
            parpool(number_of_workers);
            pctRunOnAll warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
        end

        % parfor with retry on pool crash
        max_attempts = 3;
        for attempt = 1:max_attempts
            try
                parfor n = 1:number_of_wavelengths
                    wavelength = Wavelength_scan(n);
                    kappa_n = kappa_scan(n);
                    [Jsc_array(n), EQE_array(n), IQE_array(n), J_bias_n, EQE_bias_n] = ...
                        EQE_parallel_loop(solEq, wavelength, SF_flag, kappa_n, SR_bias_vals);
                    SR_array(n, :) = J_bias_n;
                    SR_EQE_array(n, :) = EQE_bias_n;
                end
                break  % success
            catch ME_parfor
                fprintf('  parfor attempt %d failed: %s\n', attempt, ME_parfor.message);
                delete(gcp('nocreate'));
                if attempt < max_attempts
                    fprintf('  Restarting pool and retrying...\n');
                    parpool(number_of_workers);
                    pctRunOnAll warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
                else
                    rethrow(ME_parfor);
                end
            end
        end

        %% Build result structure
        thisRun = struct();

        % Parameters
        thisRun.SF_flag        = SF_flag;
        thisRun.SF_efficiency  = SF_efficiency;
        thisRun.junction_depth = junction_depth;
        thisRun.sp_r           = sp_r;
        thisRun.N0_peak        = N0_peak;    % peak doping [cm^-3]
        thisRun.EF_peak        = EF_peak;    % Fermi level at peak [eV]
        thisRun.EF_surface     = EF_surface; % Fermi level at surface [eV]
        thisRun.Phi_R          = Phi_R;
        thisRun.Phi_R_offset   = Phi_R_offset;
        thisRun.profile_type   = profile_type;
        thisRun.N_sublayers    = N_sublayers;

        % Doping profile arrays
        thisRun.profile.x_nm    = x_sublayers;   % sublayer center depths [nm]
        thisRun.profile.N_D     = N_D_profile;   % doping [cm^-3]
        thisRun.profile.EF      = EF_profile;    % Fermi level [eV]
        thisRun.profile.mu_n    = mu_n_profile;
        thisRun.profile.mu_p    = mu_p_profile;
        thisRun.profile.taup    = taup_profile;

        % JV results
        thisRun.Vapp      = Vapp_JV;
        thisRun.Jtot      = J_CV_JV.tot(:, ppos_JV);
        thisRun.Voc       = Voc_JV;
        thisRun.PCE       = eff_JV;
        thisRun.MPP_power = MPP_JV;
        thisRun.Pin       = Pin_JV;

        % EQE
        thisRun.EQE_wavelengths  = Wavelength_scan;
        thisRun.EQE              = EQE_array;
        thisRun.IQE              = IQE_array;
        thisRun.Jsc_vs_lambda    = Jsc_array;
        thisRun.kappa_scan       = kappa_scan;       % reflection at each EQE wavelength
        thisRun.reflection_model = reflection_model;  % 'spectral' | 'constant'

        % Spectral response under bias
        thisRun.SR_bias_vals = SR_bias_vals;
        thisRun.SR           = SR_array;
        thisRun.SR_EQE       = SR_EQE_array;

        % Metadata
        thisRun.timestamp = datetime('now');
        thisRun.par       = par;
        thisRun.solEq     = solEq;
        thisRun.solCV     = solCV_JV;

        %% Save
        save(filepath, 'thisRun', '-v7.3');
        fprintf('  ✓ Saved: %s\n', filename);
        n_completed_this_chunk = n_completed_this_chunk + 1;

    catch ME
        fprintf('  ⚠️ ERROR: %s\n', ME.message);
        for k = 1:length(ME.stack)
            fprintf('    %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
        end
        log_error(outputDir, filename, params, ME, profile_type);
        continue;
    end
end

fprintf('\n===========================================\n');
fprintf('Chunk %d complete!\n', chunk_id);
fprintf('  Runs in chunk: %d-%d (%d total)\n', start_idx, end_idx, end_idx - start_idx + 1);
fprintf('  Completed: %d\n', n_completed_this_chunk);
fprintf('  Skipped (already existed): %d\n', n_skipped_this_chunk);
fprintf('===========================================\n');

% Restore warning state
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames');

end


%% ========================================================================
%  HELPER FUNCTIONS
%% ========================================================================

function log_error(outputDir, filename, params, ME, profile_type)
%LOG_ERROR Write error details to log file

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    errorLog = fullfile(outputDir, 'errors.txt');
    fid = fopen(errorLog, 'a');
    fprintf(fid, '\n=== Error at %s ===\n', datestr(now));
    fprintf(fid, 'File: %s\n', filename);
    fprintf(fid, 'Params: SF=%d, eta=%.1f, jd=%d, sp=%.0e, N0=%.0e, PhiOff=%.4f, prof_idx=%d, profile=%s\n', ...
            params(1), params(2), params(3), params(4), params(5), params(8), params(9), profile_type);
    fprintf(fid, 'Error: %s\n', ME.message);
    for k = 1:length(ME.stack)
        fprintf(fid, '  %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
    end
    fclose(fid);
end


function [Jsc, EQE_calc, IQE_calc, J_at_bias, EQE_at_bias] = ...
        EQE_parallel_loop(soleq, wavelength, Tetracene_TF, kappa_lambda, SR_bias_vals)
%EQE_PARALLEL_LOOP Calculate EQE and spectral response under bias at a single wavelength
%
%   kappa_lambda: wavelength-specific reflectance from TMM (scalar, dimensionless).
%                 Interpolated from the TMM spectrum at the caller level.

    h = 6.62607015e-34;
    e = 1.602176634e-19;
    c = 2.99792458e8;

    par_temp = soleq.el.par;
    par_temp.Tetracene_TF = Tetracene_TF;
    par_temp.kappa = kappa_lambda;

    par_temp.light_source1 = 'laser';
    par_temp.laser_lambda1 = wavelength;
    par_temp.pulsepow = 1;

    par_temp.gx1 = generation(par_temp, 'laser', wavelength);

    soleq_temp = soleq.el;
    soleq_temp.par = par_temp;

    %% Narrowed voltage range -- we only need J at V=0 (and SR_bias_vals if any).
    %  With empty SR_bias_vals, min([[], -0.3]) = -0.3.
    Vmin = min([min(SR_bias_vals), -0.3]);
    Vmax = 0.1;            % Only need to interpolate at V <= 0
    scan_rate = 50e-3;
    cycles = 1;
    tpoints = 200;          % Halved from 400: narrower range needs fewer points

    try
        CVsol = doCV(soleq_temp, 1, 0, Vmax, Vmin, scan_rate, cycles, tpoints);
    catch
        soleq_temp.par.AbsTol = 1e-5;
        soleq_temp.par.RelTol = 1e-2;
        CVsol = doCV(soleq_temp, 1, 0, Vmax, Vmin, scan_rate, cycles, tpoints);
    end

    V = dfana.calcVapp(CVsol);
    J = dfana.calcJ(CVsol, "sub").tot;

    %% Extract forward sweep (0 -> Vmax) robustly
    %  doCV tri wave: V0(=0) -> Vmax -> Vmin -> V0. First segment is forward.
    [~, idx_maxV] = max(V);
    J_f = J(1:idx_maxV);
    V_f = V(1:idx_maxV);

    %% Jsc at V = 0 with validation
    if min(V_f) > 0 || max(V_f) < 0
        warning('EQE_parallel_loop: V=0 outside forward sweep range [%.4f, %.4f] at lambda=%d nm', ...
                min(V_f), max(V_f), wavelength);
        Jsc = NaN;
    else
        Jsc = abs(interp1(V_f, J_f, 0, 'linear'));
    end

    Pin = par_temp.pulsepow * 1e-3;

    %% EQE
    EQE_calc = 100 * Jsc * h * c / (Pin * wavelength * 1e-9 * e);

    %% IQE: correct for wavelength-dependent front-surface reflection
    %  IQE = EQE / (1 - R(lambda)), where R(lambda) = kappa_lambda from TMM.
    %
    %  The old code used beerlambertI with fragile wavelength-300 indexing.
    %  Parasitic absorption in Tc/ZnPc is already embedded in the generation
    %  profile via the transfer matrix model, so kappa is the only purely
    %  optical loss not captured by the electrical simulation.
    if kappa_lambda < 1
        IQE_calc = EQE_calc / (1 - kappa_lambda);
    else
        IQE_calc = NaN;  % total reflection -> IQE undefined
    end

    %% Spectral response under bias
    n_bias = length(SR_bias_vals);
    J_at_bias   = zeros(1, n_bias);
    EQE_at_bias = zeros(1, n_bias);

    for ib = 1:n_bias
        bias = SR_bias_vals(ib);
        if bias >= min(V_f) && bias <= max(V_f)
            J_at_bias(ib) = abs(interp1(V_f, J_f, bias, 'linear'));
        else
            J_at_bias(ib) = NaN;
        end
        EQE_at_bias(ib) = 100 * J_at_bias(ib) * h * c / (Pin * wavelength * 1e-9 * e);
    end
end