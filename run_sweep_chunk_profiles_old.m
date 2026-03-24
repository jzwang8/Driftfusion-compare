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
%   - Computes EQE (J at V=0) AND spectral response under bias (J at SR_bias_vals)

%% Setup
delete(gcp('nocreate'));

% Add necessary paths
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Transfer matrix modeling'));
addpath(fullfile(thisDir, 'Doping'));

file_name = 'Input_files/pn_junction_nochargedlayer.csv';

%% Define sweep parameters
% -------------------------------------------------------------------------
% MODIFY THESE TO ADD NEW VALUES - existing runs won't be repeated
% -------------------------------------------------------------------------

junction_depth_vals = [200];      % nm
sp_r_vals = [1000, 10000, 100000];  % cm/s
eta_sf_with_SF = [0, 1, 2];        % SF efficiency when SF_flag=1

% EF_vals now represents PEAK doping concentration
% (surface concentration for exponential; center peak for Gaussian)
EF_vals = [-4.05, -4.104, -4.223, -4.342];  % eV (~8e19, 1e19, 1e17, 1e15 cm^-3)

% Phi_R offset relative to the SURFACE EF (EF of outermost sublayer)
Phi_R_offset_vals = [0.054, 0, -0.054];   % eV

% Physical constraint: Phi_R cannot exceed conduction band edge
E_C = -4.05;  % Si conduction band edge (eV)

% Doping profile shapes (new sweep dimension)
% 'uniform'     - constant doping (matches previous sweep behavior)
% 'gaussian'    - ion implantation profile, peak at x0=d/2, sigma=d/5
% 'exponential' - diffusion doping, dose-matched to Gaussian
profile_shapes = {'uniform', 'gaussian', 'exponential'};

% Number of emitter sublayers for profile discretization
% (uniform uses 1 layer; gaussian/exponential use N_sublayers)
N_sublayers = 50;

% -------------------------------------------------------------------------
% Spectral response under bias settings
SR_bias_vals = [0.0];   % V (applied external bias)

% -------------------------------------------------------------------------
% Fixed parameters
reflection = 0.15;

% JV settings
light_intensity_JV = 1;
Vmin_JV = -0.3;
Vmax_JV = 0.7;
scan_rate_JV = 50e-3;
cycles_JV = 2;
tpoints_JV = 300;

% EQE settings
scan_range_1 = 405:5:600;
scan_range_2 = 610:10:800;
number_of_workers = 8;

%% Setup output directory
outputDir = 'sweep_results_profiles';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

existingFiles = dir(fullfile(outputDir, 'SF*.mat'));
fprintf('Found %d existing result files\n', length(existingFiles));

%% Pre-compute N_D for each EF value
EF_to_ND = containers.Map('KeyType', 'double', 'ValueType', 'double');

fprintf('\nFermi level (peak) -> Doping mapping:\n');
for i = 1:length(EF_vals)
    EF = EF_vals(i);
    [N_D, ~] = getDopingConcFromEf(EF);
    EF_to_ND(EF) = N_D;
    fprintf('  E_F_peak = %.4f eV -> N0_peak = %.2e cm^-3\n', EF, N_D);
end

%% Generate all parameter combinations
% combinations columns:
%   1: SF_flag
%   2: eta_sf
%   3: junction_depth (nm)
%   4: sp_r (cm/s)
%   5: N0_peak (cm^-3)  -- peak doping
%   6: EF_peak (eV)     -- Fermi level at peak doping
%   7: Phi_R (eV)       -- = EF_surface + Phi_R_offset
%   8: Phi_R_offset (eV)
%   9: profile_idx      -- index into profile_shapes cell array
%  10: EF_surface (eV)  -- Fermi level at outermost sublayer (used for Phi_R ref)

combinations = [];
combo_profile_names = {};  % profile name strings (parallel to combinations rows)
n_skipped_invalid = 0;

% --- SF_flag = 0, eta_sf = 0 only ---
for i_prof = 1:length(profile_shapes)
    profile_type = profile_shapes{i_prof};
    for i_jd = 1:length(junction_depth_vals)
        for i_sp = 1:length(sp_r_vals)
            for i_EF = 1:length(EF_vals)
                for i_offset = 1:length(Phi_R_offset_vals)

                    junction_depth = junction_depth_vals(i_jd);
                    sp_r = sp_r_vals(i_sp);
                    EF_peak = EF_vals(i_EF);
                    N0_peak = EF_to_ND(EF_peak);
                    Phi_R_offset = Phi_R_offset_vals(i_offset);

                    % Get surface EF for this profile/doping combo
                    EF_surface = get_surface_EF(N0_peak, junction_depth, profile_type, N_sublayers);
                    Phi_R = EF_surface + Phi_R_offset;

                    if Phi_R > E_C
                        n_skipped_invalid = n_skipped_invalid + 1;
                        continue;
                    end

                    combinations = [combinations; ...
                        0, 0, junction_depth, sp_r, N0_peak, EF_peak, Phi_R, Phi_R_offset, i_prof, EF_surface]; %#ok<AGROW>
                    combo_profile_names{end+1} = profile_type; %#ok<AGROW>
                end
            end
        end
    end
end

n_baseline = size(combinations, 1);

% --- SF_flag = 1, all eta_sf values ---
for i_eta = 1:length(eta_sf_with_SF)
    for i_prof = 1:length(profile_shapes)
        profile_type = profile_shapes{i_prof};
        for i_jd = 1:length(junction_depth_vals)
            for i_sp = 1:length(sp_r_vals)
                for i_EF = 1:length(EF_vals)
                    for i_offset = 1:length(Phi_R_offset_vals)

                        eta_sf = eta_sf_with_SF(i_eta);
                        junction_depth = junction_depth_vals(i_jd);
                        sp_r = sp_r_vals(i_sp);
                        EF_peak = EF_vals(i_EF);
                        N0_peak = EF_to_ND(EF_peak);
                        Phi_R_offset = Phi_R_offset_vals(i_offset);

                        EF_surface = get_surface_EF(N0_peak, junction_depth, profile_type, N_sublayers);
                        Phi_R = EF_surface + Phi_R_offset;

                        if Phi_R > E_C
                            n_skipped_invalid = n_skipped_invalid + 1;
                            continue;
                        end

                        combinations = [combinations; ...
                            1, eta_sf, junction_depth, sp_r, N0_peak, EF_peak, Phi_R, Phi_R_offset, i_prof, EF_surface]; %#ok<AGROW>
                        combo_profile_names{end+1} = profile_type; %#ok<AGROW>
                    end
                end
            end
        end
    end
end

n_with_SF = size(combinations, 1) - n_baseline;
total_combos = size(combinations, 1)

fprintf('\nParameter space:\n');
fprintf('  Profile shapes: %s\n', strjoin(profile_shapes, ', '));
fprintf('  SF_flag=0 (baseline): %d combinations\n', n_baseline);
fprintf('  SF_flag=1 (with SF): %d combinations\n', n_with_SF);
fprintf('  Total valid: %d combinations\n', total_combos);
fprintf('  Skipped (Phi_R > E_C): %d invalid combinations\n', n_skipped_invalid);
fprintf('  Spectral response bias points: %s V\n', mat2str(SR_bias_vals));

%% Determine which runs this chunk should do
runs_per_chunk = 20;
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

%% Main loop over assigned runs
n_completed_this_chunk = 0;
n_skipped_this_chunk = 0;

for iRun = start_idx:end_idx

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

    % Generate filename
    filename = generate_filename(SF_flag, SF_efficiency, junction_depth, sp_r, ...
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

        % DEBUG: print first 3 lines of temp CSV to verify header
        fdbg = fopen(temp_csv, 'r');
        fprintf('  [debug] First 3 lines of temp CSV:\n');
        for dbgi = 1:3
            ln = fgetl(fdbg);
            if ischar(ln)
                fprintf('    LINE %d: %s\n', dbgi, ln(1:min(80,end)));
            end
        end
        fclose(fdbg);

        par = pc(temp_csv);
        % delete(temp_csv);  % DEBUG: keep temp file for inspection

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

        power_JV = Vapp_JV .* transpose(J_CV_JV.tot(:, ppos_JV));
        MPP_JV = min(power_JV);
        Pin_JV = dfana.calcPin(solCV_JV);
        eff_JV = abs(100 * (MPP_JV / Pin_JV));
        Voc_JV = interp1(J_CV_JV.tot(:, ppos_JV), Vapp_JV, 0);

        fprintf('  Voc=%.5fV, PCE=%.5f%%\n', Voc_JV, eff_JV);

        %% EQE + Spectral Response under Bias
        Wavelength_scan = [scan_range_1, scan_range_2];
        number_of_wavelengths = numel(Wavelength_scan);
        n_bias = length(SR_bias_vals);

        EQE_array    = zeros(1, number_of_wavelengths);
        IQE_array    = zeros(1, number_of_wavelengths);
        Jsc_array    = zeros(1, number_of_wavelengths);
        SR_array     = zeros(number_of_wavelengths, n_bias);
        SR_EQE_array = zeros(number_of_wavelengths, n_bias);

        % Start/restart pool
        try
            pool = gcp('nocreate');
            if isempty(pool)
                parpool(number_of_workers);
            end
        catch
            delete(gcp('nocreate'));
            parpool(number_of_workers);
        end

        % parfor with retry on pool crash
        max_attempts = 3;
        for attempt = 1:max_attempts
            try
                parfor n = 1:number_of_wavelengths
                    wavelength = Wavelength_scan(n);
                    [Jsc_array(n), EQE_array(n), IQE_array(n), J_bias_n, EQE_bias_n] = ...
                        EQE_parallel_loop(solEq, wavelength, SF_flag, reflection, SR_bias_vals);
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
        thisRun.EQE_wavelengths = Wavelength_scan;
        thisRun.EQE             = EQE_array;
        thisRun.IQE             = IQE_array;
        thisRun.Jsc_vs_lambda   = Jsc_array;

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
        log_error(outputDir, filename, params, ME);
        continue;
    end
end

fprintf('\n===========================================\n');
fprintf('Chunk %d complete!\n', chunk_id);
fprintf('  Runs in chunk: %d-%d (%d total)\n', start_idx, end_idx, end_idx - start_idx + 1);
fprintf('  Completed: %d\n', n_completed_this_chunk);
fprintf('  Skipped (already existed): %d\n', n_skipped_this_chunk);
fprintf('===========================================\n');

end


%% ========================================================================
%  HELPER FUNCTIONS
%% ========================================================================

function EF_surface = get_surface_EF(N0_peak, junction_depth_nm, profile_type, N_sublayers)
%GET_SURFACE_EF Return EF at the outermost (surface) sublayer for Phi_R reference
%   Fast — does not compute full profile, just the surface value

    d = junction_depth_nm;
    sublayer_thickness = d / N_sublayers;
    x_surface = sublayer_thickness / 2;  % center of first sublayer (nm)

    switch lower(profile_type)
        case 'uniform'
            N_surface = N0_peak;
        case 'gaussian'
            x0    = d / 2;
            sigma = d / 5;
            N_surface = N0_peak * exp(-((x_surface - x0)^2) / (2 * sigma^2));
            N_surface = max(N_surface, 1e15);
        case 'exponential'
            % Surface = A = N0_peak (exponential starts at peak at surface)
            N_surface = N0_peak;
        otherwise
            N_surface = N0_peak;
    end

    EF_surface = getEfFromDopingConc(N_surface);
end




function filename = generate_filename(SF_flag, SF_efficiency, junction_depth, sp_r, N0_peak, Phi_R_offset, profile_type)
%GENERATE_FILENAME Create deterministic filename from parameters
%   Format: SF{0/1}_eta{X}_jd{X}_sp{X}_N0{X}_PhiOff{X}_prof{X}.mat

    sf_str   = sprintf('SF%d', SF_flag);
    eta_str  = sprintf('eta%.1f', SF_efficiency);
    jd_str   = sprintf('jd%d', junction_depth);
    sp_str   = sprintf('sp%.0e', sp_r);
    n0_str   = sprintf('N0%.0e', N0_peak);
    phi_str  = sprintf('PhiOff%.4f', Phi_R_offset);
    prof_str = sprintf('prof%s', profile_type);

    sp_str  = strrep(sp_str,  '+', '');
    n0_str  = strrep(n0_str,  '+', '');
    phi_str = strrep(phi_str, '-', 'n');
    phi_str = strrep(phi_str, '.', 'p');

    filename = sprintf('%s_%s_%s_%s_%s_%s_%s.mat', ...
                       sf_str, eta_str, jd_str, sp_str, n0_str, phi_str, prof_str);
end


function log_error(outputDir, filename, params, ME)
%LOG_ERROR Write error details to log file

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    errorLog = fullfile(outputDir, 'errors.txt');
    fid = fopen(errorLog, 'a');
    fprintf(fid, '\n=== Error at %s ===\n', datestr(now));
    fprintf(fid, 'File: %s\n', filename);
    fprintf(fid, 'Params: SF=%d, eta=%.1f, jd=%d, sp=%.0e, N0=%.0e, PhiOff=%.4f, prof_idx=%d\n', ...
            params(1), params(2), params(3), params(4), params(5), params(8), params(9));
    fprintf(fid, 'Error: %s\n', ME.message);
    for k = 1:length(ME.stack)
        fprintf(fid, '  %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
    end
    fclose(fid);
end


function [Jsc, EQE_calc, IQE_calc, J_at_bias, EQE_at_bias] = ...
        EQE_parallel_loop(soleq, wavelength, Tetracene_TF, reflection, SR_bias_vals)
%EQE_PARALLEL_LOOP Calculate EQE and spectral response under bias at a single wavelength

    h = 6.62607015e-34;
    e = 1.602176634e-19;
    c = 2.99792458e8;

    par_temp = soleq.el.par;
    par_temp.Tetracene_TF = Tetracene_TF;
    par_temp.kappa = reflection;

    par_temp.light_source1 = 'laser';
    par_temp.laser_lambda1 = wavelength;
    par_temp.pulsepow = 1;

    par_temp.gx1 = generation(par_temp, 'laser', wavelength);

    soleq_temp = soleq.el;
    soleq_temp.par = par_temp;

    Vmin = min([min(SR_bias_vals), -0.3]);
    Vmax = 0.6;
    scan_rate = 50e-3;
    cycles = 1;
    tpoints = 400;

    try
        CVsol = doCV(soleq_temp, 1, 0, Vmax, Vmin, scan_rate, cycles, tpoints);
    catch
        soleq_temp.par.AbsTol = 1e-5;
        soleq_temp.par.RelTol = 1e-2;
        CVsol = doCV(soleq_temp, 1, 0, Vmax, Vmin, scan_rate, cycles, tpoints);
    end

    V = dfana.calcVapp(CVsol);
    J_all = dfana.calcJ(CVsol, "sub");
    J = J_all.tot(:, end);

    idx_maxV = find(V == max(V), 1);
    J_f = J(1:idx_maxV);
    V_f = V(1:idx_maxV);

    Jsc = abs(interp1(V_f, J_f, 0, 'linear'));
    Pin = par_temp.pulsepow * 1e-3;

    I = beerlambertI(par_temp, par_temp.xx, par_temp.light_source1, wavelength, 0);
    Plost = I(1, wavelength-300);

    EQE_calc = 100 * Jsc * h * c / (Pin * wavelength * 1e-9 * e);
    IQE_calc = EQE_calc * Pin / (Pin - Plost);

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