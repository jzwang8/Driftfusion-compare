function run_sweep_chunk(chunk_id)
%RUN_SWEEP_CHUNK Run a specific chunk of the parameter sweep
%   chunk_id: which chunk to run (0-based from SLURM array)
%
%   This version:
%   - Sweeps Phi_R relative to each N_D's Fermi level
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
% ORIGINAL SWEEP VALUES
% junction_depth_vals = [200, 250, 300, 350, 400];      % nm
% sp_r_vals = [10, 100, 1000, 10000, 100000];             % cm/s
% eta_sf_with_SF = [0, 0.5, 1, 1.5, 2];                 % SF efficiency when SF_flag=1
%                                             % Example: add 1.5 later -> [0, 1, 1.5, 2]
% 
% % Sweep EF directly (easier since Phi_R is relative to EF)
% % N_D will be calculated from EF for mobility/lifetime
% EF_vals = [-4.05, -4.104, -4.223, -4.342];  % eV (corresponds to ~8e19, 1e19, 1e17, 1e15 cm^-3)
% 
% % Phi_R offset for photoinduced passivation (relative to each EF)
% % Physical meaning: positive offset = positive Q_f = field-effect passivation
% Phi_R_offset_vals = [0.0539, 0, -0.0539];   % eV (positive, neutral, negative Q_f)
% 
% % Physical constraint: Phi_R cannot exceed conduction band edge
% E_C = -4.05;  % Si conduction band edge (eV)

% NEW SWEEP VALUES (REVERSE BIAS)
junction_depth_vals = [200];      % nm
sp_r_vals = [1000, 10000, 100000];             % cm/s
eta_sf_with_SF = [0, 1, 2];                 % SF efficiency when SF_flag=1

% Sweep EF directly (easier since Phi_R is relative to EF)
% N_D will be calculated from EF for mobility/lifetime
EF_vals = [-4.05, -4.104, -4.223, -4.342];  % eV (corresponds to ~8e19, 1e19, 1e17, 1e15 cm^-3)

% Phi_R offset for photoinduced passivation (relative to each EF)
% Physical meaning: positive offset = positive Q_f = field-effect passivation
Phi_R_offset_vals = [0.0540, 0,-0.0540];   % eV (positive, neutral, negative Q_f)

% Physical constraint: Phi_R cannot exceed conduction band edge
E_C = -4.05;  % Si conduction band edge (eV)

% -------------------------------------------------------------------------
% Spectral response under bias settings
% Bias values at which to extract J(lambda) in addition to EQE (V=0)
% These are computed within each run - not separate sweep dimensions
% Must include 0 if you want to cross-check against standard EQE
SR_bias_vals = [-1.0, -0.5, 0.0];   % V (applied external bias)

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
outputDir = 'sweep_results';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Count existing results for info
existingFiles = dir(fullfile(outputDir, 'SF*.mat'));
fprintf('Found %d existing result files\n', length(existingFiles));

%% Pre-compute N_D for each EF value (for mobility/lifetime calculations)
EF_to_ND = containers.Map('KeyType', 'double', 'ValueType', 'double');

fprintf('\nFermi level -> Doping mapping:\n');
for i = 1:length(EF_vals)
    EF = EF_vals(i);
    [N_D, ~] = getDopingConcFromEf(EF);
    EF_to_ND(EF) = N_D;
    fprintf('  E_F = %.4f eV -> N_D = %.2e cm^-3\n', EF, N_D);
end

%% Generate all parameter combinations intelligently
% SF_flag=0: only eta_sf=0 (no SF contribution, so eta doesn't matter)
% SF_flag=1: all eta_sf values (meaningful when SF is on)
% Phi_R is swept relative to each EF

combinations = [];
n_skipped_invalid = 0;

% --- SF_flag = 0, eta_sf = 0 only (baseline Si cell, no SF) ---
SF_flag = 0;
eta_sf = 0;

for i_jd = 1:length(junction_depth_vals)
    for i_sp = 1:length(sp_r_vals)
        for i_EF = 1:length(EF_vals)
            for i_offset = 1:length(Phi_R_offset_vals)
                
                junction_depth = junction_depth_vals(i_jd);
                sp_r = sp_r_vals(i_sp);
                EF_emitter = EF_vals(i_EF);
                N_D = EF_to_ND(EF_emitter);
                Phi_R_offset = Phi_R_offset_vals(i_offset);
                Phi_R = EF_emitter + Phi_R_offset;
                
                % Skip if Phi_R exceeds conduction band edge (unphysical)
                if Phi_R > E_C
                    n_skipped_invalid = n_skipped_invalid + 1;
                    continue;
                end
                
                combinations = [combinations; ...
                    SF_flag, eta_sf, junction_depth, sp_r, N_D, EF_emitter, Phi_R, Phi_R_offset]; %#ok<AGROW>
            end
        end
    end
end

n_baseline = size(combinations, 1);

% --- SF_flag = 1, all eta_sf values (SF layer present) ---
SF_flag = 1;

for i_eta = 1:length(eta_sf_with_SF)
    for i_jd = 1:length(junction_depth_vals)
        for i_sp = 1:length(sp_r_vals)
            for i_EF = 1:length(EF_vals)
                for i_offset = 1:length(Phi_R_offset_vals)
                    
                    eta_sf = eta_sf_with_SF(i_eta);
                    junction_depth = junction_depth_vals(i_jd);
                    sp_r = sp_r_vals(i_sp);
                    EF_emitter = EF_vals(i_EF);
                    N_D = EF_to_ND(EF_emitter);
                    Phi_R_offset = Phi_R_offset_vals(i_offset);
                    Phi_R = EF_emitter + Phi_R_offset;
                    
                    % Skip if Phi_R exceeds conduction band edge (unphysical)
                    if Phi_R > E_C
                        n_skipped_invalid = n_skipped_invalid + 1;
                        continue;
                    end
                    
                    combinations = [combinations; ...
                        SF_flag, eta_sf, junction_depth, sp_r, N_D, EF_emitter, Phi_R, Phi_R_offset]; %#ok<AGROW>
                end
            end
        end
    end
end

n_with_SF = size(combinations, 1) - n_baseline;
total_combos = size(combinations, 1)

fprintf('\nParameter space:\n');
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
    SF_flag = logical(params(1));
    SF_efficiency = params(2);
    junction_depth = params(3);
    sp_r = params(4);
    N_D = params(5);
    EF_emitter = params(6);
    Phi_R = params(7);
    Phi_R_offset = params(8);
    
    % Generate filename based on parameters (deterministic, unique)
    filename = generate_filename(SF_flag, SF_efficiency, junction_depth, sp_r, N_D, Phi_R_offset);
    filepath = fullfile(outputDir, filename);
    
    % Skip if file already exists
    if isfile(filepath)
        fprintf('  [%d/%d] Skipping (exists): %s\n', iRun, total_combos, filename);
        n_skipped_this_chunk = n_skipped_this_chunk + 1;
        continue;
    end
    
    fprintf('\n=== Run %d/%d (Chunk %d) ===\n', iRun, total_combos, chunk_id);
    fprintf('  SF=%d, eta=%.1f, jd=%d nm, sp=%.0e, N_D=%.0e, Phi_R_offset=%.4f\n', ...
            SF_flag, SF_efficiency, junction_depth, sp_r, N_D, Phi_R_offset);
    fprintf('  -> %s\n', filename);
    
    try
        % Build parameter object with TIGHT tolerances first
        par = pc(file_name);
        par.AbsTol = 5e-6;
        par.RelTol = 5e-3;
        par.MaxStepFactor = 1;
        par.kappa = reflection;
        
        % Apply parameters
        par.Tetracene_TF = SF_flag;
        par.eta = SF_efficiency;
        par.d(1,3) = junction_depth * 1e-9 * 1e2;  % convert nm to cm
        par.sp_r = sp_r;
        par.EF0(3) = EF_emitter;
        par.Phi_right = Phi_R;
        
        % Update mobilities and lifetime based on doping
        [mu_n, mu_p] = getMobilitiesFromDopingConc(N_D);
        par.mu_n(3) = mu_n;
        par.mu_p(3) = mu_p;
        tau_p = getTaupFromDopingConc(N_D);
        par.taup(3) = tau_p;
        
        par = refresh_device(par);
        par.gx1 = generation(par, par.light_source1, par.laser_lambda1);
        
        %% Run equilibrium and JV (with fallback for tolerance)
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
        ppos_JV = getpointpos(par.d_midactive, xmesh_JV); % TODO: check mid or endpoint
        
        J_CV_JV = dfana.calcJ(solCV_JV, "sub");
        Vapp_JV = dfana.calcVapp(solCV_JV);
        
        power_JV = Vapp_JV .* transpose(J_CV_JV.tot(:, ppos_JV));
        MPP_JV = min(power_JV);
        Pin_JV = dfana.calcPin(solCV_JV);
        eff_JV = abs(100 * (MPP_JV / Pin_JV));
        Voc_JV = interp1(J_CV_JV.tot(:, ppos_JV), Vapp_JV(:), 0);
        
        fprintf('  Voc=%.5fV, PCE=%.5f%%\n', Voc_JV, eff_JV);
        
        %% Run EQE + Spectral Response under Bias
        Wavelength_scan = [scan_range_1, scan_range_2];
        number_of_wavelengths = numel(Wavelength_scan);
        n_bias = length(SR_bias_vals);
        
        % Standard EQE arrays (J extracted at V=0)
        EQE_array = zeros(1, number_of_wavelengths);
        IQE_array = zeros(1, number_of_wavelengths);
        Jsc_array = zeros(1, number_of_wavelengths);
        
        % Spectral response under bias: rows = wavelengths, cols = bias points
        SR_array = zeros(number_of_wavelengths, n_bias);   % J(lambda, bias) in A/m2
        SR_EQE_array = zeros(number_of_wavelengths, n_bias); % EQE-equivalent at each bias
        
        if isempty(gcp('nocreate'))
            parpool(number_of_workers);
        end
        
        % Reset pool if a previous run failed
        try
            % Quick health check
            parfor testIdx = 1:number_of_workers
                % no-op
            end
        catch
            delete(gcp('nocreate'));
            parpool(number_of_workers);
        end

        parfor n = 1:number_of_wavelengths
            wavelength = Wavelength_scan(n);
            [Jsc_array(n), EQE_array(n), IQE_array(n), J_bias_n, EQE_bias_n] = ...
                EQE_parallel_loop(solEq, wavelength, SF_flag, reflection, SR_bias_vals);
            SR_array(n, :) = J_bias_n;
            SR_EQE_array(n, :) = EQE_bias_n;
        end
        
        %% Build result structure
        thisRun = struct();
        
        % Parameters
        thisRun.SF_flag = SF_flag;
        thisRun.SF_efficiency = SF_efficiency;
        thisRun.junction_depth = junction_depth;
        thisRun.sp_r = sp_r;
        thisRun.N_D = N_D;
        thisRun.EF_emitter = EF_emitter;
        thisRun.Phi_R = Phi_R;
        thisRun.Phi_R_offset = Phi_R_offset;
        
        % Standard JV results
        thisRun.Vapp = Vapp_JV;
        thisRun.Jtot = J_CV_JV.tot(:, ppos_JV);
        thisRun.Voc = Voc_JV;
        thisRun.PCE = eff_JV;
        thisRun.MPP_power = MPP_JV;
        thisRun.Pin = Pin_JV;
        
        % Standard EQE (J extracted at V=0)
        thisRun.EQE_wavelengths = Wavelength_scan;
        thisRun.EQE = EQE_array;
        thisRun.IQE = IQE_array;
        thisRun.Jsc_vs_lambda = Jsc_array;
        
        % Spectral response under bias
        % SR_array:     [n_wavelengths x n_bias]  - J in A/m2 at each bias
        % SR_EQE_array: [n_wavelengths x n_bias]  - EQE-equivalent at each bias
        % SR_bias_vals: [1 x n_bias]              - applied bias values in V
        thisRun.SR_bias_vals = SR_bias_vals;
        thisRun.SR = SR_array;
        thisRun.SR_EQE = SR_EQE_array;
        
        % Metadata
        thisRun.timestamp = datetime('now');
        thisRun.par = par;
        thisRun.solEq = solEq;
        thisRun.solCV = solCV_JV;
        
        %% Save
        save(filepath, 'thisRun', '-v7.3');
        fprintf('  ✓ Saved: %s\n', filename);
        n_completed_this_chunk = n_completed_this_chunk + 1;
        
    catch ME
        fprintf('  ⚠️ ERROR: %s\n', ME.message);
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
%  ========================================================================

function filename = generate_filename(SF_flag, SF_efficiency, junction_depth, sp_r, N_D, Phi_R_offset)
%GENERATE_FILENAME Create deterministic filename from parameters
%   Format: SF{0/1}_eta{X}_jd{X}_sp{X}_ND{X}_PhiOff{X}.mat

    sf_str = sprintf('SF%d', SF_flag);
    eta_str = sprintf('eta%.1f', SF_efficiency);
    jd_str = sprintf('jd%d', junction_depth);
    sp_str = sprintf('sp%.0e', sp_r);
    nd_str = sprintf('ND%.0e', N_D);
    phi_str = sprintf('PhiOff%.4f', Phi_R_offset);
    
    sp_str = strrep(sp_str, '+', '');
    nd_str = strrep(nd_str, '+', '');
    phi_str = strrep(phi_str, '-', 'n');
    phi_str = strrep(phi_str, '.', 'p');
    
    filename = sprintf('%s_%s_%s_%s_%s_%s.mat', sf_str, eta_str, jd_str, sp_str, nd_str, phi_str);
end


function log_error(outputDir, filename, params, ME)
%LOG_ERROR Write error details to log file

    errorLog = fullfile(outputDir, 'errors.txt');
    fid = fopen(errorLog, 'a');
    fprintf(fid, '\n=== Error at %s ===\n', datestr(now));
    fprintf(fid, 'File: %s\n', filename);
    fprintf(fid, 'Params: SF=%d, eta=%.1f, jd=%d, sp=%.0e, N_D=%.0e, Phi_R_off=%.4f\n', ...
            params(1), params(2), params(3), params(4), params(5), params(8));
    fprintf(fid, 'Error: %s\n', ME.message);
    for k = 1:length(ME.stack)
        fprintf(fid, '  %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
    end
    fclose(fid);
end


function [Jsc, EQE_calc, IQE_calc, J_at_bias, EQE_at_bias] = ...
        EQE_parallel_loop(soleq, wavelength, Tetracene_TF, reflection, SR_bias_vals)
%EQE_PARALLEL_LOOP Calculate EQE and spectral response under bias at a single wavelength
%
%   Outputs:
%     Jsc        - J at V=0 (standard short-circuit current)
%     EQE_calc   - EQE at V=0 (%)
%     IQE_calc   - IQE at V=0 (%)
%     J_at_bias  - [1 x n_bias] J extracted at each bias in SR_bias_vals (A/m2)
%     EQE_at_bias- [1 x n_bias] EQE-equivalent at each bias (%)

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

    % Vmin must cover the most negative bias point
    Vmin = min([min(SR_bias_vals), -0.3]) - 0.1;  % extra 0.1V margin
    Vmax = 0.6;
    scan_rate = 50e-3;
    cycles = 1;
    tpoints = 400;
    
    % Try with current tolerances, fallback to looser if needed
    try
        CVsol = doCV(soleq_temp, 1, 0, Vmax, Vmin, scan_rate, cycles, tpoints);
    catch
        soleq_temp.par.AbsTol = 1e-5;
        soleq_temp.par.RelTol = 1e-2;
        CVsol = doCV(soleq_temp, 1, 0, Vmax, Vmin, scan_rate, cycles, tpoints);
    end

    V = dfana.calcVapp(CVsol);
    V = V(:);
    J_all = dfana.calcJ(CVsol, "sub");
    J = J_all.tot(:, end);
    J = J(:);
    
    % Forward sweep: 0 -> Vmax
    idx_maxV = find(V == max(V), 1);
    J_f = J(1:idx_maxV);
    V_f = V(1:idx_maxV);

    % Return sweep: Vmax -> Vmin (covers negative voltages)
    J_r = J(idx_maxV:end);
    V_r = V(idx_maxV:end);

    % Standard EQE: J at V=0 (from forward sweep)
    Jsc = abs(interp1(V_f, J_f, 0, 'linear'));
    Pin = par_temp.pulsepow * 1e-3;

    I = beerlambertI(par_temp, par_temp.xx, par_temp.light_source1, wavelength, 0);
    Plost = I(1, wavelength-300);

    EQE_calc = 100 * Jsc * h * c / (Pin * wavelength * 1e-9 * e);
    IQE_calc = EQE_calc * Pin / (Pin - Plost);

    % Spectral response under bias: J at each bias point
    n_bias = length(SR_bias_vals);
    J_at_bias = zeros(1, n_bias);
    EQE_at_bias = zeros(1, n_bias);

    for ib = 1:n_bias
        bias = SR_bias_vals(ib);
        if bias >= 0 && bias <= max(V_f)
            J_at_bias(ib) = abs(interp1(V_f, J_f, bias, 'linear'));
        elseif bias < 0 && bias >= min(V_r) && bias <= max(V_r)
            J_at_bias(ib) = abs(interp1(V_r, J_r, bias, 'linear'));
        else
            J_at_bias(ib) = NaN;
        end
        EQE_at_bias(ib) = 100 * J_at_bias(ib) * h * c / (Pin * wavelength * 1e-9 * e);
    end
end