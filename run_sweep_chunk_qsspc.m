function run_sweep_chunk_qsspc(chunk_id)
%RUN_SWEEP_CHUNK_QSSPC Run a specific chunk of the QSSPC lifetime sweep
%   chunk_id: which chunk to run (0-based from SLURM array)
%
%   This version:
%   - Sweeps ND (wafer doping), sp_r (chemical passivation), Vi (Qf offset)
%   - CSV is a template; per-run overrides recalculate EF0, mu, tau from ND
%   - Skips runs where identical parameters already exist (not by index)
%   - Safe to add new parameter values without re-running old combinations
%   - Uses calcQfFromVi() for Vi<->Qf conversion
%   - Uses getEfFromDopingConc, getMobilitiesFromDopingConc,
%     getTaupFromDopingConc for ND-dependent properties
%
%   Parameter definitions live in sweep_config_qsspc.m (single source of truth).
%   Combination generation and filename generation are local helper functions.

%% Setup
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');

% Add necessary paths
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Transfer matrix modeling'));
addpath(fullfile(thisDir, 'Doping'));

%% Load shared configuration (single source of truth)
config = sweep_config_qsspc();

file_name       = config.file_name;
outputDir       = config.outputDir;
runs_per_chunk  = config.runs_per_chunk;
light_intensity = config.light_intensity;
Rs_ohm_cm2      = config.Rs_ohm_cm2;
t_lighton_s     = config.t_lighton_s;
tpoints         = config.tpoints;

%% Setup output directory
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

existingFiles = dir(fullfile(outputDir, 'QSSPC*.mat'));
fprintf('Found %d existing result files\n', length(existingFiles));

%% Generate all parameter combinations (local helper)
[combinations, total_combos] = generate_combinations(config);

fprintf('\nParameter space:\n');
fprintf('  ND values: %d\n', length(config.ND_values));
fprintf('  sp_r values: %d\n', length(config.sp_r_values));
fprintf('  Vi values: %d\n', length(config.Vi_values_V));
fprintf('  Total: %d combinations\n', total_combos);

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
    ND        = params(1);
    EF0       = params(2);
    sp_r      = params(3);
    Vi_V      = params(4);
    Qf_cm2    = params(5);
    Phi_right = params(6);
    sn_r      = params(7);

    % Generate filename (local helper)
    filename = generate_filename(ND, sp_r, Vi_V);
    filepath = fullfile(outputDir, filename);

    % Skip if file already exists
    if isfile(filepath)
        fprintf('  [%d/%d] Skipping (exists): %s\n', iRun, total_combos, filename);
        n_skipped_this_chunk = n_skipped_this_chunk + 1;
        continue;
    end

    fprintf('\n=== Run %d/%d (Chunk %d) ===\n', iRun, total_combos, chunk_id);
    fprintf('  ND=%.2e, EF0=%.4f, sp_r=%.0f cm/s, Vi=%.4f V, Qf=%.2e cm^-2\n', ...
            ND, EF0, sp_r, Vi_V, Qf_cm2);
    fprintf('  -> %s\n', filename);

    try
        %% Load template CSV and override sweep parameters on par
        par = pc(file_name);

        par.AbsTol = 5e-6;
        par.RelTol = 5e-3;

        % --- ND-dependent properties ---
        [mu_n, mu_p] = getMobilitiesFromDopingConc(ND);
        taup = getTaupFromDopingConc(ND);

        % Active layer overrides (dependent properties recompute from EF0)
        n_layers = length(par.EF0);
        par.EF0(n_layers)  = EF0;
        par.mu_n(n_layers) = mu_n;
        par.mu_p(n_layers) = mu_p;
        par.taup(n_layers) = taup;

        % Left electrode (rear, passivated -- flat band at this ND's EF0)
        par.Phi_left = EF0;
        par.sn_l     = 1e7;               % majority carriers, set high
        par.sp_l     = config.template.sp_l;

        % Right electrode (Tc/ZnPc surface -- swept sp_r and Vi offset)
        par.Phi_right = Phi_right;
        par.sn_r      = sn_r;
        par.sp_r      = sp_r;

        % Illumination from Tc/ZnPc side
        par.side = 'right';

        par = refresh_device(par);

        %% Equilibrate (with fallback)
        try
            solEq = equilibrate(par);
        catch ME_eq
            fprintf('  Tight tolerances failed for equilibrate, trying looser...\n');
            par.AbsTol = 1e-5;
            par.RelTol = 1e-2;
            par = refresh_device(par);
            solEq = equilibrate(par);
        end

        %% lightonRs: open circuit at 1 sun, steady state
        try
            sol_OC = lightonRs(solEq.el, light_intensity, Rs_ohm_cm2, t_lighton_s, tpoints);
        catch ME_lr
            fprintf('  Tight tolerances failed for lightonRs, trying looser...\n');
            solEq.el.par.AbsTol = 1e-5;
            solEq.el.par.RelTol = 1e-2;
            sol_OC = lightonRs(solEq.el, light_intensity, Rs_ohm_cm2, t_lighton_s, tpoints);
        end

        %% Extract tau_eff from steady-state solution
        x = sol_OC.x;
        W = x(end) - x(1);

        % Carrier densities at final time point
        try
            [n_prof, p_prof] = dfana.calcDensities(sol_OC);
            n_final = n_prof(end, :);
            p_final = p_prof(end, :);
        catch
            n_final = sol_OC.u(end, :, 1);
            p_final = sol_OC.u(end, :, 2);
        end

        % Equilibrium hole density (depends on this run's ND)
        p_eq = config.ni_Si^2 / ND;

        % Excess minority carrier density
        Delta_p = p_final - p_eq;
        Delta_p(Delta_p < 0) = 0;
        Delta_p_avg = trapz(x, Delta_p) / W;

        % Generation rate
        try
            Gx = dfana.calcG(sol_OC);
            G_final = Gx(end, :);
        catch
            q_C = 1.602e-19;
            G_final = ones(size(x)) * 0.035 / (q_C * config.template.thickness);
            warning('calcG failed for run %d, using approximate G_avg', iRun);
        end
        G_avg = trapz(x, G_final) / W;

        % Effective lifetime
        tau_eff = Delta_p_avg / G_avg;

        % Surface hole density at Tc/ZnPc interface
        p_surface = p_final(end);

        % Effective front SRV decomposition
        %   1/tau_eff = 1/tau_bulk + (S_front + S_rear) / W
        %   tau_bulk = taup for this ND (not template value)
        S_total     = (1/tau_eff - 1/taup) * W;
        S_front_eff = S_total - config.template.sp_l;

        % QFL splitting (proxy for implied Voc)
        try
            [~, ~, Efn, Efp] = dfana.calcEnergies(sol_OC);
            delta_QFL = Efn(end, end) - Efp(end, 1);
        catch
            delta_QFL = NaN;
        end

        fprintf('  tau_eff=%.3e s (%.1f us) | S_front_eff=%.1f cm/s | Delta_p_avg=%.2e cm^-3\n', ...
                tau_eff, tau_eff * 1e6, S_front_eff, Delta_p_avg);

        %% Build result structure
        thisRun = struct();

        % Sweep parameters
        thisRun.ND         = ND;
        thisRun.EF0        = EF0;
        thisRun.mu_n       = mu_n;
        thisRun.mu_p       = mu_p;
        thisRun.taup       = taup;
        thisRun.sp_r       = sp_r;
        thisRun.sn_r       = sn_r;
        thisRun.Vi_V       = Vi_V;
        thisRun.Qf_cm2     = Qf_cm2;
        thisRun.Phi_right  = Phi_right;

        % Template parameters (for reference / postprocessing)
        thisRun.template = config.template;

        % Extracted results
        thisRun.tau_eff_s    = tau_eff;
        thisRun.Delta_p_avg  = Delta_p_avg;
        thisRun.G_avg        = G_avg;
        thisRun.p_surface    = p_surface;
        thisRun.S_front_eff  = S_front_eff;
        thisRun.S_total      = S_total;
        thisRun.delta_QFL    = delta_QFL;

        % Profiles at steady state (for diagnostics)
        thisRun.x         = x;
        thisRun.n_final   = n_final;
        thisRun.p_final   = p_final;
        thisRun.Delta_p   = Delta_p;
        thisRun.G_final   = G_final;

        % Band diagram at steady state
        try
            [Ecb, Evb, Efn, Efp] = dfana.calcEnergies(sol_OC);
            thisRun.Ecb_final = Ecb(end, :);
            thisRun.Evb_final = Evb(end, :);
            thisRun.Efn_final = Efn(end, :);
            thisRun.Efp_final = Efp(end, :);
        catch
            thisRun.Ecb_final = [];
            thisRun.Evb_final = [];
            thisRun.Efn_final = [];
            thisRun.Efp_final = [];
        end

        % Metadata
        thisRun.timestamp = datetime('now');
        thisRun.par       = par;

        %% Save
        save(filepath, 'thisRun', '-v7.3');
        fprintf('  Saved: %s\n', filename);
        n_completed_this_chunk = n_completed_this_chunk + 1;

    catch ME
        fprintf('  ERROR: %s\n', ME.message);
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

warning('on', 'MATLAB:table:ModifiedAndSavedVarnames');

end


%% ========================================================================
%  HELPER FUNCTIONS
%% ========================================================================

function log_error(outputDir, filename, params, ME)
%LOG_ERROR Write error details to log file

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    errorLog = fullfile(outputDir, 'errors.txt');
    fid = fopen(errorLog, 'a');
    fprintf(fid, '\n=== Error at %s ===\n', datestr(now));
    fprintf(fid, 'File: %s\n', filename);
    fprintf(fid, 'Params: ND=%.2e, EF0=%.4f, sp_r=%.0f, Vi=%.4f, Qf=%.2e, Phi_R=%.4f\n', ...
            params(1), params(2), params(3), params(4), params(5), params(6));
    fprintf(fid, 'Error: %s\n', ME.message);
    for k = 1:length(ME.stack)
        fprintf(fid, '  %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
    end
    fclose(fid);
end


function [combinations, n_total] = generate_combinations(config)
%GENERATE_COMBINATIONS Generate all (ND, sp_r, Vi) parameter combinations
%
%   Output: combinations [N x 7] array where each row is:
%     col 1: ND         [cm^-3]
%     col 2: EF0        [eV]     precomputed from ND
%     col 3: sp_r       [cm/s]
%     col 4: Vi_V       [V]
%     col 5: Qf_cm2     [cm^-2]  from calcQfFromVi (for labeling)
%     col 6: Phi_right   [eV]    = EF0(ND) + Vi
%     col 7: sn_r       [cm/s]   = 1e7 (majority carriers, high)

    ND_vals   = config.ND_values;
    sp_r_vals = config.sp_r_values;
    Vi_vals   = config.Vi_values_V;
    Qf_vals   = config.Qf_values_cm2;

    n_ND = length(ND_vals);
    n_sp = length(sp_r_vals);
    n_Vi = length(Vi_vals);
    n_total = n_ND * n_sp * n_Vi;

    combinations = zeros(n_total, 7);
    idx = 0;
    for iND = 1:n_ND
        EF0_this = getEfFromDopingConc(ND_vals(iND));
        for iSp = 1:n_sp
            for iVi = 1:n_Vi
                idx = idx + 1;
                combinations(idx, 1) = ND_vals(iND);
                combinations(idx, 2) = EF0_this;
                combinations(idx, 3) = sp_r_vals(iSp);
                combinations(idx, 4) = Vi_vals(iVi);
                combinations(idx, 5) = Qf_vals(iVi);
                combinations(idx, 6) = EF0_this + Vi_vals(iVi);
                combinations(idx, 7) = 1e7;           % sn_r: majority carriers, set high
            end
        end
    end

    fprintf('Generated %d combinations (%d ND x %d sp_r x %d Vi)\n', ...
        n_total, n_ND, n_sp, n_Vi);
end


function filename = generate_filename(ND, sp_r, Vi_V)
%GENERATE_FILENAME Generate unique filename for a QSSPC sweep result
%   Format: QSSPC_ND<ND>_sp<sp_r>_Vi<Vi>.mat

    nd_str = sprintf('%.1e', ND);

    if sp_r >= 1 && sp_r == round(sp_r)
        sp_str = sprintf('%.0f', sp_r);
    else
        sp_str = sprintf('%.1e', sp_r);
    end

    if Vi_V == 0
        vi_str = '0';
    elseif Vi_V > 0
        vi_str = sprintf('+%.4f', Vi_V);
    else
        vi_str = sprintf('%.4f', Vi_V);
    end

    filename = sprintf('QSSPC_ND%s_sp%s_Vi%s.mat', nd_str, sp_str, vi_str);
end