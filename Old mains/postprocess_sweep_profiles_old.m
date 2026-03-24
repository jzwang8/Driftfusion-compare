function postprocess_sweep_profiles(resultsDir, loadSolutions, extractProfiles, useParallel)
%POSTPROCESS_SWEEP_PROFILES Compile all sweep results including spatial profiles
%
%   postprocess_sweep_profiles()                      - uses 'sweep_results', default options
%   postprocess_sweep_profiles(resultsDir)            - uses specified directory
%   postprocess_sweep_profiles(resultsDir, loadSol)   - loadSol=true keeps full solutions (memory intensive!)
%   postprocess_sweep_profiles(resultsDir, loadSol, extractProf) - extractProf=false skips profiles
%   postprocess_sweep_profiles(resultsDir, loadSol, extractProf, usePar) - usePar=true uses parfor
%
%   Outputs:
%     - all_runs_table_YYYYMMDD_HHMMSS.csv  : All runs with parameters + metrics + status
%     - sweep_summary_YYYYMMDD_HHMMSS.mat   : Full data including spatial profiles for plotting
%
%   The sweep_summary.mat contains:
%     - runs: struct array with all data per run:
%         * params, JV, EQE (as before)
%         * profiles: spatial data for carrier density, charge, field, bands, recombination
%     - table: summary table matching the CSV
%     - metadata: sweep info and parameter ranges
%
%   The profiles struct contains (all at equilibrium, final time point):
%     - x_nm: position grid (nm)
%     - x_sub_nm: subgrid for recombination (nm)
%     - n, p: electron and hole density (cm^-3)
%     - rho: charge density (cm^-3)
%     - F: electric field (V/cm)
%     - V: electrostatic potential (V)
%     - Ecb, Evb: conduction and valence band edges (eV)
%     - Efn, Efp: quasi-Fermi levels (eV)
%     - r_btb, r_srh, r_vsr: recombination rates (cm^-3 s^-1)

if nargin < 1 || isempty(resultsDir)
    resultsDir = 'sweep_results';
end
if nargin < 2 || isempty(loadSolutions)
    loadSolutions = false;  % Default: don't load full solution structs
end
if nargin < 3 || isempty(extractProfiles)
    extractProfiles = true;  % Default: extract spatial profiles
end
if nargin < 4 || isempty(useParallel)
    useParallel = false;  % Default: no parallel processing
end

%% Add Driftfusion/Core to path if dfana not available
if extractProfiles && ~exist('dfana', 'class')
    % Script is in Driftfusion folder, dfana is in Driftfusion/Core
    script_dir = fileparts(mfilename('fullpath'));
    core_path = fullfile(script_dir, 'Core');
    
    if isfolder(core_path)
        addpath(genpath(core_path));
        fprintf('Added to path: %s\n', core_path);
    else
        warning('Could not find Core folder. Profile extraction may fail.');
    end
end

if loadSolutions
    fprintf('WARNING: Loading full solutions - this requires a lot of memory!\n');
end

timestamp = datestr(now, 'yyyymmdd_HHMMSS');

fprintf('===========================================================\n');
fprintf(' Postprocessing Sweep Results (with Spatial Profiles)\n');
fprintf(' %s\n', datestr(now));
fprintf('===========================================================\n');
fprintf('Results directory: %s\n', resultsDir);
fprintf('Extract profiles:  %s\n', string(extractProfiles));
fprintf('Use parallel:      %s\n', string(useParallel));
fprintf('\n');

%% Define expected parameter space (must match run_sweep_chunk.m)
expected = generate_expected_combinations();
n_expected = size(expected, 1);
fprintf('Expected combinations: %d\n', n_expected);

%% Find completed result files
files = dir(fullfile(resultsDir, 'SF*.mat'));
n_files = length(files);
fprintf('Found result files: %d\n', n_files);

%% Parse error log
errorFile = fullfile(resultsDir, 'errors.txt');
errors = parse_error_log(errorFile);
fprintf('Found logged errors: %d\n', length(errors));
fprintf('\n');

%% Initialize storage (errors struct needed for parfor)
% Table columns will be extracted from runs after processing

%% Process each expected combination
fprintf('Processing runs...\n');

% Batch size to limit memory usage
batch_size = 100;  % Process 100 runs at a time
n_batches = ceil(n_expected / batch_size);

% Pre-allocate cell arrays for results
runs_cell = cell(n_expected, 1);
Jsc_cell = cell(n_expected, 1);
Voc_cell = cell(n_expected, 1);
FF_cell = cell(n_expected, 1);
PCE_cell = cell(n_expected, 1);
status_cell = cell(n_expected, 1);
error_message_cell = cell(n_expected, 1);

if useParallel
    fprintf('  Using parallel processing (parfor) in %d batches of %d\n', n_batches, batch_size);
    
    % Start pool once before batches
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool('local', 'SpmdEnabled', false);
    end
    
    for batch = 1:n_batches
        idx_start = (batch - 1) * batch_size + 1;
        idx_end = min(batch * batch_size, n_expected);
        batch_indices = idx_start:idx_end;
        batch_len = length(batch_indices);
        
        fprintf('  Batch %d/%d (runs %d-%d)...\n', batch, n_batches, idx_start, idx_end);
        
        % Temporary cell arrays for this batch
        batch_runs = cell(batch_len, 1);
        batch_Jsc = cell(batch_len, 1);
        batch_Voc = cell(batch_len, 1);
        batch_FF = cell(batch_len, 1);
        batch_PCE = cell(batch_len, 1);
        batch_status = cell(batch_len, 1);
        batch_error = cell(batch_len, 1);
        
        parfor j = 1:batch_len
            i = batch_indices(j);
            [batch_runs{j}, batch_Jsc{j}, batch_Voc{j}, batch_FF{j}, batch_PCE{j}, ...
             batch_status{j}, batch_error{j}] = ...
                process_single_run(i, expected, resultsDir, errors, loadSolutions, extractProfiles);
        end
        
        % Copy batch results to main arrays
        for j = 1:batch_len
            i = batch_indices(j);
            runs_cell{i} = batch_runs{j};
            Jsc_cell{i} = batch_Jsc{j};
            Voc_cell{i} = batch_Voc{j};
            FF_cell{i} = batch_FF{j};
            PCE_cell{i} = batch_PCE{j};
            status_cell{i} = batch_status{j};
            error_message_cell{i} = batch_error{j};
        end
        
        % Clear batch variables to free memory
        clear batch_runs batch_Jsc batch_Voc batch_FF batch_PCE batch_status batch_error;
    end
else
    fprintf('  Using sequential processing\n');
    for i = 1:n_expected
        if mod(i, 50) == 0
            fprintf('  %d / %d\n', i, n_expected);
        end
        [runs_cell{i}, Jsc_cell{i}, Voc_cell{i}, FF_cell{i}, PCE_cell{i}, ...
         status_cell{i}, error_message_cell{i}] = ...
            process_single_run(i, expected, resultsDir, errors, loadSolutions, extractProfiles);
    end
end

% Convert cell arrays back to regular arrays
runs = [runs_cell{:}];
Jsc_mA_cm2 = cell2mat(Jsc_cell);
Voc_V = cell2mat(Voc_cell);
FF_pct = cell2mat(FF_cell);
PCE_pct = cell2mat(PCE_cell);
status = status_cell;
error_message = error_message_cell;

% Extract parameter arrays from runs for table
run_number = [runs.run_number]';
SF_flag = [runs.SF_flag]';
eta_SF = [runs.eta_SF]';
junction_depth_nm = [runs.junction_depth_nm]';
SRV_cm_s = [runs.SRV_cm_s]';
N_D_cm3 = [runs.N_D_cm3]';
EF_emitter_eV = [runs.EF_emitter_eV]';
Phi_R_eV = [runs.Phi_R_eV]';
Phi_R_offset_eV = [runs.Phi_R_offset_eV]';

fprintf('  Done.\n\n');

%% Create summary table
T = table(run_number, SF_flag, eta_SF, junction_depth_nm, SRV_cm_s, ...
          N_D_cm3, EF_emitter_eV, Phi_R_eV, Phi_R_offset_eV, ...
          Jsc_mA_cm2, Voc_V, FF_pct, PCE_pct, status, error_message);

%% Save CSV
csvFile = fullfile(resultsDir, sprintf('all_runs_table_%s.csv', timestamp));
writetable(T, csvFile);
fprintf('Saved: %s\n', csvFile);

%% Build summary struct
summary = struct();
summary.runs = runs;
summary.table = T;
summary.timestamp = timestamp;
summary.resultsDir = resultsDir;
summary.n_expected = n_expected;
summary.n_success = sum(strcmp(status, 'success'));
summary.n_failed = sum(strcmp(status, 'failed'));
summary.n_not_run = sum(strcmp(status, 'not_run'));
summary.n_load_error = sum(strcmp(status, 'load_error'));

% Parameter ranges for filtering
summary.params.junction_depth_vals = unique(junction_depth_nm);
summary.params.SRV_vals = unique(SRV_cm_s);
summary.params.EF_vals = unique(EF_emitter_eV);
summary.params.Phi_R_offset_vals = unique(Phi_R_offset_eV);
summary.params.eta_SF_vals = unique(eta_SF(SF_flag == 1));

%% Save MAT
matFile = fullfile(resultsDir, sprintf('sweep_summary_%s.mat', timestamp));
fprintf('Saving summary (this may take a moment for large sweeps)...\n');
save(matFile, 'summary', '-v7.3');
fprintf('Saved: %s\n', matFile);

%% Print summary
fprintf('\n===========================================================\n');
fprintf(' Summary\n');
fprintf('===========================================================\n');
fprintf('Total expected: %d\n', n_expected);
fprintf('  Success:    %d (%.1f%%)\n', summary.n_success, 100*summary.n_success/n_expected);
fprintf('  Failed:     %d (%.1f%%)\n', summary.n_failed, 100*summary.n_failed/n_expected);
fprintf('  Not run:    %d (%.1f%%)\n', summary.n_not_run, 100*summary.n_not_run/n_expected);
fprintf('  Load error: %d (%.1f%%)\n', summary.n_load_error, 100*summary.n_load_error/n_expected);

% Stats for successful runs
success_mask = strcmp(status, 'success');
if any(success_mask)
    fprintf('\n--- Metrics (successful runs) ---\n');
    fprintf('PCE (%%):  min=%.5f  max=%.5f  mean=%.5f\n', ...
            min(PCE_pct(success_mask)), max(PCE_pct(success_mask)), mean(PCE_pct(success_mask), 'omitnan'));
    fprintf('Voc (V):  min=%.5f  max=%.5f  mean=%.5f\n', ...
            min(Voc_V(success_mask)), max(Voc_V(success_mask)), mean(Voc_V(success_mask), 'omitnan'));
    fprintf('Jsc (mA/cm²): min=%.5f  max=%.5f  mean=%.5f\n', ...
            min(Jsc_mA_cm2(success_mask)), max(Jsc_mA_cm2(success_mask)), mean(Jsc_mA_cm2(success_mask), 'omitnan'));
    fprintf('FF (%%):  min=%.5f  max=%.5f  mean=%.5f\n', ...
            min(FF_pct(success_mask)), max(FF_pct(success_mask)), mean(FF_pct(success_mask), 'omitnan'));
    
    % Top 5 by PCE
    [~, sort_idx] = sort(PCE_pct, 'descend', 'MissingPlacement', 'last');
    fprintf('\n--- Top 5 by PCE ---\n');
    for k = 1:min(5, sum(success_mask))
        idx = sort_idx(k);
        fprintf('%d. PCE=%.5f%% | SF=%d eta=%.1f jd=%d sp=%.0e EF=%.3f PhiOff=%.4f\n', ...
                k, PCE_pct(idx), SF_flag(idx), eta_SF(idx), junction_depth_nm(idx), ...
                SRV_cm_s(idx), EF_emitter_eV(idx), Phi_R_offset_eV(idx));
    end
    
    % SF benefit
    fprintf('\n--- SF Benefit ---\n');
    baseline_mask = (SF_flag == 0) & success_mask;
    sf_eta2_mask = (SF_flag == 1) & (eta_SF == 2) & success_mask;
    if any(baseline_mask) && any(sf_eta2_mask)
        fprintf('Baseline (SF=0):     n=%d, mean PCE=%.5f%%\n', sum(baseline_mask), mean(PCE_pct(baseline_mask), 'omitnan'));
        fprintf('With SF (eta=2):     n=%d, mean PCE=%.5f%%\n', sum(sf_eta2_mask), mean(PCE_pct(sf_eta2_mask), 'omitnan'));
        fprintf('Absolute PCE gain:   %.5f%%\n', mean(PCE_pct(sf_eta2_mask), 'omitnan') - mean(PCE_pct(baseline_mask), 'omitnan'));
    end
    
    % Profile extraction stats
    n_with_profiles = sum(arrayfun(@(r) ~isempty(r.profiles.x_nm), runs));
    fprintf('\n--- Profile Extraction ---\n');
    fprintf('Runs with spatial profiles: %d / %d\n', n_with_profiles, summary.n_success);
end

fprintf('\n===========================================================\n');
fprintf(' Postprocessing Complete\n');
fprintf('===========================================================\n');

end


%% ========================================================================
%  HELPER FUNCTIONS
%% ========================================================================

function combinations = generate_expected_combinations()
%GENERATE_EXPECTED_COMBINATIONS Generate all valid parameter combinations

    junction_depth_vals = [200, 250, 300, 350, 400];
    sp_r_vals = [10, 100, 1000, 10000, 100000];
    eta_sf_with_SF = [0, 0.5, 1, 1.5, 2];
    EF_vals = [-4.05, -4.104, -4.223, -4.342];
    Phi_R_offset_vals = [0.0539, 0, -0.0539];
    E_C = -4.05;
    
    EF_to_ND = containers.Map(...
        {-4.05, -4.104, -4.223, -4.342}, ...
        {7.97e19, 9.99e18, 1.00e17, 1.00e15});
    
    combinations = [];
    
    % SF_flag = 0, eta = 0
    for jd = junction_depth_vals
        for sp = sp_r_vals
            for EF = EF_vals
                for offset = Phi_R_offset_vals
                    Phi_R = EF + offset;
                    if Phi_R <= E_C
                        N_D = EF_to_ND(EF);
                        combinations = [combinations; 0, 0, jd, sp, N_D, EF, Phi_R, offset]; %#ok<AGROW>
                    end
                end
            end
        end
    end
    
    % SF_flag = 1, all eta values
    for eta = eta_sf_with_SF
        for jd = junction_depth_vals
            for sp = sp_r_vals
                for EF = EF_vals
                    for offset = Phi_R_offset_vals
                        Phi_R = EF + offset;
                        if Phi_R <= E_C
                            N_D = EF_to_ND(EF);
                            combinations = [combinations; 1, eta, jd, sp, N_D, EF, Phi_R, offset]; %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end
end


function filename = generate_filename(SF_flag, SF_efficiency, junction_depth, sp_r, N_D, Phi_R_offset)
%GENERATE_FILENAME Create deterministic filename from parameters

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
    
    % Ensure two-digit exponent
    sp_str = regexprep(sp_str, 'e(-?)(\d)$', 'e$10$2');
    nd_str = regexprep(nd_str, 'e(-?)(\d)$', 'e$10$2');
    
    filename = sprintf('%s_%s_%s_%s_%s_%s.mat', sf_str, eta_str, jd_str, sp_str, nd_str, phi_str);
end


function errors = parse_error_log(errorFile)
%PARSE_ERROR_LOG Extract failed run info from errors.txt

    errors = struct('SF_flag', {}, 'eta', {}, 'jd', {}, 'sp', {}, ...
                    'N_D', {}, 'Phi_R_offset', {}, 'filename', {}, 'error_msg', {});
    
    if ~isfile(errorFile)
        return;
    end
    
    fid = fopen(errorFile, 'r');
    if fid == -1
        return;
    end
    
    text = fread(fid, '*char')';
    fclose(fid);
    
    blocks = strsplit(text, '=== Error');
    
    for i = 2:length(blocks)
        block = blocks{i};
        
        fileMatch = regexp(block, 'File:\s*(\S+\.mat)', 'tokens', 'once');
        if isempty(fileMatch), continue; end
        
        paramMatch = regexp(block, 'Params:\s*SF=(\d+),\s*eta=([\d.]+),\s*jd=(\d+),\s*sp=([\d.e+-]+),\s*N_D=([\d.e+-]+),\s*Phi_R_off=([\d.e+-]+)', 'tokens', 'once');
        if isempty(paramMatch), continue; end
        
        errMatch = regexp(block, 'Error:\s*(.+?)(?:\n|$)', 'tokens', 'once');
        
        err = struct();
        err.SF_flag = str2double(paramMatch{1});
        err.eta = str2double(paramMatch{2});
        err.jd = str2double(paramMatch{3});
        err.sp = str2double(paramMatch{4});
        err.N_D = str2double(paramMatch{5});
        err.Phi_R_offset = str2double(paramMatch{6});
        err.filename = fileMatch{1};
        err.error_msg = '';
        if ~isempty(errMatch)
            err.error_msg = strtrim(errMatch{1});
        end
        
        errors(end+1) = err; %#ok<AGROW>
    end
end


function idx = find_error(errors, exp)
%FIND_ERROR Find matching error entry for expected parameters

    idx = [];
    tol = 1e-3;
    
    for i = 1:length(errors)
        if errors(i).SF_flag == exp(1) && ...
           abs(errors(i).eta - exp(2)) < tol && ...
           errors(i).jd == exp(3) && ...
           abs(errors(i).sp - exp(4)) / exp(4) < tol && ...
           abs(errors(i).Phi_R_offset - exp(8)) < tol
            idx = i;
            return;
        end
    end
end


function [run_data, Jsc, Voc, FF, PCE, status, error_msg] = ...
    process_single_run(i, expected, resultsDir, errors, loadSolutions, extractProfiles)
%PROCESS_SINGLE_RUN Process a single run (parfor-compatible)

    % Expected parameters
    exp = expected(i, :);
    
    % Generate expected filename
    filename = generate_filename(exp(1), exp(2), exp(3), exp(4), exp(5), exp(8));
    filepath = fullfile(resultsDir, filename);
    
    % Initialize outputs
    Jsc = NaN;
    Voc = NaN;
    FF = NaN;
    PCE = NaN;
    status = '';
    error_msg = '';
    
    % Initialize run struct with ALL possible fields
    run_data = struct();
    run_data.run_number = i;
    run_data.SF_flag = exp(1);
    run_data.eta_SF = exp(2);
    run_data.junction_depth_nm = exp(3);
    run_data.SRV_cm_s = exp(4);
    run_data.N_D_cm3 = exp(5);
    run_data.EF_emitter_eV = exp(6);
    run_data.Phi_R_eV = exp(7);
    run_data.Phi_R_offset_eV = exp(8);
    run_data.filename = filename;
    
    % Initialize result fields
    run_data.Jsc_mA_cm2 = NaN;
    run_data.Voc_V = NaN;
    run_data.FF_pct = NaN;
    run_data.PCE_pct = NaN;
    run_data.status = '';
    run_data.error_message = '';
    
    % Initialize JV struct
    run_data.JV.Vapp = [];
    run_data.JV.Jtot = [];
    run_data.JV.MPP_power = NaN;
    run_data.JV.Pin = NaN;
    
    % Initialize EQE struct
    run_data.EQE.wavelengths = [];
    run_data.EQE.EQE = [];
    run_data.EQE.IQE = [];
    run_data.EQE.Jsc_vs_lambda = [];
    
    % Initialize profiles struct
    run_data.profiles.x_nm = [];
    run_data.profiles.x_sub_nm = [];
    run_data.profiles.n = [];
    run_data.profiles.p = [];
    run_data.profiles.rho = [];
    run_data.profiles.F = [];
    run_data.profiles.V = [];
    run_data.profiles.Ecb = [];
    run_data.profiles.Evb = [];
    run_data.profiles.Efn = [];
    run_data.profiles.Efp = [];
    run_data.profiles.r_btb = [];
    run_data.profiles.r_srh = [];
    run_data.profiles.r_vsr = [];
    
    % Initialize optional fields
    run_data.timestamp = NaT;
    if loadSolutions
        run_data.par = [];
        run_data.solEq = [];
        run_data.solCV = [];
    end
    
    % Check if file exists
    if isfile(filepath)
        try
            data = load(filepath, 'thisRun');
            r = data.thisRun;
            
            % Extract metrics
            Voc = r.Voc;
            PCE = r.PCE;
            
            % Calculate Jsc (handle duplicate Vapp values)
            [V_unique, ia] = unique(r.Vapp, 'stable');
            J_unique = r.Jtot(ia);
            if length(V_unique) >= 2 && min(V_unique) <= 0 && max(V_unique) >= 0
                Jsc = abs(interp1(V_unique, J_unique, 0, 'linear'));
            else
                [~, idx_v0] = min(abs(r.Vapp));
                Jsc = abs(r.Jtot(idx_v0));
            end
            
            % Calculate FF
            if Voc > 0 && Jsc > 0
                FF = abs(r.MPP_power) / (Voc * Jsc) * 100;
            end
            
            status = 'success';
            error_msg = '';
            
            % Fill in result fields
            run_data.Jsc_mA_cm2 = Jsc;
            run_data.Voc_V = Voc;
            run_data.FF_pct = FF;
            run_data.PCE_pct = PCE;
            run_data.status = 'success';
            run_data.error_message = '';
            
            % JV data
            run_data.JV.Vapp = r.Vapp;
            run_data.JV.Jtot = r.Jtot;
            run_data.JV.MPP_power = r.MPP_power;
            run_data.JV.Pin = r.Pin;
            
            % EQE data
            run_data.EQE.wavelengths = r.EQE_wavelengths;
            run_data.EQE.EQE = r.EQE;
            run_data.EQE.IQE = r.IQE;
            run_data.EQE.Jsc_vs_lambda = r.Jsc_vs_lambda;
            
            % ============================================================
            % EXTRACT SPATIAL PROFILES
            % ============================================================
            if extractProfiles && isfield(r, 'solEq') && isfield(r.solEq, 'el') && ~isempty(r.solEq.el)
                try
                    sol = r.solEq.el;
                    
                    % Extract using dfana
                    [~, ~, x, par, ~, n_arr, p_arr, ~, ~, V_arr] = dfana.splitsol(sol);
                    
                    rho_arr = dfana.calcrho(sol, "whole");
                    F_arr   = dfana.calcF(sol, "whole");
                    
                    [Ecb_arr, Evb_arr, Efn_arr, Efp_arr] = dfana.calcEnergies(sol);
                    
                    rsub  = dfana.calcr(sol, "sub");
                    x_sub = par.x_sub;
                    
                    % Get final time index
                    t_idx = size(n_arr, 1);
                    
                    % Store profiles (convert x from cm to nm)
                    run_data.profiles.x_nm     = x * 1e7;
                    run_data.profiles.x_sub_nm = x_sub * 1e7;
                    run_data.profiles.n        = n_arr(t_idx, :);
                    run_data.profiles.p        = p_arr(t_idx, :);
                    run_data.profiles.rho      = rho_arr(t_idx, :);
                    run_data.profiles.F        = F_arr(t_idx, :);
                    run_data.profiles.V        = V_arr(t_idx, :);
                    run_data.profiles.Ecb      = Ecb_arr(t_idx, :);
                    run_data.profiles.Evb      = Evb_arr(t_idx, :);
                    run_data.profiles.Efn      = Efn_arr(t_idx, :);
                    run_data.profiles.Efp      = Efp_arr(t_idx, :);
                    run_data.profiles.r_btb    = rsub.btb(t_idx, :);
                    run_data.profiles.r_srh    = rsub.srh(t_idx, :);
                    run_data.profiles.r_vsr    = rsub.vsr(t_idx, :);
                    
                catch
                    % Silently skip profile extraction errors in parfor
                end
            end
            % ============================================================
            
            % Full solutions (only if requested)
            if loadSolutions
                run_data.par = r.par;
                run_data.solEq = r.solEq;
                run_data.solCV = r.solCV;
            end
            if isfield(r, 'timestamp')
                run_data.timestamp = r.timestamp;
            end
            
        catch ME
            status = 'load_error';
            error_msg = ME.message;
            run_data.status = 'load_error';
            run_data.error_message = ME.message;
        end
    else
        % Check if in error log
        err_idx = find_error(errors, exp);
        if ~isempty(err_idx)
            status = 'failed';
            error_msg = errors(err_idx).error_msg;
            run_data.status = 'failed';
            run_data.error_message = errors(err_idx).error_msg;
        else
            status = 'not_run';
            error_msg = '';
            run_data.status = 'not_run';
            run_data.error_message = '';
        end
    end
end