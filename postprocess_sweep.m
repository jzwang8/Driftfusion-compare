% doesn't work right now

function postprocess_sweep(resultsDir, loadSolutions)
%POSTPROCESS_SWEEP Compile all sweep results into summary files
%
%   postprocess_sweep()                      - uses 'sweep_results', no solutions
%   postprocess_sweep(resultsDir)            - uses specified directory, no solutions
%   postprocess_sweep(resultsDir, true)      - load full solutions (memory intensive!)
%
%   Outputs:
%     - all_runs_table_YYYYMMDD_HHMMSS.csv  : All runs with parameters + metrics + status
%     - sweep_summary_YYYYMMDD_HHMMSS.mat   : Full data for plotting/analysis
%
%   The sweep_summary.mat contains:
%     - runs: struct array with all data per run (params, JV, EQE, SR, optionally solutions)
%     - table: summary table matching the CSV
%     - metadata: sweep info and parameter ranges

if nargin < 1
    resultsDir = 'sweep_results';
end
if nargin < 2
    loadSolutions = false;  % Default: don't load heavy data
end

if loadSolutions
    fprintf('WARNING: Loading full solutions - this requires a lot of memory!\n');
end

timestamp = datestr(now, 'yyyymmdd_HHMMSS');

fprintf('===========================================================\n');
fprintf(' Postprocessing Sweep Results\n');
fprintf(' %s\n', datestr(now));
fprintf('===========================================================\n');
fprintf('Results directory: %s\n\n', resultsDir);

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

%% Initialize storage
runs = struct([]);  % Will hold all run data

% Table columns
run_number = zeros(n_expected, 1);
SF_flag = zeros(n_expected, 1);
eta_SF = zeros(n_expected, 1);
junction_depth_nm = zeros(n_expected, 1);
SRV_cm_s = zeros(n_expected, 1);
N_D_cm3 = zeros(n_expected, 1);
EF_emitter_eV = zeros(n_expected, 1);
Phi_R_eV = zeros(n_expected, 1);
Phi_R_offset_eV = zeros(n_expected, 1);

Jsc_mA_cm2 = NaN(n_expected, 1);
Voc_V = NaN(n_expected, 1);
FF_pct = NaN(n_expected, 1);
PCE_pct = NaN(n_expected, 1);

status = cell(n_expected, 1);
error_message = cell(n_expected, 1);

%% Process each expected combination
fprintf('Processing runs...\n');
for i = 1:n_expected
    if mod(i, 50) == 0
        fprintf('  %d / %d\n', i, n_expected);
    end
    
    % Expected parameters
    exp = expected(i, :);
    run_number(i) = i;
    SF_flag(i) = exp(1);
    eta_SF(i) = exp(2);
    junction_depth_nm(i) = exp(3);
    SRV_cm_s(i) = exp(4);
    N_D_cm3(i) = exp(5);
    EF_emitter_eV(i) = exp(6);
    Phi_R_eV(i) = exp(7);
    Phi_R_offset_eV(i) = exp(8);
    
    % Generate expected filename
    filename = generate_filename(exp(1), exp(2), exp(3), exp(4), exp(5), exp(8));
    filepath = fullfile(resultsDir, filename);
    
    % Initialize run struct with ALL possible fields (avoids dissimilar struct error)
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
    
    % Initialize result fields (will be filled if successful)
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
    
    % Initialize SR (spectral response under bias) struct
    % SR.J:         [n_wavelengths x n_bias] - J in A/m2 at each bias
    % SR.EQE:       [n_wavelengths x n_bias] - EQE-equivalent at each bias (%)
    % SR.bias_vals: [1 x n_bias]             - applied bias values in V
    run_data.SR.bias_vals = [];
    run_data.SR.J = [];
    run_data.SR.EQE = [];
    
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
            clear data;  % Free memory
            
            % Extract metrics
            Voc_V(i) = r.Voc;
            PCE_pct(i) = r.PCE;
            
            % Calculate Jsc (handle duplicate Vapp values)
            [V_unique, ia] = unique(r.Vapp, 'stable');
            J_unique = r.Jtot(ia);
            if length(V_unique) >= 2 && min(V_unique) <= 0 && max(V_unique) >= 0
                Jsc_mA_cm2(i) = abs(interp1(V_unique, J_unique, 0, 'linear'));
            else
                [~, idx_v0] = min(abs(r.Vapp));
                Jsc_mA_cm2(i) = abs(r.Jtot(idx_v0));
            end
            
            % Calculate FF
            if Voc_V(i) > 0 && Jsc_mA_cm2(i) > 0
                FF_pct(i) = abs(r.MPP_power) / (Voc_V(i) * Jsc_mA_cm2(i)) * 100;
            end
            
            status{i} = 'success';
            error_message{i} = '';
            
            % Fill in result fields
            run_data.Jsc_mA_cm2 = Jsc_mA_cm2(i);
            run_data.Voc_V = Voc_V(i);
            run_data.FF_pct = FF_pct(i);
            run_data.PCE_pct = PCE_pct(i);
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
            
            % Spectral response under bias (only present in newer runs)
            if isfield(r, 'SR_bias_vals') && isfield(r, 'SR') && isfield(r, 'SR_EQE')
                run_data.SR.bias_vals = r.SR_bias_vals;
                run_data.SR.J = r.SR;
                run_data.SR.EQE = r.SR_EQE;
            end
            
            % Solutions (only if requested)
            if loadSolutions
                run_data.par = r.par;
                run_data.solEq = r.solEq;
                run_data.solCV = r.solCV;
            end
            if isfield(r, 'timestamp')
                run_data.timestamp = r.timestamp;
            end
            
            clear r;  % Free memory
            
        catch ME
            status{i} = 'load_error';
            error_message{i} = ME.message;
            run_data.status = 'load_error';
            run_data.error_message = ME.message;
        end
    else
        % Check if in error log
        err_idx = find_error(errors, exp);
        if ~isempty(err_idx)
            status{i} = 'failed';
            error_message{i} = errors(err_idx).error_msg;
            run_data.status = 'failed';
            run_data.error_message = errors(err_idx).error_msg;
        else
            status{i} = 'not_run';
            error_message{i} = '';
            run_data.status = 'not_run';
            run_data.error_message = '';
        end
    end
    
    % Append to runs array
    if i == 1
        runs = run_data;
    else
        runs(i) = run_data;
    end
end
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

% Check if any runs have SR data and record bias vals used
sr_mask = strcmp(status, 'success') & ...
          arrayfun(@(r) ~isempty(r.SR.bias_vals), runs)';
if any(sr_mask)
    summary.params.SR_bias_vals = runs(find(sr_mask, 1)).SR.bias_vals;
    fprintf('SR data found in %d / %d successful runs\n', sum(sr_mask), summary.n_success);
else
    summary.params.SR_bias_vals = [];
    fprintf('No SR data found (older runs without spectral response under bias)\n');
end

%% Save MAT
matFile = fullfile(resultsDir, sprintf('sweep_summary_%s.mat', timestamp));
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
    
    % SR summary (if available)
    if any(sr_mask) && ~isempty(summary.params.SR_bias_vals)
        fprintf('\n--- Spectral Response Under Bias ---\n');
        fprintf('Bias points: %s V\n', mat2str(summary.params.SR_bias_vals));
        fprintf('Runs with SR data: %d\n', sum(sr_mask));
        % Print mean EQE at each bias for SF=0 vs SF=1 eta=2
        sr_baseline = sr_mask' & (SF_flag == 0);
        sr_sf2 = sr_mask' & (SF_flag == 1) & (eta_SF == 2);
        for ib = 1:length(summary.params.SR_bias_vals)
            bias = summary.params.SR_bias_vals(ib);
            if any(sr_baseline)
                eqe_vals = arrayfun(@(r) mean(r.SR.EQE(:, ib), 'omitnan'), runs(sr_baseline));
                fprintf('  V=%.1fV | SF=0 mean EQE=%.2f%%', bias, mean(eqe_vals, 'omitnan'));
            end
            if any(sr_sf2)
                eqe_vals = arrayfun(@(r) mean(r.SR.EQE(:, ib), 'omitnan'), runs(sr_sf2));
                fprintf(' | SF=1 eta=2 mean EQE=%.2f%%', mean(eqe_vals, 'omitnan'));
            end
            fprintf('\n');
        end
    end
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
%   Must match the logic in run_sweep_chunk.m exactly

    junction_depth_vals = [200];
    sp_r_vals = [1000, 10000, 100000];
    eta_sf_with_SF = [0, 1, 2];
    EF_vals = [-4.05, -4.104, -4.223, -4.342];
    Phi_R_offset_vals = [0.054, 0, -0.054];
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