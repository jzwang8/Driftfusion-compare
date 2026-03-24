function postprocess_sweep_profiles(resultsDir, loadSolutions)
%POSTPROCESS_SWEEP_PROFILES Compile all profile-sweep results into summary files
%
%   postprocess_sweep_profiles()                      - uses default dir, no solutions
%   postprocess_sweep_profiles(resultsDir)            - uses specified directory, no solutions
%   postprocess_sweep_profiles(resultsDir, true)      - load full solutions (memory intensive!)
%
%   Parameter definitions come from sweep_config_profiles.m (single source of truth).
%   Combination generation comes from generate_sweep_combinations.m.
%   Filename generation comes from generate_sweep_filename.m.
%
%   Outputs:
%     - all_runs_table_YYYYMMDD_HHMMSS.csv  : All runs with parameters + metrics + status
%     - sweep_summary_YYYYMMDD_HHMMSS.mat   : Full data for plotting/analysis
%
%   The sweep_summary.mat contains:
%     - runs: struct array with all data per run
%         .JV       : J-V curve (Vapp, Jtot, MPP_power, Pin)
%         .EQE      : EQE spectrum (wavelengths, EQE, IQE, Jsc_vs_lambda)
%         .SR       : spectral response under bias (bias_vals, J, EQE)
%         .profile  : emitter doping profile (x_nm, N_D, EF, mu_n, mu_p, taup)
%         .N0_peak, .EF_peak, .EF_surface : scalar doping summary
%         optionally: par, solEq, solCV (if loadSolutions=true)
%     - table: summary table matching the CSV
%     - metadata: sweep info and parameter ranges

%% Load shared configuration (single source of truth)
config = sweep_config_profiles();
profile_shapes = config.profile_shapes;

if nargin < 1
    resultsDir = config.outputDir;
end
if nargin < 2
    loadSolutions = false;
end

if loadSolutions
    fprintf('WARNING: Loading full solutions - this requires a lot of memory!\n');
end

timestamp = datestr(now, 'yyyymmdd_HHMMSS');

fprintf('===========================================================\n');
fprintf(' Postprocessing Profile Sweep Results\n');
fprintf(' %s\n', datestr(now));
fprintf('===========================================================\n');
fprintf('Results directory: %s\n\n', resultsDir);

%% Generate expected parameter space (shared logic — identical to runner)
[expected, expected_profile_names] = generate_sweep_combinations(config);
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
runs = struct([]);

% Table columns
run_number        = zeros(n_expected, 1);
SF_flag           = zeros(n_expected, 1);
eta_SF            = zeros(n_expected, 1);
junction_depth_nm = zeros(n_expected, 1);
SRV_cm_s          = zeros(n_expected, 1);
N0_peak_cm3       = zeros(n_expected, 1);
EF_peak_eV        = zeros(n_expected, 1);
EF_surface_eV     = zeros(n_expected, 1);
Phi_R_eV          = zeros(n_expected, 1);
Phi_R_offset_eV   = zeros(n_expected, 1);
profile_type_col  = cell(n_expected, 1);

Jsc_mA_cm2 = NaN(n_expected, 1);
Voc_V      = NaN(n_expected, 1);
FF_pct     = NaN(n_expected, 1);
PCE_pct    = NaN(n_expected, 1);

status        = cell(n_expected, 1);
error_message = cell(n_expected, 1);

%% Process each expected combination
fprintf('Processing runs...\n');
for i = 1:n_expected
    if mod(i, 50) == 0
        fprintf('  %d / %d\n', i, n_expected);
    end

    % Expected parameters (10 columns):
    %  1:SF_flag  2:eta_sf  3:jd  4:sp  5:N0_peak  6:EF_peak
    %  7:Phi_R    8:Phi_R_offset  9:profile_idx  10:EF_surface
    exp = expected(i, :);
    prof_name = expected_profile_names{i};

    run_number(i)        = i;
    SF_flag(i)           = exp(1);
    eta_SF(i)            = exp(2);
    junction_depth_nm(i) = exp(3);
    SRV_cm_s(i)          = exp(4);
    N0_peak_cm3(i)       = exp(5);
    EF_peak_eV(i)        = exp(6);
    Phi_R_eV(i)          = exp(7);
    Phi_R_offset_eV(i)   = exp(8);
    EF_surface_eV(i)     = exp(10);
    profile_type_col{i}  = prof_name;

    % Generate expected filename (shared function — identical to runner)
    filename = generate_sweep_filename(exp(1), exp(2), exp(3), exp(4), exp(5), exp(8), prof_name);
    filepath = fullfile(resultsDir, filename);

    % -----------------------------------------------------------------
    % Initialize run struct with ALL possible fields
    % (avoids dissimilar struct error when appending)
    % -----------------------------------------------------------------
    run_data = struct();
    run_data.run_number        = i;
    run_data.SF_flag           = exp(1);
    run_data.eta_SF            = exp(2);
    run_data.junction_depth_nm = exp(3);
    run_data.SRV_cm_s          = exp(4);
    run_data.N0_peak           = exp(5);
    run_data.EF_peak           = exp(6);
    run_data.EF_surface        = exp(10);
    run_data.Phi_R             = exp(7);
    run_data.Phi_R_offset      = exp(8);
    run_data.profile_type      = prof_name;
    run_data.profile_idx       = exp(9);
    run_data.N_sublayers       = NaN;
    run_data.filename          = filename;

    % Metrics (filled on success)
    run_data.Jsc_mA_cm2   = NaN;
    run_data.Voc_V        = NaN;
    run_data.FF_pct       = NaN;
    run_data.PCE_pct      = NaN;
    run_data.status       = '';
    run_data.error_message = '';

    % JV struct
    run_data.JV.Vapp      = [];
    run_data.JV.Jtot      = [];
    run_data.JV.MPP_power = NaN;
    run_data.JV.Pin       = NaN;

    % EQE struct
    run_data.EQE.wavelengths   = [];
    run_data.EQE.EQE           = [];
    run_data.EQE.IQE           = [];
    run_data.EQE.Jsc_vs_lambda = [];

    % Spectral response under bias
    run_data.SR.bias_vals = [];
    run_data.SR.J         = [];   % [n_wavelengths x n_bias]
    run_data.SR.EQE       = [];   % [n_wavelengths x n_bias]

    % Doping profile arrays
    run_data.N0_peak_scalar  = NaN;   % redundant scalar (for quick access)
    run_data.EF_peak_scalar  = NaN;
    run_data.EF_surface_scalar = NaN;
    run_data.profile.x_nm  = [];
    run_data.profile.N_D   = [];
    run_data.profile.EF    = [];
    run_data.profile.mu_n  = [];
    run_data.profile.mu_p  = [];
    run_data.profile.taup  = [];

    % Optional heavy data
    run_data.timestamp = NaT;
    if loadSolutions
        run_data.par   = [];
        run_data.solEq = [];
        run_data.solCV = [];
    end

    % -----------------------------------------------------------------
    % Load result file
    % -----------------------------------------------------------------
    if isfile(filepath)
        try
            data = load(filepath, 'thisRun');
            r = data.thisRun;
            clear data;

            % --- Core metrics ---
            Voc_V(i)   = r.Voc;
            PCE_pct(i) = r.PCE;

            % Jsc (handle duplicate Vapp values)
            [V_unique, ia] = unique(r.Vapp, 'stable');
            J_unique = r.Jtot(ia);
            if length(V_unique) >= 2 && min(V_unique) <= 0 && max(V_unique) >= 0
                Jsc_mA_cm2(i) = abs(interp1(V_unique, J_unique, 0, 'linear'));
            else
                [~, idx_v0] = min(abs(r.Vapp));
                Jsc_mA_cm2(i) = abs(r.Jtot(idx_v0));
            end

            % FF
            if Voc_V(i) > 0 && Jsc_mA_cm2(i) > 0
                FF_pct(i) = abs(r.MPP_power) / (Voc_V(i) * Jsc_mA_cm2(i)) * 100;
            end

            status{i}        = 'success';
            error_message{i} = '';

            run_data.Jsc_mA_cm2    = Jsc_mA_cm2(i);
            run_data.Voc_V         = Voc_V(i);
            run_data.FF_pct        = FF_pct(i);
            run_data.PCE_pct       = PCE_pct(i);
            run_data.status        = 'success';
            run_data.error_message = '';

            % --- JV data ---
            run_data.JV.Vapp      = r.Vapp;
            run_data.JV.Jtot      = r.Jtot;
            run_data.JV.MPP_power = r.MPP_power;
            run_data.JV.Pin       = r.Pin;

            % --- EQE data ---
            run_data.EQE.wavelengths   = r.EQE_wavelengths;
            run_data.EQE.EQE           = r.EQE;
            run_data.EQE.IQE           = r.IQE;
            run_data.EQE.Jsc_vs_lambda = r.Jsc_vs_lambda;

            % --- Spectral response under bias ---
            if isfield(r, 'SR_bias_vals') && isfield(r, 'SR') && isfield(r, 'SR_EQE')
                run_data.SR.bias_vals = r.SR_bias_vals;
                run_data.SR.J         = r.SR;        % [n_wl x n_bias]
                run_data.SR.EQE       = r.SR_EQE;    % [n_wl x n_bias]
            end

            % --- Doping profile scalar summaries ---
            if isfield(r, 'N0_peak'),    run_data.N0_peak_scalar    = r.N0_peak;    end
            if isfield(r, 'EF_peak'),    run_data.EF_peak_scalar    = r.EF_peak;    end
            if isfield(r, 'EF_surface'), run_data.EF_surface_scalar = r.EF_surface; end

            % --- Doping profile arrays ---
            if isfield(r, 'profile')
                if isfield(r.profile, 'x_nm'),  run_data.profile.x_nm  = r.profile.x_nm;  end
                if isfield(r.profile, 'N_D'),   run_data.profile.N_D   = r.profile.N_D;   end
                if isfield(r.profile, 'EF'),    run_data.profile.EF    = r.profile.EF;    end
                if isfield(r.profile, 'mu_n'),  run_data.profile.mu_n  = r.profile.mu_n;  end
                if isfield(r.profile, 'mu_p'),  run_data.profile.mu_p  = r.profile.mu_p;  end
                if isfield(r.profile, 'taup'),  run_data.profile.taup  = r.profile.taup;  end
            end

            % --- Profile type from saved data (cross-check) ---
            if isfield(r, 'profile_type')
                run_data.profile_type = r.profile_type;
            end
            if isfield(r, 'N_sublayers')
                run_data.N_sublayers = r.N_sublayers;
            end

            % --- Solutions (only if requested) ---
            if loadSolutions
                run_data.par   = r.par;
                run_data.solEq = r.solEq;
                run_data.solCV = r.solCV;
            end
            if isfield(r, 'timestamp')
                run_data.timestamp = r.timestamp;
            end

            clear r;

        catch ME
            status{i}        = 'load_error';
            error_message{i} = ME.message;
            run_data.status        = 'load_error';
            run_data.error_message = ME.message;
        end
    else
        % Check if in error log
        err_idx = find_error(errors, exp, prof_name);
        if ~isempty(err_idx)
            status{i}        = 'failed';
            error_message{i} = errors(err_idx).error_msg;
            run_data.status        = 'failed';
            run_data.error_message = errors(err_idx).error_msg;
        else
            status{i}        = 'not_run';
            error_message{i} = '';
            run_data.status        = 'not_run';
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

%% Second pass: pick up any result files not matched by expected combinations
%  (safety net for legacy files or edge cases — should be empty now that
%   both scripts share the same combination generator)
matched_files = {runs.filename};
unmatched_files = {};
for f = 1:n_files
    fn = files(f).name;
    if ~any(strcmp(matched_files, fn))
        unmatched_files{end+1} = fn; %#ok<AGROW>
    end
end

if ~isempty(unmatched_files)
    fprintf('\nFound %d unmatched result files — loading...\n', length(unmatched_files));
    n_before = length(runs);
    for uf = 1:length(unmatched_files)
        fn = unmatched_files{uf};
        fp = fullfile(resultsDir, fn);
        try
            data = load(fp, 'thisRun');
            r = data.thisRun;
            clear data;

            idx = n_before + uf;

            % Build run_data from saved struct fields
            run_data = struct();
            run_data.run_number        = idx;
            run_data.SF_flag           = double(r.SF_flag);
            run_data.eta_SF            = r.SF_efficiency;
            run_data.junction_depth_nm = r.junction_depth;
            run_data.SRV_cm_s          = r.sp_r;
            run_data.N0_peak           = r.N0_peak;
            run_data.EF_peak           = r.EF_peak;
            run_data.EF_surface        = r.EF_surface;
            run_data.Phi_R             = r.Phi_R;
            run_data.Phi_R_offset      = r.Phi_R_offset;
            if isfield(r, 'profile_type')
                run_data.profile_type  = r.profile_type;
            else
                run_data.profile_type  = 'unknown';
            end
            % Recover profile_idx from name
            pidx = find(strcmp(profile_shapes, run_data.profile_type), 1);
            if isempty(pidx), pidx = 0; end
            run_data.profile_idx       = pidx;
            run_data.N_sublayers       = NaN;
            if isfield(r, 'N_sublayers'), run_data.N_sublayers = r.N_sublayers; end
            run_data.filename          = fn;

            % Metrics
            run_data.Voc_V   = r.Voc;
            run_data.PCE_pct = r.PCE;

            [V_unique, ia] = unique(r.Vapp, 'stable');
            J_unique = r.Jtot(ia);
            if length(V_unique) >= 2 && min(V_unique) <= 0 && max(V_unique) >= 0
                run_data.Jsc_mA_cm2 = abs(interp1(V_unique, J_unique, 0, 'linear'));
            else
                [~, idx_v0] = min(abs(r.Vapp));
                run_data.Jsc_mA_cm2 = abs(r.Jtot(idx_v0));
            end
            if r.Voc > 0 && run_data.Jsc_mA_cm2 > 0
                run_data.FF_pct = abs(r.MPP_power) / (r.Voc * run_data.Jsc_mA_cm2) * 100;
            else
                run_data.FF_pct = NaN;
            end

            run_data.status        = 'success';
            run_data.error_message = '';

            % JV
            run_data.JV.Vapp      = r.Vapp;
            run_data.JV.Jtot      = r.Jtot;
            run_data.JV.MPP_power = r.MPP_power;
            run_data.JV.Pin       = r.Pin;

            % EQE
            run_data.EQE.wavelengths   = r.EQE_wavelengths;
            run_data.EQE.EQE           = r.EQE;
            run_data.EQE.IQE           = r.IQE;
            run_data.EQE.Jsc_vs_lambda = r.Jsc_vs_lambda;

            % SR
            run_data.SR.bias_vals = [];
            run_data.SR.J         = [];
            run_data.SR.EQE       = [];
            if isfield(r, 'SR_bias_vals') && isfield(r, 'SR') && isfield(r, 'SR_EQE')
                run_data.SR.bias_vals = r.SR_bias_vals;
                run_data.SR.J         = r.SR;
                run_data.SR.EQE       = r.SR_EQE;
            end

            % Profile scalars
            run_data.N0_peak_scalar    = NaN;
            run_data.EF_peak_scalar    = NaN;
            run_data.EF_surface_scalar = NaN;
            if isfield(r, 'N0_peak'),    run_data.N0_peak_scalar    = r.N0_peak;    end
            if isfield(r, 'EF_peak'),    run_data.EF_peak_scalar    = r.EF_peak;    end
            if isfield(r, 'EF_surface'), run_data.EF_surface_scalar = r.EF_surface; end

            % Profile arrays
            run_data.profile.x_nm  = [];
            run_data.profile.N_D   = [];
            run_data.profile.EF    = [];
            run_data.profile.mu_n  = [];
            run_data.profile.mu_p  = [];
            run_data.profile.taup  = [];
            if isfield(r, 'profile')
                if isfield(r.profile, 'x_nm'),  run_data.profile.x_nm  = r.profile.x_nm;  end
                if isfield(r.profile, 'N_D'),   run_data.profile.N_D   = r.profile.N_D;   end
                if isfield(r.profile, 'EF'),    run_data.profile.EF    = r.profile.EF;    end
                if isfield(r.profile, 'mu_n'),  run_data.profile.mu_n  = r.profile.mu_n;  end
                if isfield(r.profile, 'mu_p'),  run_data.profile.mu_p  = r.profile.mu_p;  end
                if isfield(r.profile, 'taup'),  run_data.profile.taup  = r.profile.taup;  end
            end

            run_data.timestamp = NaT;
            if isfield(r, 'timestamp'), run_data.timestamp = r.timestamp; end
            if loadSolutions
                run_data.par   = r.par;
                run_data.solEq = r.solEq;
                run_data.solCV = r.solCV;
            end

            runs(idx) = run_data;

            % Extend table vectors
            run_number(idx)        = idx;
            SF_flag(idx)           = run_data.SF_flag;
            eta_SF(idx)            = run_data.eta_SF;
            junction_depth_nm(idx) = run_data.junction_depth_nm;
            SRV_cm_s(idx)          = run_data.SRV_cm_s;
            N0_peak_cm3(idx)       = run_data.N0_peak;
            EF_peak_eV(idx)        = run_data.EF_peak;
            EF_surface_eV(idx)     = run_data.EF_surface;
            Phi_R_eV(idx)          = run_data.Phi_R;
            Phi_R_offset_eV(idx)   = run_data.Phi_R_offset;
            profile_type_col{idx}  = run_data.profile_type;
            Jsc_mA_cm2(idx)        = run_data.Jsc_mA_cm2;
            Voc_V(idx)             = run_data.Voc_V;
            FF_pct(idx)            = run_data.FF_pct;
            PCE_pct(idx)           = run_data.PCE_pct;
            status{idx}            = 'success';
            error_message{idx}     = '';

            fprintf('  Loaded unmatched file: %s (PCE=%.4f%%)\n', fn, run_data.PCE_pct);
            clear r;
        catch ME
            fprintf('  Failed to load unmatched file %s: %s\n', fn, ME.message);
        end
    end
    n_expected = length(runs);  % update total count
    fprintf('  Total runs after second pass: %d\n', n_expected);
end

%% Create summary table
T = table(run_number, SF_flag, eta_SF, junction_depth_nm, SRV_cm_s, ...
          N0_peak_cm3, EF_peak_eV, EF_surface_eV, Phi_R_eV, Phi_R_offset_eV, ...
          profile_type_col, ...
          Jsc_mA_cm2, Voc_V, FF_pct, PCE_pct, status, error_message);

%% Save CSV
csvFile = fullfile(resultsDir, sprintf('all_runs_table_%s.csv', timestamp));
writetable(T, csvFile);
fprintf('Saved: %s\n', csvFile);

%% Build summary struct
summary = struct();
summary.runs       = runs;
summary.table      = T;
summary.timestamp  = timestamp;
summary.resultsDir = resultsDir;
summary.n_expected   = n_expected;
summary.n_success    = sum(strcmp(status, 'success'));
summary.n_failed     = sum(strcmp(status, 'failed'));
summary.n_not_run    = sum(strcmp(status, 'not_run'));
summary.n_load_error = sum(strcmp(status, 'load_error'));

% Parameter ranges for filtering
summary.params.junction_depth_vals = unique(junction_depth_nm);
summary.params.SRV_vals            = unique(SRV_cm_s);
summary.params.EF_peak_vals        = unique(EF_peak_eV);
summary.params.Phi_R_offset_vals   = unique(Phi_R_offset_eV);
summary.params.eta_SF_vals         = unique(eta_SF(SF_flag == 1));
summary.params.profile_shapes      = profile_shapes;

% --- SR data availability ---
sr_mask = strcmp(status, 'success') & ...
          arrayfun(@(r) ~isempty(r.SR.bias_vals), runs)';
if any(sr_mask)
    summary.params.SR_bias_vals = runs(find(sr_mask, 1)).SR.bias_vals;
    fprintf('SR data found in %d / %d successful runs\n', sum(sr_mask), summary.n_success);
else
    summary.params.SR_bias_vals = [];
    fprintf('No SR data found (older runs without spectral response under bias)\n');
end

% --- Doping profile availability ---
profile_mask = strcmp(status, 'success') & ...
               arrayfun(@(r) ~isempty(r.profile.N_D), runs)';
if any(profile_mask)
    fprintf('Doping profile data found in %d / %d successful runs\n', sum(profile_mask), summary.n_success);
else
    fprintf('No doping profile data found\n');
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

success_mask = strcmp(status, 'success');
if any(success_mask)
    fprintf('\n--- Metrics (successful runs) ---\n');
    fprintf('PCE (%%):  min=%.5f  max=%.5f  mean=%.5f\n', ...
            min(PCE_pct(success_mask)), max(PCE_pct(success_mask)), mean(PCE_pct(success_mask), 'omitnan'));
    fprintf('Voc (V):  min=%.5f  max=%.5f  mean=%.5f\n', ...
            min(Voc_V(success_mask)), max(Voc_V(success_mask)), mean(Voc_V(success_mask), 'omitnan'));
    fprintf('Jsc (mA/cm2): min=%.5f  max=%.5f  mean=%.5f\n', ...
            min(Jsc_mA_cm2(success_mask)), max(Jsc_mA_cm2(success_mask)), mean(Jsc_mA_cm2(success_mask), 'omitnan'));
    fprintf('FF (%%):  min=%.5f  max=%.5f  mean=%.5f\n', ...
            min(FF_pct(success_mask)), max(FF_pct(success_mask)), mean(FF_pct(success_mask), 'omitnan'));

    % Top 5 by PCE
    [~, sort_idx] = sort(PCE_pct, 'descend', 'MissingPlacement', 'last');
    fprintf('\n--- Top 5 by PCE ---\n');
    for k = 1:min(5, sum(success_mask))
        idx = sort_idx(k);
        fprintf('%d. PCE=%.5f%% | SF=%d eta=%.1f jd=%d sp=%.0e N0=%.0e prof=%s EF_pk=%.3f EF_sf=%.4f PhiOff=%.4f\n', ...
                k, PCE_pct(idx), SF_flag(idx), eta_SF(idx), junction_depth_nm(idx), ...
                SRV_cm_s(idx), N0_peak_cm3(idx), profile_type_col{idx}, ...
                EF_peak_eV(idx), EF_surface_eV(idx), Phi_R_offset_eV(idx));
    end

    % SF benefit (per profile shape)
    fprintf('\n--- SF Benefit by Profile Shape ---\n');
    for ip = 1:length(profile_shapes)
        pname = profile_shapes{ip};
        prof_mask = strcmp(profile_type_col, pname);
        baseline_mask = (SF_flag == 0) & success_mask & prof_mask;
        sf_eta2_mask  = (SF_flag == 1) & (eta_SF == 2) & success_mask & prof_mask;
        if any(baseline_mask) && any(sf_eta2_mask)
            bl_pce = mean(PCE_pct(baseline_mask), 'omitnan');
            sf_pce = mean(PCE_pct(sf_eta2_mask), 'omitnan');
            fprintf('  %s: Baseline n=%d mean=%.5f%% | SF(eta=2) n=%d mean=%.5f%% | gain=%.5f%%\n', ...
                    pname, sum(baseline_mask), bl_pce, sum(sf_eta2_mask), sf_pce, sf_pce - bl_pce);
        end
    end

    % SR summary (if available)
    if any(sr_mask) && ~isempty(summary.params.SR_bias_vals)
        fprintf('\n--- Spectral Response Under Bias ---\n');
        fprintf('Bias points: %s V\n', mat2str(summary.params.SR_bias_vals));
        fprintf('Runs with SR data: %d\n', sum(sr_mask));
        sr_baseline = sr_mask' & (SF_flag == 0);
        sr_sf2      = sr_mask' & (SF_flag == 1) & (eta_SF == 2);
        for ib = 1:length(summary.params.SR_bias_vals)
            bias = summary.params.SR_bias_vals(ib);
            if any(sr_baseline)
                eqe_vals = arrayfun(@(r) mean(r.SR.EQE(:, ib), 'omitnan'), runs(sr_baseline));
                fprintf('  V=%.2fV | SF=0 mean EQE=%.2f%%', bias, mean(eqe_vals, 'omitnan'));
            end
            if any(sr_sf2)
                eqe_vals = arrayfun(@(r) mean(r.SR.EQE(:, ib), 'omitnan'), runs(sr_sf2));
                fprintf(' | SF=1 eta=2 mean EQE=%.2f%%', mean(eqe_vals, 'omitnan'));
            end
            fprintf('\n');
        end
    end

    % Doping profile summary (if available)
    if any(profile_mask)
        fprintf('\n--- Doping Profile Summary ---\n');
        fprintf('Runs with profile data: %d\n', sum(profile_mask));
        n0_vals    = arrayfun(@(r) r.N0_peak_scalar,    runs(profile_mask));
        ef_pk_vals = arrayfun(@(r) r.EF_peak_scalar,    runs(profile_mask));
        ef_sf_vals = arrayfun(@(r) r.EF_surface_scalar, runs(profile_mask));
        fprintf('N0_peak    [cm^-3]: min=%.3e  max=%.3e\n', min(n0_vals),    max(n0_vals));
        fprintf('EF_peak    [eV]:    min=%.4f  max=%.4f\n', min(ef_pk_vals), max(ef_pk_vals));
        fprintf('EF_surface [eV]:    min=%.4f  max=%.4f\n', min(ef_sf_vals), max(ef_sf_vals));
        for ip = 1:length(profile_shapes)
            pname = profile_shapes{ip};
            pm = profile_mask & strcmp(profile_type_col, pname);
            fprintf('  %s: %d runs with profiles\n', pname, sum(pm));
        end
    end
end

fprintf('\n===========================================================\n');
fprintf(' Postprocessing Complete\n');
fprintf('===========================================================\n');

end


%% ========================================================================
%  HELPER FUNCTIONS (postprocess-specific)
%% ========================================================================

function errors = parse_error_log(errorFile)
%PARSE_ERROR_LOG Extract failed run info from errors.txt
%   Updated regex to match run_sweep_chunk_profiles log format:
%   Params: SF=%d, eta=%.1f, jd=%d, sp=%.0e, N0=%.0e, PhiOff=%.4f, prof_idx=%d

    errors = struct('SF_flag', {}, 'eta', {}, 'jd', {}, 'sp', {}, ...
                    'N0', {}, 'Phi_R_offset', {}, 'prof_idx', {}, ...
                    'filename', {}, 'error_msg', {});

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

        % Try new format first (with prof_idx)
        paramMatch = regexp(block, ...
            'Params:\s*SF=(\d+),\s*eta=([\d.]+),\s*jd=(\d+),\s*sp=([\d.e+-]+),\s*N0=([\d.e+-]+),\s*PhiOff=([\d.e+-]+),\s*prof_idx=(\d+)', ...
            'tokens', 'once');
        if isempty(paramMatch)
            % Fall back to old format (without prof_idx)
            paramMatch = regexp(block, ...
                'Params:\s*SF=(\d+),\s*eta=([\d.]+),\s*jd=(\d+),\s*sp=([\d.e+-]+),\s*N_D=([\d.e+-]+),\s*Phi_R_off=([\d.e+-]+)', ...
                'tokens', 'once');
            if isempty(paramMatch), continue; end
            % Append default prof_idx = 0 (unknown)
            paramMatch{end+1} = '0';
        end

        errMatch = regexp(block, 'Error:\s*(.+?)(?:\n|$)', 'tokens', 'once');

        err = struct();
        err.SF_flag      = str2double(paramMatch{1});
        err.eta          = str2double(paramMatch{2});
        err.jd           = str2double(paramMatch{3});
        err.sp           = str2double(paramMatch{4});
        err.N0           = str2double(paramMatch{5});
        err.Phi_R_offset = str2double(paramMatch{6});
        err.prof_idx     = str2double(paramMatch{7});
        err.filename     = fileMatch{1};
        err.error_msg    = '';
        if ~isempty(errMatch)
            err.error_msg = strtrim(errMatch{1});
        end

        errors(end+1) = err; %#ok<AGROW>
    end
end


function idx = find_error(errors, exp, profile_name)
%FIND_ERROR Find matching error entry for expected parameters
%   exp columns: 1:SF 2:eta 3:jd 4:sp 5:N0 6:EF 7:Phi_R 8:offset 9:prof_idx 10:EF_surf

    idx = [];
    tol = 1e-3;

    for i = 1:length(errors)
        if errors(i).SF_flag == exp(1) && ...
           abs(errors(i).eta - exp(2)) < tol && ...
           errors(i).jd == exp(3) && ...
           abs(errors(i).sp - exp(4)) / max(exp(4), 1) < tol && ...
           abs(errors(i).Phi_R_offset - exp(8)) < tol && ...
           errors(i).prof_idx == exp(9)
            idx = i;
            return;
        end
    end

    % Fallback: match by filename (handles edge cases)
    if isempty(idx)
        expected_fn = generate_sweep_filename(exp(1), exp(2), exp(3), exp(4), exp(5), exp(8), profile_name);
        for i = 1:length(errors)
            if strcmp(errors(i).filename, expected_fn)
                idx = i;
                return;
            end
        end
    end
end
