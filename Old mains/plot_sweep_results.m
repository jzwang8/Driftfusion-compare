function plot_sweep_results(summaryFile)
%PLOT_SWEEP_RESULTS Plot EQE and JV curves from sweep results
%
%   plot_sweep_results(summaryFile)
%
%   Configure the filtering options in the CONFIGURATION section below.
%
%   Two modes:
%     1. Run Index Mode: Plot specific runs by their index
%     2. Parameter Mode: Filter by parameter values, vary one parameter
%
%   Examples:
%     plot_sweep_results('sweep_summary_20260113.mat')

close all;

%% ========================================================================
%  CONFIGURATION - Modify this section to control what gets plotted
%% ========================================================================

% ---- Mode 1: Manual run-index mode ----
useRunIndexMode = false;            % Set true to plot specific runs by index
runIdxList      = [50, 150, 250];   % Run indices to plot (from CSV run_number column)

% ---- Mode 2: Parameter-based selection ----
% Choose ONE parameter to report as "varying" in the summary:
%   'SF_flag', 'eta_SF', 'junction_depth', 'SRV', 'EF', 'Phi_R_offset'
% (Coloring is now fixed by SF_flag & eta_SF; varyParam no longer controls colors.)
varyParam = 'eta_SF';

% ---- Plot options ----
plotType       = 'EQE';    % 'EQE', 'JV', or 'both'
lambdaRange    = [400 800]; % Wavelength range for EQE [min max], or [] for full
plotThinCloud  = true;      % Plot thin, light individual curves
plotMeanCurve  = false;     % Plot thick mean curve for each SF/eta category

% ---- Fixed-parameter filters ----
% Set .use = true to lock that parameter to .value
% Set .use = false to allow all values
cfg.SF_flag.use         = false;   cfg.SF_flag.value         = 1;
cfg.eta_SF.use          = true;  cfg.eta_SF.value          = 2;
cfg.junction_depth.use  = false;  cfg.junction_depth.value  = 300;      % nm
cfg.SRV.use             = false;  cfg.SRV.value             = 1000;     % cm/s
cfg.EF.use              = false;  cfg.EF.value              = -4.104;   % eV
cfg.Phi_R_offset.use    = false;  cfg.Phi_R_offset.value    = 0;        % eV

%% ========================================================================
%  END CONFIGURATION
%% ========================================================================

%% Load data
fprintf('Loading %s...\n', summaryFile);
data    = load(summaryFile, 'summary');
summary = data.summary;
runs    = summary.runs;
n_runs  = length(runs);
fprintf('Loaded %d runs\n', n_runs);

%% Select runs based on mode
if useRunIndexMode
    % Mode 1: Use specific run indices
    selected_idx = runIdxList(:);
    selected_idx = selected_idx(selected_idx >= 1 & selected_idx <= n_runs);
    fprintf('Run Index Mode: Selected %d runs\n', length(selected_idx));
else
    % Mode 2: Parameter-based filtering
    mask = true(n_runs, 1);
    
    % Only include successful runs
    for i = 1:n_runs
        if ~strcmp(runs(i).status, 'success')
            mask(i) = false;
        end
    end
    
    % Apply fixed parameter filters
    tol = 1e-3;  % Tolerance for floating point comparison
    
    if cfg.SF_flag.use
        for i = 1:n_runs
            if mask(i) && runs(i).SF_flag ~= cfg.SF_flag.value
                mask(i) = false;
            end
        end
    end
    
    if cfg.eta_SF.use
        for i = 1:n_runs
            if mask(i) && abs(runs(i).eta_SF - cfg.eta_SF.value) > tol
                mask(i) = false;
            end
        end
    end
    
    if cfg.junction_depth.use
        for i = 1:n_runs
            if mask(i) && runs(i).junction_depth_nm ~= cfg.junction_depth.value
                mask(i) = false;
            end
        end
    end
    
    if cfg.SRV.use
        for i = 1:n_runs
            if mask(i) && abs(runs(i).SRV_cm_s - cfg.SRV.value) / cfg.SRV.value > tol
                mask(i) = false;
            end
        end
    end
    
    if cfg.EF.use
        for i = 1:n_runs
            if mask(i) && abs(runs(i).EF_emitter_eV - cfg.EF.value) > 0.01
                mask(i) = false;
            end
        end
    end
    
    if cfg.Phi_R_offset.use
        for i = 1:n_runs
            if mask(i) && abs(runs(i).Phi_R_offset_eV - cfg.Phi_R_offset.value) > tol
                mask(i) = false;
            end
        end
    end
    
    selected_idx = find(mask);
    fprintf('Parameter Mode: Selected %d runs (summary varying %s)\n', ...
            length(selected_idx), varyParam);
end

if isempty(selected_idx)
    warning('No runs match the specified filters');
    return;
end

%% Color assignment based on (SF_flag, eta_SF)
[colors, cat_idx, legend_labels] = assign_sf_eta_colors(runs, selected_idx);

%% Plot
if strcmp(plotType, 'EQE') || strcmp(plotType, 'both')
    plot_EQE_curves(runs, selected_idx, colors, cat_idx, legend_labels, ...
                    lambdaRange, plotThinCloud, plotMeanCurve);
end

if strcmp(plotType, 'JV') || strcmp(plotType, 'both')
    plot_JV_curves(runs, selected_idx, colors, cat_idx, legend_labels, ...
                   plotThinCloud, plotMeanCurve);
end

%% Print summary
print_summary(runs, selected_idx, varyParam);

end


%% ========================================================================
%  HELPER FUNCTIONS
%% ========================================================================

function [colors, cat_idx, legend_labels] = assign_sf_eta_colors(runs, selected_idx)
%ASSIGN_SF_ETA_COLORS
%   Categories & fixed colors:
%     1: SF=0                    → blue
%     2: SF=1, eta_SF=0          → red
%     3: SF=1, eta_SF=1          → yellow
%     4: SF=1, eta_SF=2          → green

    n       = length(selected_idx);
    cat_idx = zeros(n,1);
    colors  = zeros(n,3);
    
    % Base colors (RGB)
    base_colors = [ ...
        0.0  0.0  1.0;    % 1: blue
        1.0  0.0  0.0;    % 2: red
        1.0  0.85 0.0;    % 3: yellow-ish
        0.0  0.6  0.0];   % 4: green
    
    % Legend labels in category order
    legend_labels = { ...
        'SF = 0', ...
        'SF = 1, \eta_{SF} = 0', ...
        'SF = 1, \eta_{SF} = 1', ...
        'SF = 1, \eta_{SF} = 2'};
    
    for i = 1:n
        r = runs(selected_idx(i));
        
        if r.SF_flag == 0
            cat_idx(i) = 1;
        elseif r.SF_flag == 1 && r.eta_SF == 0
            cat_idx(i) = 2;
        elseif r.SF_flag == 1 && r.eta_SF == 1
            cat_idx(i) = 3;
        elseif r.SF_flag == 1 && r.eta_SF == 2
            cat_idx(i) = 4;
        else
            cat_idx(i) = 0; % unexpected combination
        end
        
        if cat_idx(i) > 0
            colors(i,:) = base_colors(cat_idx(i),:);
        else
            colors(i,:) = [0 0 0]; % fallback (black)
        end
    end
end


function plot_EQE_curves(runs, selected_idx, colors, cat_idx, legend_labels, ...
                         lambdaRange, plotThinCloud, plotMeanCurve)
%PLOT_EQE_CURVES Plot EQE curves with SF/eta-based colors

    % Styling constants
    AXIS_FS    = 16;
    LABEL_FS   = 18;
    LEGEND_FS  = 14;
    LW_THIN    = 0.6;
    LW_MEAN    = 2.5;
    ALPHA_THIN = 0.5;   % "opacity" of line color vs white

    AXIS_BG = [1 1 1];   % assume white background
    
    figure('Name', 'EQE Comparison', 'Position', [100 100 800 600], 'Color', 'w');
    ax = axes('Parent', gcf);
    hold(ax, 'on');
    
    n        = length(selected_idx);
    cat_used = cat_idx;
    
    % Plot individual curves (lightened color for "semi-transparency")
    for i = 1:n
        r = runs(selected_idx(i));
        if ~isfield(r, 'EQE') || isempty(r.EQE.wavelengths)
            continue;
        end
        
        wl  = r.EQE.wavelengths;
        eqe = r.EQE.EQE;
        
        if plotThinCloud
            base_col  = colors(i,:);
            faded_col = ALPHA_THIN*base_col + (1-ALPHA_THIN)*AXIS_BG;
            plot(ax, wl, eqe, '-', 'Color', faded_col, 'LineWidth', LW_THIN);
        end
    end
    
    % Mean curves per category (if requested)
    base_colors = [ ...
        0.0  0.0  1.0;    % blue
        1.0  0.0  0.0;    % red
        1.0  0.85 0.0;    % yellow-ish
        0.0  0.6  0.0];   % green
    
    h_leg       = gobjects(4,1);
    use_legend  = false(4,1);
    
    if plotMeanCurve
        % Determine common wavelength grid from first valid run
        wl_ref = [];
        for i = 1:n
            r = runs(selected_idx(i));
            if isfield(r, 'EQE') && ~isempty(r.EQE.wavelengths)
                wl_ref = r.EQE.wavelengths(:).';
                break;
            end
        end
        
        if ~isempty(wl_ref)
            eqe_group = cell(4,1);
            for k = 1:4
                eqe_group{k} = [];
            end
            
            for i = 1:n
                if cat_used(i) <= 0, continue; end
                r = runs(selected_idx(i));
                if ~isfield(r, 'EQE') || isempty(r.EQE.wavelengths)
                    continue;
                end
                wl  = r.EQE.wavelengths(:).';
                eqe = r.EQE.EQE(:).';
                if ~isequal(wl, wl_ref)
                    % Interpolate to common grid
                    eqe = interp1(wl, eqe, wl_ref, 'linear', NaN);
                end
                eqe_group{cat_used(i)} = [eqe_group{cat_used(i)}; eqe];
            end
            
            % Plot mean per category
            for k = 1:4
                if ~isempty(eqe_group{k})
                    mean_eqe = mean(eqe_group{k}, 1, 'omitnan');
                    h_leg(k) = plot(ax, wl_ref, mean_eqe, '-', ...
                        'Color', base_colors(k,:), 'LineWidth', LW_MEAN);
                    use_legend(k) = true;
                end
            end
        end
    end
    
    % If no mean curves, still create legend entries for categories present
    if ~plotMeanCurve
        present_cats = unique(cat_used(cat_used > 0));
        for k = present_cats(:).'
            h_leg(k)     = plot(ax, NaN, NaN, '-', ...
                                'Color', base_colors(k,:), 'LineWidth', LW_MEAN);
            use_legend(k) = true;
        end
    end
    
    % Labels, limits, formatting
    xlabel(ax, 'Wavelength (nm)', 'FontSize', LABEL_FS);
    ylabel(ax, 'EQE (%)',        'FontSize', LABEL_FS);
    
    if ~isempty(lambdaRange)
        xlim(ax, lambdaRange);
    end
    ylim(ax, [10 110]);
    
    set(ax, 'FontSize', AXIS_FS, 'LineWidth', 1.2, 'Box', 'on');
    grid(ax, 'on');
    
    % Legend with white background
    if any(use_legend)
        labels_used = legend_labels(use_legend);
        h_used      = h_leg(use_legend);
        lgd = legend(ax, h_used, labels_used, ...
                     'Location', 'southeast', 'FontSize', LEGEND_FS);
        set(lgd, 'Color', 'w', 'Box', 'on');  % white legend background
    end
    
    hold(ax, 'off');
end


function plot_JV_curves(runs, selected_idx, colors, cat_idx, legend_labels, ...
                        plotThinCloud, plotMeanCurve)
%PLOT_JV_CURVES Plot J-V curves with SF/eta-based colors

    % Styling constants
    AXIS_FS    = 16;
    LABEL_FS   = 18;
    LEGEND_FS  = 14;
    LW_THIN    = 0.6;
    LW_MEAN    = 2.5;
    ALPHA_THIN = 0.5;
    
    AXIS_BG = [1 1 1];
    
    figure('Name', 'J-V Comparison', 'Position', [150 150 800 600], 'Color', 'w');
    ax = axes('Parent', gcf);
    hold(ax, 'on');
    
    n        = length(selected_idx);
    cat_used = cat_idx;
    
    % Plot individual curves
    for i = 1:n
        r = runs(selected_idx(i));
        if ~isfield(r, 'JV') || isempty(r.JV.Vapp)
            continue;
        end
        
        % Forward sweep only
        V = r.JV.Vapp;
        J = r.JV.Jtot;
        [~, idx_max] = max(V);
        V = V(1:idx_max);
        J = J(1:idx_max);
        
        if plotThinCloud
            base_col  = colors(i,:);
            faded_col = ALPHA_THIN*base_col + (1-ALPHA_THIN)*AXIS_BG;
            plot(ax, V, J, '-', 'Color', faded_col, 'LineWidth', LW_THIN);
        end
    end
    
    % Base colors for categories
    base_colors = [ ...
        0.0  0.0  1.0;    % blue
        1.0  0.0  0.0;    % red
        1.0  0.85 0.0;    % yellow-ish
        0.0  0.6  0.0];   % green
    
    h_leg      = gobjects(4,1);
    use_legend = false(4,1);
    
    % Mean curves (optional)
    if plotMeanCurve
        V_ref = [];
        JV_group = cell(4,1);
        for k = 1:4
            JV_group{k} = [];
        end
        
        for i = 1:n
            if cat_used(i) <= 0, continue; end
            r = runs(selected_idx(i));
            if ~isfield(r, 'JV') || isempty(r.JV.Vapp)
                continue;
            end
            V = r.JV.Vapp;
            J = r.JV.Jtot;
            [~, idx_max] = max(V);
            V = V(1:idx_max);
            J = J(1:idx_max);
            
            if isempty(V_ref)
                V_ref = V(:).';
            end
            J_interp = interp1(V(:).', J(:).', V_ref, 'linear', NaN);
            JV_group{cat_used(i)} = [JV_group{cat_used(i)}; J_interp];
        end
        
        for k = 1:4
            if ~isempty(JV_group{k})
                mean_J = mean(JV_group{k}, 1, 'omitnan');
                h_leg(k) = plot(ax, V_ref, mean_J, '-', ...
                    'Color', base_colors(k,:), 'LineWidth', LW_MEAN);
                use_legend(k) = true;
            end
        end
    end
    
    % If no mean curves, create dummy legend lines for used categories
    if ~plotMeanCurve
        present_cats = unique(cat_used(cat_used > 0));
        for k = present_cats(:).'
            h_leg(k)     = plot(ax, NaN, NaN, '-', ...
                                'Color', base_colors(k,:), 'LineWidth', LW_MEAN);
            use_legend(k) = true;
        end
    end
    
    % Zero-current line
    x_limits = xlim(ax);
    plot(ax, x_limits, [0 0], 'k--', 'LineWidth', 0.8);
    
    % Labels / formatting
    xlabel(ax, 'Voltage (V)',               'FontSize', LABEL_FS);
    ylabel(ax, 'Current density (mA/cm^2)', 'FontSize', LABEL_FS);
    
    set(ax, 'FontSize', AXIS_FS, 'LineWidth', 1.2, 'Box', 'on');
    grid(ax, 'on');
    
    % Legend with white background
    if any(use_legend)
        labels_used = legend_labels(use_legend);
        h_used      = h_leg(use_legend);
        lgd = legend(ax, h_used, labels_used, ...
                     'Location', 'best', 'FontSize', LEGEND_FS);
        set(lgd, 'Color', 'w', 'Box', 'on');
    end
    
    hold(ax, 'off');
end


function print_summary(runs, selected_idx, varyParam)
%PRINT_SUMMARY Print summary table of selected runs

    fprintf('\n');
    fprintf('===============================================================\n');
    fprintf(' Selected Runs Summary (varying %s)\n', varyParam);
    fprintf('===============================================================\n');
    fprintf('%-5s %-4s %-5s %-4s %-9s %-8s %-9s %-9s %-8s %-8s\n', ...
            'Run', 'SF', 'eta', 'jd', 'SRV', 'EF', 'PhiOff', 'PCE%', 'Voc', 'Jsc');
    fprintf('---------------------------------------------------------------\n');
    
    n = length(selected_idx);
    for i = 1:min(30, n)
        r = runs(selected_idx(i));
        fprintf('%-5d %-4d %-5.1f %-4d %-9.0e %-8.3f %-9.4f %-9.4f %-8.4f %-8.4f\n', ...
                r.run_number, r.SF_flag, r.eta_SF, r.junction_depth_nm, ...
                r.SRV_cm_s, r.EF_emitter_eV, r.Phi_R_offset_eV, ...
                r.PCE_pct, r.Voc_V, r.Jsc_mA_cm2);
    end
    
    if n > 30
        fprintf('... and %d more runs\n', n - 30);
    end
    fprintf('===============================================================\n');
end