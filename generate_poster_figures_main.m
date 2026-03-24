function generate_poster_figures_main(summaryFile)
%GENERATE_POSTER_FIGURES Create publication-quality figures for research poster
%
%   generate_poster_figures('sweep_summary_20260113.mat')
%
%   Generates selected figures based on the FIGURES_TO_PLOT configuration below.
%
%   Available figures:
%     1  - EQE comparison showing >100% region (SF signature)
%     2  - Heatmap: PCE vs (SRV, E_F)
%     3  - Bar chart: PCE gain from SF at different SRV
%     4  - Line plot: PCE vs junction depth
%     5  - Bar chart: Integrated blue EQE (405-500nm)
%     6  - 3-panel EQE: effect of eta, SRV, Phi_R
%     7  - 3D surface: PCE vs (SRV, junction depth)
%     8  - Box plot: PCE distribution by SF/eta
%     9  - Correlation matrix: parameters vs metrics
%     10 - Parallel coordinates: top vs bottom performers
%     11 - Scatter matrix: PCE vs key parameters
%     12 - Heatmap: Delta PCE (SF improvement)
%     13 - Radar chart: best vs worst config
%     14 - Stacked area: EQE decomposition
%     15 - Summary table
%     16 - Best device EQE
%     17 - Best device J-V
%     18 - Best device vs baseline comparison

close all;

%% ========================================================================
%  CONFIGURATION - Select which figures to generate
%% ========================================================================

% Set to 'all' to plot everything, or specify figure numbers as array
% Examples:
%   FIGURES_TO_PLOT = 'all';           % Plot all 18 figures
%   FIGURES_TO_PLOT = [1, 16, 17, 18]; % Plot only EQE comparison and best device
%   FIGURES_TO_PLOT = [1:6];           % Plot figures 1-6
%   FIGURES_TO_PLOT = [16, 17, 18];    % Plot only best device figures

FIGURES_TO_PLOT = [6];
% FIGURES_TO_PLOT = 'all';

%% ========================================================================
%  END CONFIGURATION
%% ========================================================================

% Process figure selection
if ischar(FIGURES_TO_PLOT) && strcmp(FIGURES_TO_PLOT, 'all')
    figures_to_plot = 1:18;
else
    figures_to_plot = FIGURES_TO_PLOT;
end

fprintf('Will generate figures: %s\n', mat2str(figures_to_plot));

%% ========================================================================
%  COLORBLIND-FRIENDLY PALETTE (used throughout)
%% ========================================================================
% Two-class: baseline vs SF
CB_BLUE   = [0.12 0.47 0.71];   % baseline / "good" end
CB_ORANGE = [0.85 0.37 0.01];   % SF / "high" end

% 5-level diverging ramp  (blue -> grey -> orange)
CB_RAMP5 = [ ...
    0.12 0.47 0.71;   % level 1 (best / lowest)
    0.50 0.70 0.86;   % level 2
    0.60 0.60 0.60;   % level 3 (neutral)
    0.93 0.60 0.25;   % level 4
    0.85 0.37 0.01];  % level 5 (worst / highest)

% 4-level categorical (for box-plot scatter etc.)
CB_CAT4 = [ ...
    0.50 0.50 0.50;   % baseline (grey)
    0.85 0.37 0.01;   % SF eta=0 (orange)
    0.55 0.34 0.64;   % SF eta=1 (purple)
    0.12 0.47 0.71];  % SF eta=2 (blue)

%% Load data
fprintf('Loading %s...\n', summaryFile);
data = load(summaryFile, 'summary');
summary = data.summary;
runs = summary.runs;
T = summary.table;
fprintf('Loaded %d runs (%d successful)\n', length(runs), summary.n_success);

% Filter to successful runs only
success_mask = strcmp(T.status, 'success');
T_success = T(success_mask, :);

% ------------------------------------------------------------------------
% Shared sweep settings used across multiple figures
% ------------------------------------------------------------------------
optimal_SRV   = 10000;        % Best surface passivation
optimal_EF    = -4.104;    % ~1e19 doping
optimal_PhiOff = 0;
optimal_jd    = 200;
tol           = 0.01;

% SRV levels used in multiple figures
SRV_list = [10, 100, 1000, 10000, 100000];

%% ========================================================================
%  FIGURE 1: EQE Comparison - SF Signature (EQE > 100%)
%% ========================================================================
if ismember(1, figures_to_plot)
figure('Name', 'EQE_SF_Comparison', 'Position', [50 50 900 600], 'Color', 'w');

% Find matching runs
for eta_val = [0, 2]  % Compare eta=0 (no SF benefit) vs eta=2 (ideal SF)
    mask = (T_success.SF_flag == 1) & ...
           (abs(T_success.eta_SF - eta_val) < tol) & ...
           (T_success.SRV_cm_s == optimal_SRV) & ...
           (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
           (abs(T_success.Phi_R_offset_eV - optimal_PhiOff) < tol) & ...
           (T_success.junction_depth_nm == optimal_jd);
    
    idx = find(mask);
    if ~isempty(idx)
        run_idx = T_success.run_number(idx(1));
        r = runs(run_idx);
        
        if eta_val == 0
            plot(r.EQE.wavelengths, r.EQE.EQE, '-', 'LineWidth', 2.5, ...
                 'Color', CB_ORANGE, ...
                 'DisplayName', 'SF layer, \eta_{SF} = 0');
        else
            plot(r.EQE.wavelengths, r.EQE.EQE, '-', 'LineWidth', 2.5, ...
                 'Color', CB_BLUE, ...
                 'DisplayName', 'SF layer, \eta_{SF} = 2');
        end
        hold on;
    end
end

% Also plot baseline (SF=0)
mask = (T_success.SF_flag == 0) & ...
       (T_success.SRV_cm_s == optimal_SRV) & ...
       (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
       (abs(T_success.Phi_R_offset_eV - optimal_PhiOff) < tol) & ...
       (T_success.junction_depth_nm == optimal_jd);
idx = find(mask);
if ~isempty(idx)
    run_idx = T_success.run_number(idx(1));
    r = runs(run_idx);
    plot(r.EQE.wavelengths, r.EQE.EQE, '-', 'LineWidth', 2.5, ...
         'Color', [0.50 0.50 0.50], ...
         'DisplayName', 'No SF layer (baseline)');
end

% Add 100% threshold line
xlims = [400 650];
plot(xlims, [100 100], 'k--', 'LineWidth', 1.5, 'DisplayName', 'EQE = 100%');

% Shade the >100% region
fill([xlims(1) xlims(2) xlims(2) xlims(1)], [100 100 150 150], ...
     [0.85 0.92 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.3, ...
     'HandleVisibility', 'off');

xlabel('Wavelength (nm)', 'FontSize', 18);
ylabel('EQE (%)', 'FontSize', 18);
title('Singlet Fission Enhancement of EQE', 'FontSize', 20);
xlim(xlims);
ylim([0 150]);
set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'on');
grid on;
legend('Location', 'northeast', 'FontSize', 12);

% Add annotation
text(420, 130, 'SF benefit region', 'FontSize', 12, 'Color', CB_BLUE);

hold off;
end

%% ========================================================================
%  FIGURE 2: Heatmap - PCE vs (SRV, N_D) for ideal SF
%% ========================================================================
if ismember(2, figures_to_plot)
figure('Name', 'Heatmap_PCE_SRV_ND', 'Position', [100 100 600 450], 'Color', 'w');

% Filter for SF=1, eta=2, fixed jd and Phi_R_offset
mask = (T_success.SF_flag == 1) & ...
       (abs(T_success.eta_SF - 2) < tol) & ...
       (T_success.junction_depth_nm == 300) & ...
       (abs(T_success.Phi_R_offset_eV - 0) < tol);

T_sub = T_success(mask, :);

% Get unique values (still in EF space)
SRV_vals = unique(T_sub.SRV_cm_s);
EF_vals  = unique(T_sub.EF_emitter_eV);

% --- Map E_F values to donor concentration N_D for axis labels ---
EF_map = [-4.05, -4.10, -4.22, -4.34];
ND_map = [8.06e19, 1e19, 1e17, 1e15];

ND_axis = NaN(size(EF_vals));
for ii = 1:numel(EF_vals)
    [~, idx_map] = min(abs(EF_map - EF_vals(ii)));
    ND_axis(ii) = ND_map(idx_map);
end

% Build PCE matrix
PCE_matrix = NaN(length(EF_vals), length(SRV_vals));
for i = 1:length(EF_vals)
    for j = 1:length(SRV_vals)
        idx = find(abs(T_sub.EF_emitter_eV - EF_vals(i)) < tol & ...
                   abs(T_sub.SRV_cm_s - SRV_vals(j))/SRV_vals(j) < tol, 1);
        if ~isempty(idx)
            PCE_matrix(i, j) = T_sub.PCE_pct(idx);
        end
    end
end

% Plot heatmap
h = imagesc(1:length(SRV_vals), 1:length(EF_vals), PCE_matrix);
colormap(blueorange(256));
set(h, 'AlphaData', ~isnan(PCE_matrix));
set(gca, 'Color', 'w');

cb = colorbar;
ylabel(cb, 'PCE (%)', 'FontSize', 18);
set(gca, 'YDir', 'normal');

% Axis tick labels: SRV on x, N_D on y
set(gca, 'XTick', 1:length(SRV_vals), ...
         'XTickLabel', arrayfun(@(x) sprintf('%.0e', x), SRV_vals, 'UniformOutput', false));

set(gca, 'YTick', 1:length(EF_vals), ...
         'YTickLabel', arrayfun(@(x) sprintf('%.1e', x), ND_axis, 'UniformOutput', false));

xlabel('Surface recombination velocity s_p (cm/s)', 'FontSize', 20);
ylabel('n-emitter doping concentration N_D (cm^{-3})', 'FontSize', 20);
set(gca, 'FontSize', 16, 'LineWidth', 1.2);

% Add value annotations
for i = 1:length(EF_vals)
    for j = 1:length(SRV_vals)
        if ~isnan(PCE_matrix(i,j))
            text(j, i, sprintf('%.2f', PCE_matrix(i,j)), ...
                 'HorizontalAlignment', 'center', 'FontSize', 16);
        end
    end
end

end

%% ========================================================================
%  FIGURE 3: Line Plot - SF Benefit by Surface Recombination Velocity
%% ========================================================================
if ismember(3, figures_to_plot)
figure('Name', 'SF_Benefit_vs_SRV', 'Position', [150 150 450 500], 'Color', 'w');

% Compare PCE for matched configurations: SF=0 vs SF=1 (eta=2)
nSRV = numel(SRV_list);
PCE_baseline = zeros(nSRV, 1);
PCE_SF       = zeros(nSRV, 1);

for k = 1:nSRV
    % Baseline (SF=0)
    mask = (T_success.SF_flag == 0) & ...
           (T_success.SRV_cm_s == SRV_list(k)) & ...
           (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
           (T_success.junction_depth_nm == 300) & ...
           (abs(T_success.Phi_R_offset_eV - 0) < tol);
    if any(mask)
        PCE_baseline(k) = mean(T_success.PCE_pct(mask), 'omitnan');
    else
        PCE_baseline(k) = NaN;
    end
    
    % With SF (eta=2)
    mask = (T_success.SF_flag == 1) & ...
           (abs(T_success.eta_SF - 2) < tol) & ...
           (T_success.SRV_cm_s == SRV_list(k)) & ...
           (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
           (T_success.junction_depth_nm == 300) & ...
           (abs(T_success.Phi_R_offset_eV - 0) < tol);
    if any(mask)
        PCE_SF(k) = mean(T_success.PCE_pct(mask), 'omitnan');
    else
        PCE_SF(k) = NaN;
    end
end

% ---- Line plot (log-x for SRV) ----
hold on;

h1 = semilogx(SRV_list, PCE_baseline, '-o', ...
              'LineWidth', 2.5, ...
              'MarkerSize', 10, ...
              'MarkerFaceColor', [0.50 0.50 0.50], ...
              'Color', [0.50 0.50 0.50], ...
              'DisplayName', 'Baseline');

h2 = semilogx(SRV_list, PCE_SF, '-s', ...
              'LineWidth', 2.5, ...
              'MarkerSize', 10, ...
              'MarkerFaceColor', CB_BLUE, ...
              'Color', CB_BLUE, ...
              'DisplayName', 'SF (\eta=2)');

ylabel('PCE (%)', 'FontSize', 20);
xlabel('Surface recombination velocity s_p (cm/s)', 'FontSize', 20);

legend([h1 h2], {'Baseline', 'SF (\eta=2)'}, ...
       'Location', 'northeast', 'FontSize', 16);

% Log x-axis, ticks at your SRV values only
set(gca, 'XScale', 'log', 'XTick', SRV_list);

% Nicer tick labels as 10^n
set(gca, 'TickLabelInterpreter', 'tex');
xticklabels(arrayfun(@(x) sprintf('10^{%d}', round(log10(x))), ...
                     SRV_list, 'UniformOutput', false));

% Remove minor ticks / minor grid on x-axis
ax = gca;
ax.XMinorTick = 'off';
ax.XMinorGrid = 'off';

grid on;

xlim([min(SRV_list)*0.8, max(SRV_list)*1.2]);

% Axes styling
set(gca, 'Box', 'on', ...
         'LineWidth', 1.5, ...
         'FontSize', 18);

% ---- Add value labels near markers ----
ymin = nanmin([PCE_baseline(:); PCE_SF(:)]);
ymax = nanmax([PCE_baseline(:); PCE_SF(:)]);
dy   = 0.02*(ymax - ymin);

for k = 1:nSRV
    if ~isnan(PCE_baseline(k))
        text(SRV_list(k), PCE_baseline(k) + dy, ...
             sprintf('%.2f%%', PCE_baseline(k)), ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment',   'bottom', ...
             'FontSize', 16, ...
             'Color', [0.30 0.30 0.30]);
    end

    if ~isnan(PCE_SF(k))
        text(SRV_list(k), PCE_SF(k) + 2*dy, ...
             sprintf('%.2f%%', PCE_SF(k)), ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment',   'bottom', ...
             'FontSize', 16, ...
             'Color', CB_BLUE);
    end
end

hold off;
end

%% ========================================================================
%  FIGURE 4: Line Plot - PCE vs Junction Depth for different eta
%% ========================================================================
if ismember(4, figures_to_plot)
figure('Name', 'PCE_vs_JD', 'Position', [200 200 600 400], 'Color', 'w');

jd_vals   = [200, 250, 300, 350, 400];
eta_vals  = [2];
markers    = {'o'};

hold on;

% ---------- SF (eta = 2) ----------
PCE_vs_jd = NaN(length(jd_vals), 1);
for j = 1:length(jd_vals)
    mask = (T_success.SF_flag == 1) & ...
           (abs(T_success.eta_SF - eta_vals(1)) < tol) & ...
           (T_success.junction_depth_nm == jd_vals(j)) & ...
           (T_success.SRV_cm_s == optimal_SRV) & ...
           (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
           (abs(T_success.Phi_R_offset_eV - 0) < tol);
    if any(mask)
        PCE_vs_jd(j) = mean(T_success.PCE_pct(mask), 'omitnan');
    end
end

h_sf = plot(jd_vals, PCE_vs_jd, ['-' markers{1}], ...
        'Color',       CB_BLUE, ...
        'LineWidth',   2, ...
        'MarkerSize',  8, ...
        'MarkerFaceColor', CB_BLUE, ...
        'DisplayName', '\eta_{SF} = 2');

% ---------- Baseline (SF = 0) ----------
PCE_baseline_jd = NaN(length(jd_vals), 1);
for j = 1:length(jd_vals)
    mask = (T_success.SF_flag == 0) & ...
           (T_success.junction_depth_nm == jd_vals(j)) & ...
           (T_success.SRV_cm_s == optimal_SRV) & ...
           (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
           (abs(T_success.Phi_R_offset_eV - 0) < tol);
    if any(mask)
        PCE_baseline_jd(j) = mean(T_success.PCE_pct(mask), 'omitnan');
    end
end

h_base = plot(jd_vals, PCE_baseline_jd, '--^', ...
           'Color', [0.50 0.50 0.50], ...
           'LineWidth', 2, ...
           'MarkerSize', 8, ...
           'MarkerFaceColor', [0.50 0.50 0.50], ...
           'DisplayName', 'Baseline (no SF)');

% ---------- Linear fits & slope comparison ----------
valid_sf    = ~isnan(PCE_vs_jd);
valid_base  = ~isnan(PCE_baseline_jd);

p_sf   = polyfit(jd_vals(valid_sf),   PCE_vs_jd(valid_sf),       1);
p_base = polyfit(jd_vals(valid_base), PCE_baseline_jd(valid_base), 1);

slope_sf   = p_sf(1);
slope_base = p_base(1);

fprintf('FIG 4 slopes: baseline = %.3g %%/nm,  eta=2 = %.3g %%/nm\n', ...
        slope_base, slope_sf);

x_txt = 205;
y_txt = min([PCE_vs_jd; PCE_baseline_jd]) + 0.5 * ...
        (max([PCE_vs_jd; PCE_baseline_jd]) - min([PCE_vs_jd; PCE_baseline_jd]));
txt_str = sprintf(['Slope (per 100nm):\\newline', ...
                   'Baseline: %.4f%% / 100nm\\newline', ...
                   '\\eta_{SF}=2: %.4f%% / 100nm'], ...
                   100*slope_base, 100*slope_sf);
text(x_txt, y_txt, txt_str, ...
     'FontSize', 14, ...
     'BackgroundColor', 'w', ...
     'EdgeColor', [0.5 0.5 0.5], ...
     'Margin', 6);

hold off;

xlabel('Junction depth (nm)', 'FontSize', 20);
ylabel('PCE (%)', 'FontSize', 20);
legend([h_sf, h_base], 'Location', 'best', 'FontSize', 16);
set(gca, 'FontSize', 16, 'LineWidth', 1.2, 'Box', 'on');
xlim([180 420]);
grid on;

end

%% ========================================================================
%  FIGURE 4: Line Plot - EQE(520 nm) vs Junction Depth for baseline and η=2
%% ========================================================================
if ismember(4, figures_to_plot)
figure('Name', 'EQE520_vs_JD', 'Position', [200 200 600 400], 'Color', 'w');

jd_vals   = [200, 250, 300, 350, 400];
eta_val   = 2;
lambda_pk = 520;

EQE520_sf    = NaN(length(jd_vals), 1);
EQE520_base  = NaN(length(jd_vals), 1);

get_eqe_at_lambda = @(run_idx_list, lambda_target) ...
    local_mean_eqe_at_lambda(runs, run_idx_list, lambda_target);

for j = 1:length(jd_vals)
    jd = jd_vals(j);

    mask_base = (T_success.SF_flag == 0) & ...
                (T_success.junction_depth_nm == jd) & ...
                (T_success.SRV_cm_s == optimal_SRV) & ...
                (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
                (abs(T_success.Phi_R_offset_eV - 0) < tol);

    if any(mask_base)
        run_idx_list = T_success.run_number(mask_base);
        EQE520_base(j) = get_eqe_at_lambda(run_idx_list, lambda_pk);
    end

    mask_sf = (T_success.SF_flag == 1) & ...
              (abs(T_success.eta_SF - eta_val) < tol) & ...
              (T_success.junction_depth_nm == jd) & ...
              (T_success.SRV_cm_s == optimal_SRV) & ...
              (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
              (abs(T_success.Phi_R_offset_eV - 0) < tol);

    if any(mask_sf)
        run_idx_list = T_success.run_number(mask_sf);
        EQE520_sf(j) = get_eqe_at_lambda(run_idx_list, lambda_pk);
    end
end

hold on;

h_sf = plot(jd_vals, EQE520_sf, '-o', ...
    'Color', CB_BLUE, ...
    'LineWidth', 2, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', CB_BLUE, ...
    'DisplayName', '\eta_{SF} = 2');

h_base = plot(jd_vals, EQE520_base, '--^', ...
    'Color', [0.50 0.50 0.50], ...
    'LineWidth', 2, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', [0.50 0.50 0.50], ...
    'DisplayName', 'Baseline (no SF)');

xlabel('Junction depth (nm)', 'FontSize', 20);
ylabel('EQE at 520nm (%)', 'FontSize', 20);
legend([h_sf, h_base], 'Location', 'best', 'FontSize', 16);

set(gca, 'FontSize', 16, 'LineWidth', 1.2, 'Box', 'on');
xlim([180 420]);
grid on;

all_vals = [EQE520_sf(:); EQE520_base(:)];
all_vals = all_vals(~isnan(all_vals));
if ~isempty(all_vals)
    y_min = min(all_vals);
    y_max = max(all_vals);
    yrng  = y_max - y_min;
    if yrng == 0
        yrng = max(1, 0.1*y_max);
    end
    ylim([y_min - 0.15*yrng, y_max + 0.15*yrng]);
end

valid_sf   = ~isnan(EQE520_sf);
valid_base = ~isnan(EQE520_base);

if sum(valid_sf) >= 2 && sum(valid_base) >= 2
    p_sf   = polyfit(jd_vals(valid_sf),   EQE520_sf(valid_sf),      1);
    p_base = polyfit(jd_vals(valid_base), EQE520_base(valid_base),  1);

    slope_sf   = p_sf(1);
    slope_base = p_base(1);

    fprintf('FIG 4 (EQE 520 nm) slopes: baseline = %.3g %%/nm,  eta=2 = %.3g %%/nm\n', ...
            slope_base, slope_sf);

    x_txt = 205;
    y_txt = y_min + 0.5 * (y_max - y_min);
    txt_str = sprintf(['Slope:\\newline', ...
                       '\\eta_{SF}=2: %.2f%% / 100nm\\newline',...
                        'Baseline: %.2f%% / 100nm'], ...
                       100*slope_sf, 100*slope_base);

    text(x_txt, y_txt, txt_str, ...
         'FontSize', 14, ...
         'BackgroundColor', 'w', ...
         'EdgeColor', [0.5 0.5 0.5], ...
         'Margin', 6);
end

hold off;
end

%% ========================================================================
%  FIGURE 5: Integrated EQE in Blue Region (405-500nm) - SF Signature
%% ========================================================================
if ismember(5, figures_to_plot)
figure('Name', 'Blue_EQE_Line', 'Position', [250 250 700 500], 'Color', 'w');

eta_vals   = [0, 0.5, 1, 1.5, 2];
n_eta      = numel(eta_vals);

blue_EQE = zeros(1 + n_eta, 1);

r_sample = runs(find(strcmp({runs.status}, 'success'), 1));
wl = r_sample.EQE.wavelengths;
blue_mask_wl = (wl >= 405) & (wl <= 500);

param_mask_base = @(Ttbl) (Ttbl.SF_flag == 0) & ...
                          (Ttbl.SRV_cm_s == optimal_SRV) & ...
                          (abs(Ttbl.EF_emitter_eV - optimal_EF) < tol) & ...
                          (Ttbl.junction_depth_nm == optimal_jd) & ...
                          (abs(Ttbl.Phi_R_offset_eV - 0) < tol);

mask = param_mask_base(T_success);
if any(mask)
    run_idx = T_success.run_number(find(mask, 1));
    blue_EQE(1) = mean(runs(run_idx).EQE.EQE(blue_mask_wl), 'omitnan');
else
    blue_EQE(1) = NaN;
end

for e = 1:n_eta
    eta_val = eta_vals(e);
    mask = (T_success.SF_flag == 1) & ...
           (abs(T_success.eta_SF - eta_val) < tol) & ...
           (T_success.SRV_cm_s == optimal_SRV) & ...
           (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
           (T_success.junction_depth_nm == optimal_jd) & ...
           (abs(T_success.Phi_R_offset_eV - 0) < tol);
    if any(mask)
        run_idx = T_success.run_number(find(mask, 1));
        blue_EQE(e+1) = mean(runs(run_idx).EQE.EQE(blue_mask_wl), 'omitnan');
    else
        blue_EQE(e+1) = NaN;
    end
end

hold on;

if ~isnan(blue_EQE(1))
    h_base = yline(blue_EQE(1), '--', 'Baseline', ...
                   'LineWidth', 2, 'Color', [0.50 0.50 0.50]);
else
    h_base = [];
end

sf_vals = blue_EQE(2:end);
valid_sf = ~isnan(sf_vals);

h_sf = plot(eta_vals(valid_sf), sf_vals(valid_sf), '-o', ...
            'LineWidth', 2.5, ...
            'MarkerSize', 8, ...
            'Color', CB_BLUE, ...
            'MarkerFaceColor', [0.50 0.70 0.86], ...
            'DisplayName', 'SF device');

h_100 = yline(100, 'k--', 'EQE = 100%', 'LineWidth', 1.5);

xlabel('\eta_{SF}', 'FontSize', 16);
ylabel('Mean EQE in Blue Region (405-500nm) (%)', 'FontSize', 16);
title('Blue-Region EQE Enhancement vs \eta_{SF}', 'FontSize', 18);

xlim([min(eta_vals) - 0.1, max(eta_vals) + 0.1]);

valid_vals = blue_EQE(~isnan(blue_EQE));
if ~isempty(valid_vals)
    ylim([0, max(valid_vals)*1.15]);
end

set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'on');
grid on;

if ~isempty(h_base)
    legend([h_sf, h_base, h_100], ...
           {'SF device', 'Baseline', 'EQE = 100%'}, ...
           'Location', 'best', 'FontSize', 12);
else
    legend([h_sf, h_100], ...
           {'SF device', 'EQE = 100%'}, ...
           'Location', 'best', 'FontSize', 12);
end

for k = 1:n_eta
    if ~isnan(sf_vals(k))
        text(eta_vals(k), sf_vals(k) + 3, ...
             sprintf('%.1f%%', sf_vals(k)), ...
             'HorizontalAlignment', 'center', ...
             'FontSize', 12, 'FontWeight', 'bold');
    end
end

hold off;
end

%% ========================================================================
%  FIGURE 6: Multi-panel EQE showing effect of each parameter
%% ========================================================================
if ismember(6, figures_to_plot)
figure('Name', 'EQE_Parameter_Effects', 'Position', [50 50 1000 800], 'Color', 'w');

%% ------------------------------------------------------------------------
% Panel A: Effect of eta_SF
%% ------------------------------------------------------------------------
subplot(2, 2, 1);
hold on;

eta_vals   = [2, 1.5, 1, 0.5, 0];

% Colorblind-friendly: blue (best) -> orange (worst)
eta_colors = [ ...
    0.12 0.47 0.71;   % eta = 2.0  - blue (best)
    0.50 0.70 0.86;   % eta = 1.5  - light blue
    0.60 0.60 0.60;   % eta = 1.0  - grey
    0.93 0.60 0.25;   % eta = 0.5  - light orange
    0.85 0.37 0.01];  % eta = 0.0  - orange (worst)

h_eta = gobjects(numel(eta_vals), 1);

for i = 1:numel(eta_vals)
    eta_val = eta_vals(i);

    mask = (T_success.SF_flag == 1) & ...
           (abs(T_success.eta_SF - eta_val) < tol) & ...
           (T_success.SRV_cm_s == optimal_SRV) & ...
           (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
           (T_success.junction_depth_nm == optimal_jd) & ...
           (abs(T_success.Phi_R_offset_eV - 0) < tol);

    if any(mask)
        run_idx = T_success.run_number(find(mask, 1));
        r = runs(run_idx);

        h_eta(i) = plot(r.EQE.wavelengths, r.EQE.EQE, '-', ...
                        'LineWidth', 2, ...
                        'Color', eta_colors(i, :));
    end
end

% 100% line (black dashed instead of blue to avoid confusion with data)
plot([400 800], [100 100], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

xlabel('Wavelength (nm)', 'FontSize', 20);
ylabel('EQE (%)', 'FontSize', 20);
xlim([400 800]);
ylim([50 120]);
set(gca, 'FontSize', 13, 'Box', 'on');
grid on;

leg_labels_eta = arrayfun(@(v) sprintf('\\eta_{SF} = %s', num2str(v)), ...
                          eta_vals, 'UniformOutput', false);
legend(h_eta, leg_labels_eta, 'Location', 'southeast', 'FontSize', 13);

hold off;

%% ------------------------------------------------------------------------
% Panel B: Effect of SRV
%% ------------------------------------------------------------------------
subplot(2, 2, 2);
hold on;

% Low SRV = blue (good), high SRV = orange (bad)
SRV_colors = CB_RAMP5;

h_srv = gobjects(numel(SRV_list), 1);

for s = 1:length(SRV_list)
    mask = (T_success.SF_flag == 1) & ...
           (abs(T_success.eta_SF - 2) < tol) & ...
           (T_success.SRV_cm_s == SRV_list(s)) & ...
           (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
           (T_success.junction_depth_nm == optimal_jd) & ...
           (abs(T_success.Phi_R_offset_eV - 0) < tol);
    if any(mask)
        run_idx = T_success.run_number(find(mask, 1));
        r = runs(run_idx);
        h_srv(s) = plot(r.EQE.wavelengths, r.EQE.EQE, '-', 'LineWidth', 2, ...
                        'Color', SRV_colors(s, :), ...
                        'DisplayName', sprintf('s_p = %.0d', SRV_list(s)));
    end
end

plot([400 800], [100 100], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Wavelength (nm)', 'FontSize', 20);
ylabel('EQE (%)', 'FontSize', 20);
xlim([400 800]);
ylim([50 120]);
legend(h_srv, 'Location', 'southeast', 'FontSize', 13);
set(gca, 'FontSize', 13, 'Box', 'on');
grid on;
hold off;

%% ------------------------------------------------------------------------
% Panel C: Effect of junction depth j_d
%% ------------------------------------------------------------------------
subplot(2, 2, 3);
hold on;

jd_vals_panel = [200, 250, 300, 350, 400];

% Shallow = blue (good), deep = orange (bad)
jd_colors = CB_RAMP5;

h_jd = gobjects(numel(jd_vals_panel), 1);

for jj = 1:numel(jd_vals_panel)
    jd_val = jd_vals_panel(jj);

    mask = (T_success.SF_flag == 1) & ...
           (abs(T_success.eta_SF - 2) < tol) & ...
           (T_success.SRV_cm_s == 10000) & ...
           (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
           (T_success.junction_depth_nm == jd_val) & ...
           (abs(T_success.Phi_R_offset_eV - 0) < tol);

    if any(mask)
        run_idx = T_success.run_number(find(mask, 1));
        r = runs(run_idx);
        h_jd(jj) = plot(r.EQE.wavelengths, r.EQE.EQE, '-', ...
                        'LineWidth', 2, ...
                        'Color', jd_colors(jj, :));
    end
end

plot([400 800], [100 100], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

xlabel('Wavelength (nm)', 'FontSize', 20);
ylabel('EQE (%)', 'FontSize', 20);
xlim([400 800]);
ylim([50 120]);
set(gca, 'FontSize', 13, 'Box', 'on');
grid on;

leg_labels_jd = arrayfun(@(d) sprintf('j_d = %d nm', d), ...
                         jd_vals_panel, 'UniformOutput', false);
legend(h_jd, leg_labels_jd, 'Location', 'southeast', 'FontSize', 13);

hold off;

%% ------------------------------------------------------------------------
% Panel D: Effect of doping concentration N_D
%% ------------------------------------------------------------------------
subplot(2, 2, 4);
hold on;

ND_vals_panel = [8e19, 1e19, 1e17, 1e15];
EF_vals_panel = [-4.05, -4.10, -4.22, -4.34];

% High doping = orange (overloaded), low doping = blue
ND_colors = [ ...
    0.85 0.37 0.01;   % 8e19  - orange (highest doping)
    0.93 0.60 0.25;   % 1e19  - light orange
    0.50 0.70 0.86;   % 1e17  - light blue
    0.12 0.47 0.71];  % 1e15  - blue (lowest doping)

h_nd = gobjects(numel(ND_vals_panel), 1);

for nn = 1:numel(ND_vals_panel)
    mask = (T_success.SF_flag == 1) & ...
           (abs(T_success.eta_SF - 2) < tol) & ...
           (T_success.SRV_cm_s == optimal_SRV) & ...
           (abs(T_success.EF_emitter_eV - EF_vals_panel(nn)) < tol) & ...
           (T_success.junction_depth_nm == optimal_jd) & ...
           (abs(T_success.Phi_R_offset_eV - 0) < tol);

    if any(mask)
        run_idx = T_success.run_number(find(mask, 1));
        r = runs(run_idx);
        h_nd(nn) = plot(r.EQE.wavelengths, r.EQE.EQE, '-', ...
                        'LineWidth', 2, ...
                        'Color', ND_colors(nn, :));
    end
end

plot([400 800], [100 100], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

xlabel('Wavelength (nm)', 'FontSize', 20);
ylabel('EQE (%)', 'FontSize', 20);
xlim([400 800]);
ylim([50 120]);
set(gca, 'FontSize', 13, 'Box', 'on');
grid on;

leg_labels_nd = arrayfun(@(d) sprintf('N_D = %.0e cm^{-3}', d), ...
                         ND_vals_panel, 'UniformOutput', false);
legend(h_nd, leg_labels_nd, 'Location', 'southeast', 'FontSize', 13);

hold off;

end  % end of multi-panel figure 6

%% ========================================================================
%  Extra Figure: Effect of \Phi_R offset at s_p = 10^4 cm/s
%% ========================================================================
if ismember(6, figures_to_plot)
figure('Name', 'EQE_vs_PhiRoffset_sp1e4', ...
       'Position', [200 200 430 400], 'Color', 'w');
hold on;

PhiOff_vals   = [0.0539, 0, -0.0539];
PhiOff_colors = [CB_BLUE; 0.50 0.50 0.50; CB_ORANGE];
PhiOff_labels = {'Q_f>0', 'Q_f=0', 'Q_f<0'};

SRV_plot = 1e4;

nPhi = numel(PhiOff_vals);
h_phi = gobjects(nPhi,1);

for p = 1:nPhi
    mask = (T_success.SF_flag == 1) & ...
           (abs(T_success.eta_SF - 2) < tol) & ...
           (T_success.SRV_cm_s == SRV_plot) & ...
           (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
           (T_success.junction_depth_nm == optimal_jd) & ...
           (abs(T_success.Phi_R_offset_eV - PhiOff_vals(p)) < tol);

    if any(mask)
        run_idx = T_success.run_number(find(mask, 1));
        r = runs(run_idx);
        h_phi(p) = plot(r.EQE.wavelengths, r.EQE.EQE, '-', ...
                        'LineWidth', 2, ...
                        'Color', PhiOff_colors(p, :), ...
                        'DisplayName', PhiOff_labels{p});
    end
end

plot([400 800], [100 100], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

xlabel('Wavelength (nm)', 'FontSize', 20);
ylabel('EQE (%)', 'FontSize', 20);
xlim([400 800]);
ylim([40 120]);
set(gca, 'FontSize', 15, 'Box', 'on');
grid on;

legend(h_phi, PhiOff_labels, 'Location', 'southeast', 'FontSize', 15);

hold off;
end

%% ========================================================================
%  FIGURE 8: Box Plot - PCE Distribution by eta_SF
%% ========================================================================
if ismember(8, figures_to_plot)
figure('Name', 'PCE_Distribution', 'Position', [350 150 500 400], 'Color', 'w');

PCE_by_cat = cell(4, 1);

mask = (T_success.SF_flag == 0);
PCE_by_cat{1} = T_success.PCE_pct(mask);

for e = 0:2
    mask = (T_success.SF_flag == 1) & (abs(T_success.eta_SF - e) < tol);
    PCE_by_cat{e+2} = T_success.PCE_pct(mask);
end

all_PCE    = [];
all_groups = [];
for k = 1:4
    all_PCE    = [all_PCE;    PCE_by_cat{k}];
    all_groups = [all_groups; k * ones(length(PCE_by_cat{k}), 1)];
end

boxplot(all_PCE, all_groups, 'Colors', 'k', 'Widths', 0.6);
hold on;

set(gca, 'XTick', 1:4, ...
         'XTickLabel', {'Baseline', 'SF (\eta=0)', 'SF (\eta=1)', 'SF (\eta=2)'}, ...
         'TickLabelInterpreter', 'tex');

% Scatter with colorblind-friendly categorical colors
for k = 1:4
    x_jitter = k + 0.15 * (rand(length(PCE_by_cat{k}), 1) - 0.5);
    scatter(x_jitter, PCE_by_cat{k}, 30, CB_CAT4(k,:), ...
            'filled', 'MarkerFaceAlpha', 0.5);
end

ylim([8, 21]);
ylabel('PCE (%)', 'FontSize', 20);
set(gca, 'FontSize', 18, 'LineWidth', 1.2, 'Box', 'on');
grid on;

med_vals = zeros(1,4);
for k = 1:4
    med_vals(k) = median(PCE_by_cat{k}, 'omitnan');
end

yl = ylim;
y_top_vals = yl(2) - 0.08*(yl(2)-yl(1));
y_top_text = yl(2) - 0.01*(yl(2)-yl(1));

for k = 1:4
    if ~isnan(med_vals(k))
        text(k, y_top_vals, ...
             sprintf('%.2f%%', med_vals(k)), ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'top', ...
             'FontSize', 16, ...
             'Color', 'k', ...
             'Clipping', 'off');
    end
end

text(2.5, y_top_text, 'Median PCE', ...
     'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'top', ...
     'FontSize', 14, ...
     'Clipping', 'off');

hold off;
end

%% ========================================================================
%  FIGURE 9: Correlation Matrix - All Parameters vs Metrics
%% ========================================================================
if ismember(9, figures_to_plot)
figure('Name', 'Correlation_Matrix', 'Position', [400 200 700 600], 'Color', 'w');

vars = {'SF', '\eta_{SF}', 'j_d', 'log(SRV)', 'E_F', '\Phi_R off', 'PCE', 'V_{oc}', 'J_{sc}', 'FF'};
data_matrix = [T_success.SF_flag, T_success.eta_SF, T_success.junction_depth_nm, ...
               log10(T_success.SRV_cm_s), T_success.EF_emitter_eV, T_success.Phi_R_offset_eV, ...
               T_success.PCE_pct, T_success.Voc_V, T_success.Jsc_mA_cm2, T_success.FF_pct];

R = corrcoef(data_matrix, 'Rows', 'complete');

imagesc(R);
colormap(blueorange(256));
caxis([-1 1]);
cb = colorbar;
ylabel(cb, 'Correlation', 'FontSize', 12);

set(gca, 'XTick', 1:length(vars), 'XTickLabel', vars, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:length(vars), 'YTickLabel', vars);
title('Parameter-Metric Correlation Matrix', 'FontSize', 16);
set(gca, 'FontSize', 11);

for i = 1:size(R, 1)
    for j = 1:size(R, 2)
        if abs(R(i,j)) > 0.3
            text(j, i, sprintf('%.2f', R(i,j)), 'HorizontalAlignment', 'center', ...
                 'FontSize', 9, 'Color', 'k', 'FontWeight', 'bold');
        end
    end
end

end

%% ========================================================================
%  FIGURE 10: Parallel Coordinates - Top vs Bottom Performers
%% ========================================================================
if ismember(10, figures_to_plot)
figure('Name', 'Parallel_Coords', 'Position', [450 250 800 400], 'Color', 'w');

n_select = max(10, round(0.1 * height(T_success)));
[~, sort_idx] = sort(T_success.PCE_pct, 'descend');
top_idx = sort_idx(1:n_select);
bottom_idx = sort_idx(end-n_select+1:end);

param_names = {'SF', '\eta_{SF}', 'j_d (nm)', 'SRV', 'E_F (eV)', '\Phi_R off'};
data_norm = zeros(height(T_success), 6);

data_norm(:, 1) = T_success.SF_flag;
data_norm(:, 2) = T_success.eta_SF / 2;
data_norm(:, 3) = (T_success.junction_depth_nm - 200) / 200;
data_norm(:, 4) = (log10(T_success.SRV_cm_s) - 1) / 4;
data_norm(:, 5) = (T_success.EF_emitter_eV - (-4.35)) / 0.3;
data_norm(:, 6) = (T_success.Phi_R_offset_eV - (-0.06)) / 0.12;

hold on;
for i = bottom_idx'
    plot(1:6, data_norm(i, :), '-', 'Color', [CB_ORANGE 0.3], 'LineWidth', 1);
end
for i = top_idx'
    plot(1:6, data_norm(i, :), '-', 'Color', [CB_BLUE 0.6], 'LineWidth', 1.5);
end

h1 = plot(NaN, NaN, '-', 'Color', CB_BLUE, 'LineWidth', 2);
h2 = plot(NaN, NaN, '-', 'Color', CB_ORANGE, 'LineWidth', 2);

set(gca, 'XTick', 1:6, 'XTickLabel', param_names);
xlim([0.5 6.5]);
ylim([-0.1 1.1]);
ylabel('Normalized Value', 'FontSize', 20);
legend([h1 h2], {'Top 10% PCE', 'Bottom 10% PCE'}, 'Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 18, 'LineWidth', 1.2, 'Box', 'on');
grid on;
hold off;

end

%% ========================================================================
%  FIGURE 11: Scatter Matrix - Key Parameters vs PCE
%% ========================================================================
if ismember(11, figures_to_plot)
figure('Name', 'Scatter_Matrix', 'Position', [100 50 1000 800], 'Color', 'w');

subplot(2, 2, 1);
scatter(T_success.eta_SF, T_success.PCE_pct, 40, T_success.SF_flag, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('\eta_{SF}', 'FontSize', 14);
ylabel('PCE (%)', 'FontSize', 14);
title('PCE vs \eta_{SF}', 'FontSize', 14);
colormap(gca, [CB_ORANGE; CB_BLUE]);
set(gca, 'FontSize', 12, 'Box', 'on');
grid on;

subplot(2, 2, 2);
scatter(log10(T_success.SRV_cm_s), T_success.PCE_pct, 40, T_success.eta_SF, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('log_{10}(SRV)', 'FontSize', 14);
ylabel('PCE (%)', 'FontSize', 14);
title('PCE vs SRV', 'FontSize', 14);
colormap(gca, blueorange(256));
cb = colorbar; ylabel(cb, '\eta_{SF}');
set(gca, 'FontSize', 12, 'Box', 'on');
grid on;

subplot(2, 2, 3);
scatter(T_success.junction_depth_nm, T_success.PCE_pct, 40, T_success.eta_SF, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('Junction Depth (nm)', 'FontSize', 14);
ylabel('PCE (%)', 'FontSize', 14);
title('PCE vs Junction Depth', 'FontSize', 14);
colormap(gca, blueorange(256));
cb = colorbar; ylabel(cb, '\eta_{SF}');
set(gca, 'FontSize', 12, 'Box', 'on');
grid on;

subplot(2, 2, 4);
scatter(T_success.EF_emitter_eV, T_success.PCE_pct, 40, log10(T_success.SRV_cm_s), 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('E_F (eV)', 'FontSize', 14);
ylabel('PCE (%)', 'FontSize', 14);
title('PCE vs Fermi Level', 'FontSize', 14);
colormap(gca, blueorange(256));
cb = colorbar; ylabel(cb, 'log_{10}(SRV)');
set(gca, 'FontSize', 12, 'Box', 'on');
grid on;

end

%% ========================================================================
%  FIGURE 12: Delta PCE Heatmap (SF benefit = PCE_SF - PCE_baseline)
%% ========================================================================
if ismember(12, figures_to_plot)
figure('Name', 'Delta_PCE_Heatmap', 'Position', [150 100 600 450], 'Color', 'w');

jd_fixed = 300;
EF_vals  = unique(T_success.EF_emitter_eV);
SRV_vals = unique(T_success.SRV_cm_s);

% --- Map E_F values to donor concentration N_D for axis labels ---
EF_map = [-4.05, -4.10, -4.22, -4.34];
ND_map = [8.06e19, 1e19, 1e17, 1e15];

ND_axis = NaN(size(EF_vals));
for ii = 1:numel(EF_vals)
    [~, idx_map] = min(abs(EF_map - EF_vals(ii)));
    ND_axis(ii) = ND_map(idx_map);
end

delta_PCE = NaN(length(EF_vals), length(SRV_vals));

for i = 1:length(EF_vals)
    for j = 1:length(SRV_vals)
        mask_base = (T_success.SF_flag == 0) & ...
                    (T_success.junction_depth_nm == jd_fixed) & ...
                    (T_success.SRV_cm_s == SRV_vals(j)) & ...
                    (abs(T_success.EF_emitter_eV - EF_vals(i)) < tol) & ...
                    (abs(T_success.Phi_R_offset_eV - 0) < tol);
        
        mask_sf = (T_success.SF_flag == 1) & ...
                  (abs(T_success.eta_SF - 2) < tol) & ...
                  (T_success.junction_depth_nm == jd_fixed) & ...
                  (T_success.SRV_cm_s == SRV_vals(j)) & ...
                  (abs(T_success.EF_emitter_eV - EF_vals(i)) < tol) & ...
                  (abs(T_success.Phi_R_offset_eV - 0) < tol);
        
        if any(mask_base) && any(mask_sf)
            PCE_base = mean(T_success.PCE_pct(mask_base), 'omitnan');
            PCE_sf   = mean(T_success.PCE_pct(mask_sf),   'omitnan');
            delta_PCE(i, j) = PCE_sf - PCE_base;
        end
    end
end

h = imagesc(1:length(SRV_vals), 1:length(EF_vals), delta_PCE);
colormap(blueorange(256));
set(h, 'AlphaData', ~isnan(delta_PCE));
set(gca, 'Color', 'w');

cb = colorbar;
ylabel(cb, '\Delta PCE (%) [SF - Baseline]', 'FontSize', 18);
set(gca, 'YDir', 'normal');

set(gca, 'XTick', 1:length(SRV_vals), ...
         'XTickLabel', arrayfun(@(x) sprintf('%.0e', x), SRV_vals, 'UniformOutput', false));
set(gca, 'YTick', 1:length(EF_vals), ...
         'YTickLabel', arrayfun(@(x) sprintf('%.1e', x), ND_axis, 'UniformOutput', false));

xlabel('Surface recombination velocity s_p (cm/s)', 'FontSize', 20);
ylabel('n-emitter doping concentration N_D (cm^{-3})', 'FontSize', 20);
set(gca, 'FontSize', 16);

for i = 1:length(EF_vals)
    for j = 1:length(SRV_vals)
        if ~isnan(delta_PCE(i,j))
            
            if delta_PCE(i,j) > 0
                text(j, i, sprintf('+%.2f', delta_PCE(i,j)), ...
                     'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', 'k');
            else
                text(j, i, sprintf('%.2f', delta_PCE(i,j)), ...
                     'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', 'w');
            end
        end
    end
end

end

%% ========================================================================
%  FIGURE 13: Radar Chart - Best vs Worst Configuration
%% ========================================================================
if ismember(13, figures_to_plot)
figure('Name', 'Radar_Chart', 'Position', [200 150 600 550], 'Color', 'w');

[~, best_idx] = max(T_success.PCE_pct);
[~, worst_idx] = min(T_success.PCE_pct);

metrics = {'PCE', 'V_{oc}', 'J_{sc}', 'FF', 'Blue EQE'};
n_metrics = length(metrics);

best_run = runs(T_success.run_number(best_idx));
worst_run = runs(T_success.run_number(worst_idx));

wl = best_run.EQE.wavelengths;
blue_mask = (wl >= 405) & (wl <= 500);
best_blue_eqe = mean(best_run.EQE.EQE(blue_mask), 'omitnan');
worst_blue_eqe = mean(worst_run.EQE.EQE(blue_mask), 'omitnan');

PCE_range = [min(T_success.PCE_pct), max(T_success.PCE_pct)];
Voc_range = [min(T_success.Voc_V), max(T_success.Voc_V)];
Jsc_range = [min(T_success.Jsc_mA_cm2), max(T_success.Jsc_mA_cm2)];
FF_range = [min(T_success.FF_pct), max(T_success.FF_pct)];
BlueEQE_range = [50, 130];

best_norm = [(T_success.PCE_pct(best_idx) - PCE_range(1)) / diff(PCE_range), ...
             (T_success.Voc_V(best_idx) - Voc_range(1)) / diff(Voc_range), ...
             (T_success.Jsc_mA_cm2(best_idx) - Jsc_range(1)) / diff(Jsc_range), ...
             (T_success.FF_pct(best_idx) - FF_range(1)) / diff(FF_range), ...
             (best_blue_eqe - BlueEQE_range(1)) / diff(BlueEQE_range)];

worst_norm = [(T_success.PCE_pct(worst_idx) - PCE_range(1)) / diff(PCE_range), ...
              (T_success.Voc_V(worst_idx) - Voc_range(1)) / diff(Voc_range), ...
              (T_success.Jsc_mA_cm2(worst_idx) - Jsc_range(1)) / diff(Jsc_range), ...
              (T_success.FF_pct(worst_idx) - FF_range(1)) / diff(FF_range), ...
              (worst_blue_eqe - BlueEQE_range(1)) / diff(BlueEQE_range)];

best_norm = [best_norm, best_norm(1)];
worst_norm = [worst_norm, worst_norm(1)];
angles = linspace(0, 2*pi, n_metrics + 1);

polarplot(angles, best_norm, '-', 'LineWidth', 2.5, 'Color', CB_BLUE, ...
          'DisplayName', 'Best config');
hold on;
polarplot(angles, worst_norm, '-', 'LineWidth', 2.5, 'Color', CB_ORANGE, ...
          'DisplayName', 'Worst config');

ax = gca;
ax.ThetaTick = rad2deg(angles(1:end-1));
ax.ThetaTickLabel = metrics;
ax.RLim = [0 1.1];
ax.FontSize = 12;

title('Performance Comparison: Best vs Worst Configuration', 'FontSize', 14);
legend('Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 12);
hold off;

end

%% ========================================================================
%  FIGURE 14: Stacked Area - EQE Decomposition
%% ========================================================================
if ismember(14, figures_to_plot)
figure('Name', 'EQE_Stacked', 'Position', [250 200 800 500], 'Color', 'w');

mask_base = (T_success.SF_flag == 0) & ...
            (T_success.SRV_cm_s == optimal_SRV) & ...
            (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
            (T_success.junction_depth_nm == optimal_jd) & ...
            (abs(T_success.Phi_R_offset_eV - 0) < tol);

mask_sf = (T_success.SF_flag == 1) & ...
          (abs(T_success.eta_SF - 2) < tol) & ...
          (T_success.SRV_cm_s == optimal_SRV) & ...
          (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
          (T_success.junction_depth_nm == optimal_jd) & ...
          (abs(T_success.Phi_R_offset_eV - 0) < tol);

if any(mask_base) && any(mask_sf)
    base_run = runs(T_success.run_number(find(mask_base, 1)));
    sf_run = runs(T_success.run_number(find(mask_sf, 1)));
    
    wl = base_run.EQE.wavelengths;
    eqe_base = base_run.EQE.EQE;
    eqe_sf = sf_run.EQE.EQE;
    
    area(wl, eqe_base, 'FaceColor', [0.50 0.70 0.86], 'FaceAlpha', 0.6, ...
         'EdgeColor', 'none', 'DisplayName', 'Baseline Si');
    hold on;
    
    eqe_diff = max(0, eqe_sf - eqe_base);
    area(wl, eqe_base + eqe_diff, 'FaceColor', CB_ORANGE, 'FaceAlpha', 0.6, ...
         'EdgeColor', 'none', 'DisplayName', 'SF Enhancement');
    
    plot(wl, eqe_sf, '-', 'LineWidth', 2, 'Color', CB_ORANGE, ...
         'DisplayName', 'Total (SF, \eta=2)');
    plot([400 800], [100 100], 'k--', 'LineWidth', 1.5, 'DisplayName', '100% EQE');
    
    xlabel('Wavelength (nm)', 'FontSize', 16);
    ylabel('EQE (%)', 'FontSize', 16);
    title('EQE Decomposition: Baseline + SF Enhancement', 'FontSize', 18);
    xlim([400 700]);
    ylim([0 150]);
    legend('Location', 'northeast', 'FontSize', 12);
    set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'on');
    grid on;
    hold off;
end

end

%% ========================================================================
%  FIGURE 15: Summary Table as Figure
%% ========================================================================
if ismember(15, figures_to_plot)
figure('Name', 'Summary_Table', 'Position', [300 250 700 400], 'Color', 'w');
axis off;

[max_PCE, best_idx] = max(T_success.PCE_pct);
best = T_success(best_idx, :);

summary_text = {
    'BEST CONFIGURATION:'
    sprintf('  PCE = %.3f %%', best.PCE_pct)
    sprintf('  Voc = %.4f V', best.Voc_V)
    sprintf('  Jsc = %.4f mA/cm^2', best.Jsc_mA_cm2)
    sprintf('  FF = %.2f %%', best.FF_pct)
    ''
    'PARAMETERS:'
    sprintf('  SF flag = %d', best.SF_flag)
    sprintf('  eta_SF = %.1f', best.eta_SF)
    sprintf('  Junction depth = %d nm', best.junction_depth_nm)
    sprintf('  SRV = %.0e cm/s', best.SRV_cm_s)
    sprintf('  E_F = %.3f eV', best.EF_emitter_eV)
    sprintf('  Phi_R offset = %.4f eV', best.Phi_R_offset_eV)
    ''
    'SWEEP STATISTICS:'
    sprintf('  Total runs: %d', summary.n_expected)
    sprintf('  Successful: %d (%.1f%%)', summary.n_success, 100*summary.n_success/summary.n_expected)
    sprintf('  Mean PCE (all): %.3f %%', mean(T_success.PCE_pct, 'omitnan'))
    sprintf('  Mean PCE (SF=1, eta=2): %.3f %%', mean(T_success.PCE_pct(T_success.SF_flag==1 & T_success.eta_SF==2), 'omitnan'))
};

text(0.1, 0.95, summary_text, 'VerticalAlignment', 'top', 'FontSize', 11, ...
     'FontName', 'FixedWidth');

title('Parameter Sweep Summary', 'FontSize', 18);

end

%% ========================================================================
%  FIGURE 16: Best Device - EQE
%% ========================================================================
if ismember(16, figures_to_plot)
figure('Name', 'Best_Device_EQE', 'Position', [350 300 800 550], 'Color', 'w');

[~, best_idx] = max(T_success.PCE_pct);
best = T_success(best_idx, :);
best_run = runs(best.run_number);

wl = best_run.EQE.wavelengths;
eqe = best_run.EQE.EQE;
iqe = best_run.EQE.IQE;

plot(wl, eqe, '-', 'LineWidth', 2.5, 'Color', CB_BLUE, 'DisplayName', 'EQE');
hold on;
plot(wl, iqe, '--', 'LineWidth', 2, 'Color', CB_ORANGE, 'DisplayName', 'IQE');

plot([400 800], [100 100], 'k--', 'LineWidth', 1.5, 'DisplayName', '100%');

% Shade SF-active region (400-550nm)
fill([400 550 550 400], [0 0 160 160], [0.85 0.92 1.0], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off');

% Shade >100% region
idx_above = eqe > 100;
if any(idx_above)
    wl_above = wl(idx_above);
    eqe_above = eqe(idx_above);
    fill([wl_above, fliplr(wl_above)], [100*ones(size(wl_above)), fliplr(eqe_above)], ...
         [0.85 0.92 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
end

xlabel('Wavelength (nm)', 'FontSize', 18);
ylabel('Quantum Efficiency (%)', 'FontSize', 18);
title(sprintf('Best Device EQE (PCE = %.3f%%)', best.PCE_pct), 'FontSize', 20);
xlim([400 800]);
ylim([0 160]);
legend('Location', 'northeast', 'FontSize', 14);
set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'on');
grid on;

[max_eqe, max_idx] = max(eqe);
text(wl(max_idx), max_eqe + 5, sprintf('Peak: %.1f%%', max_eqe), ...
     'FontSize', 12, 'HorizontalAlignment', 'center', 'Color', CB_BLUE);

param_str = sprintf(['SF = %d, \\eta_{SF} = %.0f\n' ...
                     'j_d = %d nm\n' ...
                     'SRV = %.0e cm/s\n' ...
                     'E_F = %.3f eV'], ...
                    best.SF_flag, best.eta_SF, best.junction_depth_nm, ...
                    best.SRV_cm_s, best.EF_emitter_eV);
annotation('textbox', [0.15, 0.65, 0.25, 0.2], 'String', param_str, ...
           'FontSize', 11, 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
           'FitBoxToText', 'on', 'Interpreter', 'tex');

hold off;

end

%% ========================================================================
%  FIGURE 17: Best Device - J-V Curve
%% ========================================================================
if ismember(17, figures_to_plot)
figure('Name', 'Best_Device_JV', 'Position', [400 350 800 550], 'Color', 'w');

V = best_run.JV.Vapp;
J = best_run.JV.Jtot;

[~, idx_max] = max(V);
V_fwd = V(1:idx_max);
J_fwd = J(1:idx_max);

plot(V_fwd, J_fwd, '-', 'LineWidth', 2.5, 'Color', CB_BLUE);
hold on;

plot(xlim, [0 0], 'k--', 'LineWidth', 1);

% Jsc
[~, idx_v0] = min(abs(V_fwd));
Jsc_plot = J_fwd(idx_v0);
plot(0, Jsc_plot, 'o', 'MarkerSize', 10, 'MarkerFaceColor', CB_ORANGE, ...
     'MarkerEdgeColor', CB_ORANGE);
text(0.02, Jsc_plot, sprintf('J_{sc} = %.4f mA/cm^2', abs(Jsc_plot)), ...
     'FontSize', 11, 'VerticalAlignment', 'bottom');

% Voc
Voc_plot = best.Voc_V;
plot(Voc_plot, 0, 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0.55 0.34 0.64], ...
     'MarkerEdgeColor', [0.55 0.34 0.64]);
text(Voc_plot, 0.02, sprintf('V_{oc} = %.4f V', Voc_plot), ...
     'FontSize', 11, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% MPP
power = V_fwd .* J_fwd;
[MPP, idx_mpp] = min(power);
V_mpp = V_fwd(idx_mpp(1));
J_mpp = J_fwd(idx_mpp(1));
plot(V_mpp, J_mpp, 's', 'MarkerSize', 12, 'MarkerFaceColor', [0.55 0.34 0.64], ...
     'MarkerEdgeColor', [0.55 0.34 0.64]);
text(V_mpp + 0.02, J_mpp, sprintf('MPP\nP = %.4f mW/cm^2', abs(MPP(1))), ...
     'FontSize', 10, 'VerticalAlignment', 'middle');

% MPP rectangle
fill([0 V_mpp V_mpp 0], [0 0 J_mpp J_mpp], [0.85 0.92 1.0], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.3);

xlabel('Voltage (V)', 'FontSize', 18);
ylabel('Current Density (mA/cm^2)', 'FontSize', 18);
title(sprintf('Best Device J-V Curve (PCE = %.3f%%, FF = %.2f%%)', ...
              best.PCE_pct, best.FF_pct), 'FontSize', 20);
set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'on');
grid on;

metrics_str = sprintf(['PCE = %.3f %%\n' ...
                       'V_{oc} = %.4f V\n' ...
                       'J_{sc} = %.4f mA/cm^2\n' ...
                       'FF = %.2f %%'], ...
                      best.PCE_pct, best.Voc_V, best.Jsc_mA_cm2, best.FF_pct);
annotation('textbox', [0.15, 0.15, 0.25, 0.2], 'String', metrics_str, ...
           'FontSize', 12, 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
           'FitBoxToText', 'on', 'Interpreter', 'tex');

hold off;

end

%% ========================================================================
%  FIGURE 18: Best Device vs Baseline Comparison (side by side)
%% ========================================================================
if ismember(18, figures_to_plot)
figure('Name', 'Best_vs_Baseline', ...
       'Position', [100 100 1400 500], ...
       'Color', 'w');

[~, best_idx] = max(T_success.PCE_pct);
best      = T_success(best_idx,:);
best_run  = runs(best.run_number);

mask_baseline = (T_success.SF_flag == 0) & ...
                (T_success.junction_depth_nm == best.junction_depth_nm) & ...
                (T_success.SRV_cm_s           == best.SRV_cm_s) & ...
                (abs(T_success.EF_emitter_eV      - best.EF_emitter_eV) < 0.01) & ...
                (abs(T_success.Phi_R_offset_eV    - best.Phi_R_offset_eV) < 0.01);

if any(mask_baseline)
    baseline     = T_success(find(mask_baseline, 1), :);
    baseline_run = runs(baseline.run_number);

    %% -------- Panel A: EQE comparison -----------------------------------
    subplot(1, 2, 1);
    hold on;
    
    % SF device EQE
    h_sf_eqe = plot(best_run.EQE.wavelengths, best_run.EQE.EQE, ...
                    '-', 'Color', CB_ORANGE, 'LineWidth', 2.5);
    
    % Baseline EQE
    h_bl_eqe = plot(baseline_run.EQE.wavelengths, baseline_run.EQE.EQE, ...
                    '-', 'Color', CB_BLUE, 'LineWidth', 2.5);
    
    plot([400 800], [100 100], 'k--', 'LineWidth', 1.5, ...
         'HandleVisibility', 'off');
    
    xlabel('Wavelength (nm)', 'FontSize', 24);
    ylabel('EQE (%)', 'FontSize', 24);
    xlim([400 800]);
    ylim([80 110]);
    set(gca, 'FontSize', 18, 'LineWidth', 1.2, 'Box', 'on');
    grid on;
    
    [eqe_max, idx_max_EQE] = max(best_run.EQE.EQE);
    lambda_max = best_run.EQE.wavelengths(idx_max_EQE);
    
    plot(lambda_max, eqe_max, 'o', ...
         'MarkerSize', 8, ...
         'MarkerFaceColor', CB_ORANGE, ...
         'MarkerEdgeColor', [0.65 0.25 0.00], ...
         'HandleVisibility', 'off');
    
    yl = ylim;
    y_offset = 0.02 * (yl(2) - yl(1));
    
    text(lambda_max + 5, eqe_max + y_offset, ...
         sprintf('%.1f%% @ %dnm', eqe_max, round(lambda_max)), ...
         'FontSize', 18, ...
         'Color', [0.65 0.25 0.00], ...
         'HorizontalAlignment', 'left', ...
         'VerticalAlignment', 'bottom');
    
    legend([h_bl_eqe, h_sf_eqe], ...
           {'Baseline', sprintf('SF (\\eta=%.0f)', best.eta_SF)}, ...
           'Location', 'northeast', 'FontSize', 16);
    
    hold off;

    %% -------- Panel B: J-V comparison -----------------------------------
    ax2 = subplot(1, 2, 2); hold on;

    V_best = best_run.JV.Vapp;
    J_best = best_run.JV.Jtot;
    [~, idx_max_best] = max(V_best);
    h_sf = plot(V_best(1:idx_max_best), J_best(1:idx_max_best), ...
                '-', 'Color', CB_ORANGE, 'LineWidth', 2.5);

    V_base = baseline_run.JV.Vapp;
    J_base = baseline_run.JV.Jtot;
    [~, idx_max_base] = max(V_base);
    h_bl = plot(V_base(1:idx_max_base), J_base(1:idx_max_base), ...
                '-', 'Color', CB_BLUE, 'LineWidth', 2.5);

    plot(xlim, [0 0], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    yl = ylim;
    plot([0 0], yl, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

    xlabel('Voltage (V)', 'FontSize', 24);
    ylabel('Current Density (mA/cm^2)', 'FontSize', 24);
    xlim([-0.3, 0.8]);
    ylim([-0.05, 0.1]);
    set(gca, 'FontSize', 18, 'LineWidth', 1.2, 'Box', 'on');
    grid on;

    base_line1 = sprintf('Baseline: PCE = %.3f%%            V_{oc} = %.4f V', ...
                         baseline.PCE_pct, baseline.Voc_V);
    base_line2 = sprintf('                J_{sc} = %.4f mA/cm^2   FF = %.2f%%', ...
                         baseline.Jsc_mA_cm2, baseline.FF_pct);

    sf_line1   = sprintf('SF (\\eta=%.0f): PCE = %.3f%%            V_{oc} = %.4f V', ...
                         best.eta_SF, best.PCE_pct, best.Voc_V);
    sf_line2   = sprintf('                J_{sc} = %.4f mA/cm^2   FF = %.2f%%', ...
                         best.Jsc_mA_cm2, best.FF_pct);

    h_dummy1 = plot(NaN, NaN, 'w', 'LineStyle', 'none', 'Marker', 'none');
    h_dummy2 = plot(NaN, NaN, 'w', 'LineStyle', 'none', 'Marker', 'none');

    legend([h_bl, h_dummy1, h_sf, h_dummy2], ...
           {base_line1, base_line2, sf_line1, sf_line2}, ...
           'Location', 'northwest', ...
           'FontSize', 16, ...
           'Interpreter', 'tex');

    hold off;
else
    text(0.5, 0.5, 'No matching baseline found', 'FontSize', 16, ...
         'HorizontalAlignment', 'center');
end

end

%% Summary
fprintf('\n========================================\n');
fprintf('Generated %d figures!\n', length(figures_to_plot));
fprintf('========================================\n');

end

%% ========================================================================
%  UTILITY FUNCTIONS
%% ========================================================================

function c = redblue(m)
%REDBLUE Red-white-blue colormap (kept for compatibility if needed)
    if nargin < 1, m = 256; end
    
    r = [linspace(0, 1, m/2), ones(1, m/2)];
    g = [linspace(0, 1, m/2), linspace(1, 0, m/2)];
    b = [ones(1, m/2), linspace(1, 0, m/2)];
    
    c = [r(:), g(:), b(:)];
end

function c = blueorange(m)
%BLUEORANGE  Colorblind-friendly diverging colormap (blue-white-orange)
%   Replaces redgreen for accessibility.
%   Low values  -> blue   [0.12  0.47  0.71]
%   Mid values  -> white  [1.00  1.00  1.00]
%   High values -> orange [0.85  0.37  0.01]

    if nargin < 1
        m = 256;
    end

    n1 = floor(m/2);
    n2 = m - n1;

    % Blue to white
    r1 = linspace(0.12, 1.0, n1);
    g1 = linspace(0.47, 1.0, n1);
    b1 = linspace(0.71, 1.0, n1);

    % White to orange
    r2 = linspace(1.0, 0.85, n2);
    g2 = linspace(1.0, 0.37, n2);
    b2 = linspace(1.0, 0.01, n2);

    r = [r1, r2];
    g = [g1, g2];
    b = [b1, b2];

    c = [r(:), g(:), b(:)];
end

function c = redgreen(m)
%REDGREEN  (DEPRECATED - kept for backward compatibility, use blueorange instead)
    c = blueorange(m);
end

function mEQE = local_mean_eqe_at_lambda(runs, run_idx_list, lambda_target)
%LOCAL_MEAN_EQE_AT_LAMBDA  Mean EQE at a target wavelength over a set of runs.
    eqe_vals = NaN(numel(run_idx_list), 1);
    for ii = 1:numel(run_idx_list)
        r  = runs(run_idx_list(ii));
        wl = r.EQE.wavelengths;
        [~, idx] = min(abs(wl - lambda_target));
        eqe_vals(ii) = r.EQE.EQE(idx);
    end
    mEQE = mean(eqe_vals, 'omitnan');
end