function generate_poster_figures_reversebias(summaryFile)
%GENERATE_POSTER_FIGURES_RB Create publication-quality figures for research poster
%   Modified for the reverse-bias sweep (run_sweep_chunk with SR data)
%
%   generate_poster_figures_rb('sweep_summary_rb.mat')
%
%   Generates selected figures based on the FIGURES_TO_PLOT configuration below.
%
%   Available figures:
%     1  - EQE comparison showing >100% region (SF signature)
%     2  - Heatmap: PCE vs (SRV, N_D)
%     3  - Line plot: PCE vs SRV (baseline vs SF)
%     4  - Line plot: PCE vs junction depth  [DISABLED - only 1 jd in this sweep]
%     5  - Line plot: Integrated blue EQE vs eta_SF
%     6  - 3-panel EQE: effect of eta, SRV, Phi_R offset
%     7  - [reserved]
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
%   NEW (reverse-bias spectral response):
%     19 - SR_EQE under bias: effect of bias voltage (single device)
%     20 - SR_EQE under bias: SF vs baseline at each bias point
%     21 - SR_EQE ratio (reverse-bias / short-circuit) vs wavelength
%     22 - Heatmap: integrated blue SR_EQE vs (SRV, bias)
%     23 - SR_EQE under bias: effect of SRV at each bias point
%     24 - SR_EQE under bias: effect of eta_SF at each bias point
%     25 - Delta EQE (reverse-bias minus short-circuit) decomposition

close all;

%% ========================================================================
%  CONFIGURATION - Select which figures to generate
%% ========================================================================

% Set to 'all' to plot everything, or specify figure numbers as array
% FIGURES_TO_PLOT = [19,20,21,22,23,24,25];
FIGURES_TO_PLOT = 'all';

%% ========================================================================
%  END CONFIGURATION
%% ========================================================================

% Process figure selection
if ischar(FIGURES_TO_PLOT) && strcmp(FIGURES_TO_PLOT, 'all')
    figures_to_plot = [1:3, 5:6, 8:25];  % skip 4 (only 1 jd) and 7
else
    figures_to_plot = FIGURES_TO_PLOT;
end

fprintf('Will generate figures: %s\n', mat2str(figures_to_plot));

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
% Shared sweep settings — UPDATED for reverse-bias sweep
% ------------------------------------------------------------------------
% New sweep parameter space:
%   junction_depth_vals = [200]
%   sp_r_vals = [1000, 10000, 100000]
%   eta_sf_with_SF = [0, 1, 2]
%   EF_vals = [-4.05, -4.104, -4.223, -4.342]
%   Phi_R_offset_vals = [0.0540, 0, -0.0540]
%   SR_bias_vals = [-1.0, -0.5, 0.0]

optimal_SRV   = 1000;       % Best SRV in this sweep
optimal_EF    = -4.104;     % ~1e19 doping
optimal_PhiOff = 0;
optimal_jd    = 200;         % Only junction depth in this sweep
tol           = 0.01;

% SRV levels in this sweep
SRV_list = [10, 100, 1000, 10000, 100000];

% Bias voltages used in spectral response
SR_bias_vals = [-1.0, -0.5, 0.0];

% EF -> N_D mapping (unchanged)
EF_map = [-4.05, -4.10, -4.22, -4.34];
ND_map = [8.06e19, 1e19, 1e17, 1e15];

%% ========================================================================
%  FIGURE 1: EQE Comparison - SF Signature (EQE > 100%)
%% ========================================================================
if ismember(1, figures_to_plot)
figure('Name', 'EQE_SF_Comparison', 'Position', [50 50 900 600], 'Color', 'w');

for eta_val = [0, 2]
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
            plot(r.EQE.wavelengths, r.EQE.EQE, 'r-', 'LineWidth', 2.5, ...
                 'DisplayName', 'SF layer, \eta_{SF} = 0');
        else
            plot(r.EQE.wavelengths, r.EQE.EQE, 'g-', 'LineWidth', 2.5, ...
                 'DisplayName', 'SF layer, \eta_{SF} = 2');
        end
        hold on;
    end
end

% Baseline (SF=0)
mask = (T_success.SF_flag == 0) & ...
       (T_success.SRV_cm_s == optimal_SRV) & ...
       (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
       (abs(T_success.Phi_R_offset_eV - optimal_PhiOff) < tol) & ...
       (T_success.junction_depth_nm == optimal_jd);
idx = find(mask);
if ~isempty(idx)
    run_idx = T_success.run_number(idx(1));
    r = runs(run_idx);
    plot(r.EQE.wavelengths, r.EQE.EQE, 'b-', 'LineWidth', 2.5, ...
         'DisplayName', 'No SF layer (baseline)');
end

xlims = [400 650];
plot(xlims, [100 100], 'k--', 'LineWidth', 1.5, 'DisplayName', 'EQE = 100%');
fill([xlims(1) xlims(2) xlims(2) xlims(1)], [100 100 150 150], ...
     [0.9 1 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'HandleVisibility', 'off');

xlabel('Wavelength (nm)', 'FontSize', 18);
ylabel('EQE (%)', 'FontSize', 18);
title('Singlet Fission Enhancement of EQE', 'FontSize', 20);
xlim(xlims); ylim([0 150]);
set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'on');
grid on;
legend('Location', 'northeast', 'FontSize', 12);
text(420, 130, 'SF benefit region', 'FontSize', 12, 'Color', [0 0.5 0]);
hold off;
end

%% ========================================================================
%  FIGURE 2: Heatmap - PCE vs (SRV, N_D) for ideal SF
%% ========================================================================
if ismember(2, figures_to_plot)
figure('Name', 'Heatmap_PCE_SRV_ND', 'Position', [100 100 600 450], 'Color', 'w');

mask = (T_success.SF_flag == 1) & ...
       (abs(T_success.eta_SF - 2) < tol) & ...
       (T_success.junction_depth_nm == optimal_jd) & ...
       (abs(T_success.Phi_R_offset_eV - 0) < tol);

T_sub = T_success(mask, :);
SRV_vals = unique(T_sub.SRV_cm_s);
EF_vals  = unique(T_sub.EF_emitter_eV);

ND_axis = NaN(size(EF_vals));
for ii = 1:numel(EF_vals)
    [~, idx_map] = min(abs(EF_map - EF_vals(ii)));
    ND_axis(ii) = ND_map(idx_map);
end

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

h = imagesc(1:length(SRV_vals), 1:length(EF_vals), PCE_matrix);
colormap(redgreen(256));
set(h, 'AlphaData', ~isnan(PCE_matrix));
set(gca, 'Color', 'w');

cb = colorbar;
ylabel(cb, 'PCE (%)', 'FontSize', 18);
set(gca, 'YDir', 'normal');

set(gca, 'XTick', 1:length(SRV_vals), ...
         'XTickLabel', arrayfun(@(x) sprintf('%.0e', x), SRV_vals, 'UniformOutput', false));
set(gca, 'YTick', 1:length(EF_vals), ...
         'YTickLabel', arrayfun(@(x) sprintf('%.1e', x), ND_axis, 'UniformOutput', false));

xlabel('Surface recombination velocity s_p (cm/s)', 'FontSize', 20);
ylabel('n-emitter doping concentration N_D (cm^{-3})', 'FontSize', 20);
set(gca, 'FontSize', 16, 'LineWidth', 1.2);

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

nSRV = numel(SRV_list);
PCE_baseline = zeros(nSRV, 1);
PCE_SF       = zeros(nSRV, 1);

for k = 1:nSRV
    mask = (T_success.SF_flag == 0) & ...
           (T_success.SRV_cm_s == SRV_list(k)) & ...
           (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
           (T_success.junction_depth_nm == optimal_jd) & ...
           (abs(T_success.Phi_R_offset_eV - 0) < tol);
    if any(mask)
        PCE_baseline(k) = mean(T_success.PCE_pct(mask), 'omitnan');
    else
        PCE_baseline(k) = NaN;
    end
    
    mask = (T_success.SF_flag == 1) & ...
           (abs(T_success.eta_SF - 2) < tol) & ...
           (T_success.SRV_cm_s == SRV_list(k)) & ...
           (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
           (T_success.junction_depth_nm == optimal_jd) & ...
           (abs(T_success.Phi_R_offset_eV - 0) < tol);
    if any(mask)
        PCE_SF(k) = mean(T_success.PCE_pct(mask), 'omitnan');
    else
        PCE_SF(k) = NaN;
    end
end

hold on;
h1 = semilogx(SRV_list, PCE_baseline, '-o', ...
              'LineWidth', 2.5, 'MarkerSize', 10, ...
              'MarkerFaceColor', [0.27 0.33 0.80], ...
              'Color', [0.27 0.33 0.80], 'DisplayName', 'Baseline');

h2 = semilogx(SRV_list, PCE_SF, '-s', ...
              'LineWidth', 2.5, 'MarkerSize', 10, ...
              'MarkerFaceColor', [0.13 0.55 0.13], ...
              'Color', [0.13 0.55 0.13], 'DisplayName', 'SF (\eta=2)');

ylabel('PCE (%)', 'FontSize', 20);
xlabel('Surface recombination velocity s_p (cm/s)', 'FontSize', 20);
legend([h1 h2], {'Baseline', 'SF (\eta=2)'}, 'Location', 'northeast', 'FontSize', 16);

set(gca, 'XScale', 'log', 'XTick', SRV_list);
set(gca, 'TickLabelInterpreter', 'tex');
xticklabels(arrayfun(@(x) sprintf('10^{%d}', round(log10(x))), SRV_list, 'UniformOutput', false));
ax = gca; ax.XMinorTick = 'off'; ax.XMinorGrid = 'off';
grid on;
xlim([min(SRV_list)*0.8, max(SRV_list)*1.2]);
set(gca, 'Box', 'on', 'LineWidth', 1.5, 'FontSize', 18);

ymin = nanmin([PCE_baseline(:); PCE_SF(:)]);
ymax = nanmax([PCE_baseline(:); PCE_SF(:)]);
dy = 0.02*(ymax - ymin);

for k = 1:nSRV
    if ~isnan(PCE_baseline(k))
        text(SRV_list(k), PCE_baseline(k) + dy, sprintf('%.2f%%', PCE_baseline(k)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontSize', 16, 'Color', [0.1 0.1 0.4]);
    end
    if ~isnan(PCE_SF(k))
        text(SRV_list(k), PCE_SF(k) + 2*dy, sprintf('%.2f%%', PCE_SF(k)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontSize', 16, 'Color', [0.0 0.3 0.0]);
    end
end
hold off;
end

%% ========================================================================
%  FIGURE 4: SKIPPED (only 1 junction depth in this sweep)
%% ========================================================================
if ismember(4, figures_to_plot)
    fprintf('Figure 4 skipped: only one junction depth (200 nm) in this sweep.\n');
end

%% ========================================================================
%  FIGURE 5: Integrated EQE in Blue Region (405-500nm) - SF Signature
%% ========================================================================
if ismember(5, figures_to_plot)
figure('Name', 'Blue_EQE_Line', 'Position', [250 250 700 500], 'Color', 'w');

eta_vals = [0, 1, 2];  % eta values in this sweep
n_eta = numel(eta_vals);
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
    h_base = yline(blue_EQE(1), '--', 'Baseline', 'LineWidth', 2, 'Color', [0.30 0.30 0.80]);
else
    h_base = [];
end

sf_vals = blue_EQE(2:end);
valid_sf = ~isnan(sf_vals);
h_sf = plot(eta_vals(valid_sf), sf_vals(valid_sf), '-o', ...
            'LineWidth', 2.5, 'MarkerSize', 8, ...
            'Color', [0.00 0.45 0.00], 'MarkerFaceColor', [0.30 0.75 0.30], ...
            'DisplayName', 'SF device');

h_100 = yline(100, 'k--', 'EQE = 100%', 'LineWidth', 1.5);

xlabel('\eta_{SF}', 'FontSize', 16);
ylabel('Mean EQE in Blue Region (405-500nm) (%)', 'FontSize', 16);
title('Blue-Region EQE Enhancement vs \eta_{SF}', 'FontSize', 18);
xlim([min(eta_vals) - 0.1, max(eta_vals) + 0.1]);

valid_vals = blue_EQE(~isnan(blue_EQE));
if ~isempty(valid_vals), ylim([0, max(valid_vals)*1.15]); end

set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'on');
grid on;

if ~isempty(h_base)
    legend([h_sf, h_base, h_100], {'SF device', 'Baseline', 'EQE = 100%'}, ...
           'Location', 'best', 'FontSize', 12);
else
    legend([h_sf, h_100], {'SF device', 'EQE = 100%'}, 'Location', 'best', 'FontSize', 12);
end

for k = 1:n_eta
    if ~isnan(sf_vals(k))
        text(eta_vals(k), sf_vals(k) + 3, sprintf('%.1f%%', sf_vals(k)), ...
             'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
    end
end
hold off;
end

%% ========================================================================
%  FIGURE 6: Multi-panel EQE showing effect of each parameter
%% ========================================================================
if ismember(6, figures_to_plot)
figure('Name', 'EQE_Parameter_Effects', 'Position', [50 50 1300 400], 'Color', 'w');

%% Panel A: Effect of eta_SF
subplot(1, 3, 1);
hold on;

eta_vals = [2, 1, 0];
eta_colors = [0.00 0.60 0.00; 0.90 0.80 0.20; 0.80 0.00 0.00];
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
                        'LineWidth', 2, 'Color', eta_colors(i, :));
    end
end

plot([400 800], [100 100], 'b--', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Wavelength (nm)', 'FontSize', 24);
ylabel('EQE (%)', 'FontSize', 24);
xlim([400 800]); ylim([50 120]);
set(gca, 'FontSize', 14, 'Box', 'on');
grid on;

leg_labels_eta = arrayfun(@(v) sprintf('\\eta_{SF} = %s', num2str(v)), eta_vals, 'UniformOutput', false);
valid = isgraphics(h_eta);
legend(h_eta(valid), leg_labels_eta(valid), 'Location', 'southeast', 'FontSize', 14);
hold off;

%% Panel B: Effect of SRV
subplot(1, 3, 2);
hold on;

SRV_colors = [0.00 0.45 0.00; 0.90 0.80 0.20; 0.80 0.00 0.00];
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
                        'DisplayName', sprintf('s_p = %.0d cm/s', SRV_list(s)));
    end
end

plot([400 800], [100 100], 'b--', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Wavelength (nm)', 'FontSize', 24);
ylabel('EQE (%)', 'FontSize', 24);
xlim([400 800]); ylim([50 120]);
valid = isgraphics(h_srv);
legend(h_srv(valid), 'Location', 'southeast', 'FontSize', 14);
set(gca, 'FontSize', 14, 'Box', 'on');
grid on;
hold off;

%% Panel C: Effect of Phi_R offset
subplot(1, 3, 3);
hold on;

PhiOff_vals   = [0.0540, 0, -0.0540];
PhiOff_colors = [0 0.5 0; 0.5 0.5 0.5; 0.8 0 0];
PhiOff_labels = {'Q_f > 0', 'Q_f = 0', 'Q_f < 0'};

SRV_plot = 1e4;
nPhi = numel(PhiOff_vals);
h_phi = gobjects(nPhi, 1);

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
                        'LineWidth', 2, 'Color', PhiOff_colors(p, :), ...
                        'DisplayName', PhiOff_labels{p});
    end
end

plot([400 800], [100 100], 'b--', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Wavelength (nm)', 'FontSize', 24);
ylabel('EQE (%)', 'FontSize', 24);
xlim([400 800]); ylim([40 120]);
set(gca, 'FontSize', 14, 'Box', 'on');
grid on;
valid = isgraphics(h_phi);
legend(h_phi(valid), PhiOff_labels(valid), 'Location', 'southeast', 'FontSize', 14);
hold off;
end

%% ========================================================================
%  FIGURE 8: Box Plot - PCE Distribution by eta_SF
%% ========================================================================
if ismember(8, figures_to_plot)
figure('Name', 'PCE_Distribution', 'Position', [350 150 500 400], 'Color', 'w');

% Categories: Baseline, SF eta=0, SF eta=1, SF eta=2
PCE_by_cat = cell(4, 1);

mask = (T_success.SF_flag == 0);
PCE_by_cat{1} = T_success.PCE_pct(mask);

for idx_e = 1:3
    eta_val = [0, 1, 2];
    mask = (T_success.SF_flag == 1) & (abs(T_success.eta_SF - eta_val(idx_e)) < tol);
    PCE_by_cat{idx_e+1} = T_success.PCE_pct(mask);
end

all_PCE = []; all_groups = [];
for k = 1:4
    all_PCE    = [all_PCE;    PCE_by_cat{k}];
    all_groups = [all_groups; k * ones(length(PCE_by_cat{k}), 1)];
end

boxplot(all_PCE, all_groups, 'Colors', 'k', 'Widths', 0.6);
hold on;

set(gca, 'XTick', 1:4, ...
         'XTickLabel', {'Baseline', 'SF (\eta=0)', 'SF (\eta=1)', 'SF (\eta=2)'}, ...
         'TickLabelInterpreter', 'tex');

colors_cat = [0 0 0.8; 0.8 0 0; 1 0.6 0; 0 0.6 0];
for k = 1:4
    x_jitter = k + 0.15 * (rand(length(PCE_by_cat{k}), 1) - 0.5);
    scatter(x_jitter, PCE_by_cat{k}, 30, colors_cat(k,:), 'filled', 'MarkerFaceAlpha', 0.5);
end

ylabel('PCE (%)', 'FontSize', 20);
set(gca, 'FontSize', 18, 'LineWidth', 1.2, 'Box', 'on');
grid on;

med_vals = zeros(1,4);
for k = 1:4, med_vals(k) = median(PCE_by_cat{k}, 'omitnan'); end

yl = ylim;
y_top_vals = yl(2) - 0.08*(yl(2)-yl(1));
y_top_text = yl(2) - 0.01*(yl(2)-yl(1));

for k = 1:4
    if ~isnan(med_vals(k))
        text(k, y_top_vals, sprintf('%.2f%%', med_vals(k)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
             'FontSize', 16, 'Color', 'k', 'Clipping', 'off');
    end
end
text(2.5, y_top_text, 'Median PCE', 'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'top', 'FontSize', 14, 'Clipping', 'off');
hold off;
end

%% ========================================================================
%  FIGURE 9: Correlation Matrix
%% ========================================================================
if ismember(9, figures_to_plot)
figure('Name', 'Correlation_Matrix', 'Position', [400 200 700 600], 'Color', 'w');

vars = {'SF', '\eta_{SF}', 'log(SRV)', 'E_F', '\Phi_R off', 'PCE', 'V_{oc}', 'J_{sc}', 'FF'};
data_matrix = [T_success.SF_flag, T_success.eta_SF, ...
               log10(T_success.SRV_cm_s), T_success.EF_emitter_eV, T_success.Phi_R_offset_eV, ...
               T_success.PCE_pct, T_success.Voc_V, T_success.Jsc_mA_cm2, T_success.FF_pct];

R = corrcoef(data_matrix, 'Rows', 'complete');
imagesc(R);
colormap(redgreen(256));
caxis([-1 1]);
cb = colorbar; ylabel(cb, 'Correlation', 'FontSize', 12);
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
top_idx = sort_idx(1:min(n_select, height(T_success)));
bottom_idx = sort_idx(max(1, end-n_select+1):end);

% No junction depth axis (only 1 value), replace with Phi_R offset
param_names = {'SF', '\eta_{SF}', 'SRV', 'E_F (eV)', '\Phi_R off'};
data_norm = zeros(height(T_success), 5);

data_norm(:, 1) = T_success.SF_flag;
data_norm(:, 2) = T_success.eta_SF / 2;
data_norm(:, 3) = (log10(T_success.SRV_cm_s) - 3) / 2;  % log10(1e3)=3, log10(1e5)=5
data_norm(:, 4) = (T_success.EF_emitter_eV - (-4.35)) / 0.3;
data_norm(:, 5) = (T_success.Phi_R_offset_eV - (-0.06)) / 0.12;

hold on;
for i = bottom_idx'
    plot(1:5, data_norm(i, :), '-', 'Color', [1 0.5 0.5 0.3], 'LineWidth', 1);
end
for i = top_idx'
    plot(1:5, data_norm(i, :), '-', 'Color', [0 0.6 0 0.6], 'LineWidth', 1.5);
end

h1 = plot(NaN, NaN, '-', 'Color', [0 0.6 0], 'LineWidth', 2);
h2 = plot(NaN, NaN, '-', 'Color', [1 0.5 0.5], 'LineWidth', 2);

set(gca, 'XTick', 1:5, 'XTickLabel', param_names);
xlim([0.5 5.5]); ylim([-0.1 1.1]);
ylabel('Normalized Value', 'FontSize', 20);
legend([h1 h2], {'Top 10% PCE', 'Bottom 10% PCE'}, 'Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 18, 'LineWidth', 1.2, 'Box', 'on');
grid on;
hold off;
end

%% ========================================================================
%  FIGURE 11: Scatter Matrix
%% ========================================================================
if ismember(11, figures_to_plot)
figure('Name', 'Scatter_Matrix', 'Position', [100 50 1000 600], 'Color', 'w');

subplot(1, 3, 1);
scatter(T_success.eta_SF, T_success.PCE_pct, 40, T_success.SF_flag, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('\eta_{SF}', 'FontSize', 14); ylabel('PCE (%)', 'FontSize', 14);
title('PCE vs \eta_{SF}', 'FontSize', 14);
colormap(gca, [0 0 0.8; 0 0.6 0]);
set(gca, 'FontSize', 12, 'Box', 'on'); grid on;

subplot(1, 3, 2);
scatter(log10(T_success.SRV_cm_s), T_success.PCE_pct, 40, T_success.eta_SF, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('log_{10}(SRV)', 'FontSize', 14); ylabel('PCE (%)', 'FontSize', 14);
title('PCE vs SRV', 'FontSize', 14);
colormap(gca, redgreen(256));
cb = colorbar; ylabel(cb, '\eta_{SF}');
set(gca, 'FontSize', 12, 'Box', 'on'); grid on;

subplot(1, 3, 3);
scatter(T_success.EF_emitter_eV, T_success.PCE_pct, 40, log10(T_success.SRV_cm_s), 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('E_F (eV)', 'FontSize', 14); ylabel('PCE (%)', 'FontSize', 14);
title('PCE vs Fermi Level', 'FontSize', 14);
colormap(gca, redgreen(256));
cb = colorbar; ylabel(cb, 'log_{10}(SRV)');
set(gca, 'FontSize', 12, 'Box', 'on'); grid on;
end

%% ========================================================================
%  FIGURE 12: Delta PCE Heatmap (SF benefit = PCE_SF - PCE_baseline)
%% ========================================================================
if ismember(12, figures_to_plot)
figure('Name', 'Delta_PCE_Heatmap', 'Position', [150 100 600 450], 'Color', 'w');

EF_vals_sweep = unique(T_success.EF_emitter_eV);
SRV_vals_sweep = unique(T_success.SRV_cm_s);

ND_axis = NaN(size(EF_vals_sweep));
for ii = 1:numel(EF_vals_sweep)
    [~, idx_map] = min(abs(EF_map - EF_vals_sweep(ii)));
    ND_axis(ii) = ND_map(idx_map);
end

delta_PCE = NaN(length(EF_vals_sweep), length(SRV_vals_sweep));

for i = 1:length(EF_vals_sweep)
    for j = 1:length(SRV_vals_sweep)
        mask_base = (T_success.SF_flag == 0) & ...
                    (T_success.junction_depth_nm == optimal_jd) & ...
                    (T_success.SRV_cm_s == SRV_vals_sweep(j)) & ...
                    (abs(T_success.EF_emitter_eV - EF_vals_sweep(i)) < tol) & ...
                    (abs(T_success.Phi_R_offset_eV - 0) < tol);

        mask_sf = (T_success.SF_flag == 1) & ...
                  (abs(T_success.eta_SF - 2) < tol) & ...
                  (T_success.junction_depth_nm == optimal_jd) & ...
                  (T_success.SRV_cm_s == SRV_vals_sweep(j)) & ...
                  (abs(T_success.EF_emitter_eV - EF_vals_sweep(i)) < tol) & ...
                  (abs(T_success.Phi_R_offset_eV - 0) < tol);

        if any(mask_base) && any(mask_sf)
            PCE_base = mean(T_success.PCE_pct(mask_base), 'omitnan');
            PCE_sf   = mean(T_success.PCE_pct(mask_sf),   'omitnan');
            delta_PCE(i, j) = PCE_sf - PCE_base;
        end
    end
end

h = imagesc(1:length(SRV_vals_sweep), 1:length(EF_vals_sweep), delta_PCE);
colormap(redgreen(256));
set(h, 'AlphaData', ~isnan(delta_PCE));
set(gca, 'Color', 'w');

cb = colorbar;
ylabel(cb, '\Delta PCE (%) [SF - Baseline]', 'FontSize', 18);
set(gca, 'YDir', 'normal');

set(gca, 'XTick', 1:length(SRV_vals_sweep), ...
         'XTickLabel', arrayfun(@(x) sprintf('%.0e', x), SRV_vals_sweep, 'UniformOutput', false));
set(gca, 'YTick', 1:length(EF_vals_sweep), ...
         'YTickLabel', arrayfun(@(x) sprintf('%.1e', x), ND_axis, 'UniformOutput', false));

xlabel('Surface recombination velocity s_p (cm/s)', 'FontSize', 20);
ylabel('n-emitter doping concentration N_D (cm^{-3})', 'FontSize', 20);
set(gca, 'FontSize', 16);

for i = 1:length(EF_vals_sweep)
    for j = 1:length(SRV_vals_sweep)
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
%  FIGURE 13: Radar Chart - Best vs Worst
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

best_norm = [(T_success.PCE_pct(best_idx) - PCE_range(1)) / max(diff(PCE_range), 1e-6), ...
             (T_success.Voc_V(best_idx) - Voc_range(1)) / max(diff(Voc_range), 1e-6), ...
             (T_success.Jsc_mA_cm2(best_idx) - Jsc_range(1)) / max(diff(Jsc_range), 1e-6), ...
             (T_success.FF_pct(best_idx) - FF_range(1)) / max(diff(FF_range), 1e-6), ...
             (best_blue_eqe - BlueEQE_range(1)) / diff(BlueEQE_range)];

worst_norm = [(T_success.PCE_pct(worst_idx) - PCE_range(1)) / max(diff(PCE_range), 1e-6), ...
              (T_success.Voc_V(worst_idx) - Voc_range(1)) / max(diff(Voc_range), 1e-6), ...
              (T_success.Jsc_mA_cm2(worst_idx) - Jsc_range(1)) / max(diff(Jsc_range), 1e-6), ...
              (T_success.FF_pct(worst_idx) - FF_range(1)) / max(diff(FF_range), 1e-6), ...
              (worst_blue_eqe - BlueEQE_range(1)) / diff(BlueEQE_range)];

best_norm = [best_norm, best_norm(1)];
worst_norm = [worst_norm, worst_norm(1)];
angles = linspace(0, 2*pi, n_metrics + 1);

polarplot(angles, best_norm, 'g-', 'LineWidth', 2.5, 'DisplayName', 'Best config');
hold on;
polarplot(angles, worst_norm, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Worst config');
ax = gca;
ax.ThetaTick = rad2deg(angles(1:end-1));
ax.ThetaTickLabel = metrics;
ax.RLim = [0 1.1];
ax.FontSize = 12;
title('Best vs Worst Configuration', 'FontSize', 14);
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

    area(wl, eqe_base, 'FaceColor', [0.3 0.3 0.8], 'FaceAlpha', 0.6, ...
         'EdgeColor', 'none', 'DisplayName', 'Baseline Si');
    hold on;
    eqe_diff = max(0, eqe_sf - eqe_base);
    area(wl, eqe_base + eqe_diff, 'FaceColor', [0.2 0.7 0.2], 'FaceAlpha', 0.6, ...
         'EdgeColor', 'none', 'DisplayName', 'SF Enhancement');
    plot(wl, eqe_sf, 'g-', 'LineWidth', 2, 'DisplayName', 'Total (SF, \eta=2)');
    plot([400 800], [100 100], 'k--', 'LineWidth', 1.5, 'DisplayName', '100% EQE');

    xlabel('Wavelength (nm)', 'FontSize', 16);
    ylabel('EQE (%)', 'FontSize', 16);
    title('EQE Decomposition: Baseline + SF Enhancement', 'FontSize', 18);
    xlim([400 700]); ylim([0 150]);
    legend('Location', 'northeast', 'FontSize', 12);
    set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'on');
    grid on; hold off;
end
end

%% ========================================================================
%  FIGURE 15: Summary Table
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
    sprintf('  Mean PCE (SF=1, eta=2): %.3f %%', ...
            mean(T_success.PCE_pct(T_success.SF_flag==1 & T_success.eta_SF==2), 'omitnan'))
};

text(0.1, 0.95, summary_text, 'VerticalAlignment', 'top', 'FontSize', 11, 'FontName', 'FixedWidth');
title('Parameter Sweep Summary (Reverse Bias)', 'FontSize', 18);
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

plot(wl, eqe, 'b-', 'LineWidth', 2.5, 'DisplayName', 'EQE');
hold on;
plot(wl, iqe, 'r--', 'LineWidth', 2, 'DisplayName', 'IQE');
plot([400 800], [100 100], 'k--', 'LineWidth', 1.5, 'DisplayName', '100%');

fill([400 550 550 400], [0 0 160 160], [0.9 0.95 1], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off');

idx_above = eqe > 100;
if any(idx_above)
    wl_above = wl(idx_above);
    eqe_above = eqe(idx_above);
    fill([wl_above, fliplr(wl_above)], [100*ones(size(wl_above)), fliplr(eqe_above)], ...
         [0.8 1 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
end

xlabel('Wavelength (nm)', 'FontSize', 18);
ylabel('Quantum Efficiency (%)', 'FontSize', 18);
title(sprintf('Best Device EQE (PCE = %.3f%%)', best.PCE_pct), 'FontSize', 20);
xlim([400 800]); ylim([0 160]);
legend('Location', 'northeast', 'FontSize', 14);
set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'on');
grid on;

[max_eqe, max_idx] = max(eqe);
text(wl(max_idx), max_eqe + 5, sprintf('Peak: %.1f%%', max_eqe), ...
     'FontSize', 12, 'HorizontalAlignment', 'center', 'Color', [0 0 0.7]);

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

[~, best_idx] = max(T_success.PCE_pct);
best = T_success(best_idx, :);
best_run = runs(best.run_number);

V = best_run.JV.Vapp;
J = best_run.JV.Jtot;
[~, idx_max] = max(V);
V_fwd = V(1:idx_max);
J_fwd = J(1:idx_max);

plot(V_fwd, J_fwd, 'b-', 'LineWidth', 2.5);
hold on;
plot(xlim, [0 0], 'k--', 'LineWidth', 1);

[~, idx_v0] = min(abs(V_fwd));
Jsc_plot = J_fwd(idx_v0);
plot(0, Jsc_plot, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
text(0.02, Jsc_plot, sprintf('J_{sc} = %.4f mA/cm^2', abs(Jsc_plot)), ...
     'FontSize', 11, 'VerticalAlignment', 'bottom');

Voc_plot = best.Voc_V;
plot(Voc_plot, 0, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
text(Voc_plot, 0.02, sprintf('V_{oc} = %.4f V', Voc_plot), ...
     'FontSize', 11, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

power = V_fwd .* J_fwd;
[MPP, idx_mpp] = min(power);
V_mpp = V_fwd(idx_mpp(1));
J_mpp = J_fwd(idx_mpp(1));
plot(V_mpp, J_mpp, 'ms', 'MarkerSize', 12, 'MarkerFaceColor', 'm');
text(V_mpp + 0.02, J_mpp, sprintf('MPP\nP = %.4f mW/cm^2', abs(MPP(1))), ...
     'FontSize', 10, 'VerticalAlignment', 'middle');

fill([0 V_mpp V_mpp 0], [0 0 J_mpp J_mpp], [1 0.9 0.9], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.3);

xlabel('Voltage (V)', 'FontSize', 18);
ylabel('Current Density (mA/cm^2)', 'FontSize', 18);
title(sprintf('Best Device J-V (PCE = %.3f%%, FF = %.2f%%)', best.PCE_pct, best.FF_pct), 'FontSize', 20);
set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'on');
grid on;

metrics_str = sprintf(['PCE = %.3f %%\nV_{oc} = %.4f V\n' ...
                       'J_{sc} = %.4f mA/cm^2\nFF = %.2f %%'], ...
                      best.PCE_pct, best.Voc_V, best.Jsc_mA_cm2, best.FF_pct);
annotation('textbox', [0.15, 0.15, 0.25, 0.2], 'String', metrics_str, ...
           'FontSize', 12, 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
           'FitBoxToText', 'on', 'Interpreter', 'tex');
hold off;
end

%% ========================================================================
%  FIGURE 18: Best Device vs Baseline
%% ========================================================================
if ismember(18, figures_to_plot)
figure('Name', 'Best_vs_Baseline', 'Position', [100 100 1400 500], 'Color', 'w');

[~, best_idx] = max(T_success.PCE_pct);
best = T_success(best_idx,:);
best_run = runs(best.run_number);

mask_baseline = (T_success.SF_flag == 0) & ...
                (T_success.junction_depth_nm == best.junction_depth_nm) & ...
                (T_success.SRV_cm_s == best.SRV_cm_s) & ...
                (abs(T_success.EF_emitter_eV - best.EF_emitter_eV) < tol) & ...
                (abs(T_success.Phi_R_offset_eV - best.Phi_R_offset_eV) < tol);

if any(mask_baseline)
    baseline = T_success(find(mask_baseline, 1), :);
    baseline_run = runs(baseline.run_number);

    subplot(1, 2, 1); hold on;
    h_sf_eqe = plot(best_run.EQE.wavelengths, best_run.EQE.EQE, '-', 'Color', [0 0.6 0], 'LineWidth', 2.5);
    h_bl_eqe = plot(baseline_run.EQE.wavelengths, baseline_run.EQE.EQE, 'b-', 'LineWidth', 2.5);
    plot([400 800], [100 100], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlabel('Wavelength (nm)', 'FontSize', 24);
    ylabel('EQE (%)', 'FontSize', 24);
    xlim([400 800]); ylim([80 110]);
    set(gca, 'FontSize', 18, 'LineWidth', 1.2, 'Box', 'on'); grid on;
    [eqe_max, idx_max_EQE] = max(best_run.EQE.EQE);
    lambda_max = best_run.EQE.wavelengths(idx_max_EQE);
    plot(lambda_max, eqe_max, 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0 0.6 0], ...
         'MarkerEdgeColor', [0 0.3 0], 'HandleVisibility', 'off');
    yl = ylim; y_offset = 0.02 * (yl(2) - yl(1));
    text(lambda_max + 5, eqe_max + y_offset, sprintf('%.1f%% @ %dnm', eqe_max, round(lambda_max)), ...
         'FontSize', 18, 'Color', [0 0.4 0], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    legend([h_bl_eqe, h_sf_eqe], {'Baseline', sprintf('SF (\\eta=%.0f)', best.eta_SF)}, ...
           'Location', 'northeast', 'FontSize', 16);
    hold off;

    subplot(1, 2, 2); hold on;
    V_best = best_run.JV.Vapp; J_best = best_run.JV.Jtot;
    [~, idx_max_best] = max(V_best);
    h_sf = plot(V_best(1:idx_max_best), J_best(1:idx_max_best), '-', 'Color', [0 0.6 0], 'LineWidth', 2.5);
    V_base = baseline_run.JV.Vapp; J_base = baseline_run.JV.Jtot;
    [~, idx_max_base] = max(V_base);
    h_bl = plot(V_base(1:idx_max_base), J_base(1:idx_max_base), 'b-', 'LineWidth', 2.5);
    plot(xlim, [0 0], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlabel('Voltage (V)', 'FontSize', 24);
    ylabel('Current Density (mA/cm^2)', 'FontSize', 24);
    xlim([-0.3, 0.8]); ylim([-0.05, 0.1]);
    set(gca, 'FontSize', 18, 'LineWidth', 1.2, 'Box', 'on'); grid on;

    base_line1 = sprintf('Baseline: PCE = %.3f%%   V_{oc} = %.4f V', baseline.PCE_pct, baseline.Voc_V);
    base_line2 = sprintf('         J_{sc} = %.4f mA/cm^2   FF = %.2f%%', baseline.Jsc_mA_cm2, baseline.FF_pct);
    sf_line1 = sprintf('SF (\\eta=%.0f): PCE = %.3f%%   V_{oc} = %.4f V', best.eta_SF, best.PCE_pct, best.Voc_V);
    sf_line2 = sprintf('         J_{sc} = %.4f mA/cm^2   FF = %.2f%%', best.Jsc_mA_cm2, best.FF_pct);
    h_d1 = plot(NaN, NaN, 'w', 'LineStyle', 'none', 'Marker', 'none');
    h_d2 = plot(NaN, NaN, 'w', 'LineStyle', 'none', 'Marker', 'none');
    legend([h_bl, h_d1, h_sf, h_d2], {base_line1, base_line2, sf_line1, sf_line2}, ...
           'Location', 'northwest', 'FontSize', 16, 'Interpreter', 'tex');
    hold off;
else
    text(0.5, 0.5, 'No matching baseline found', 'FontSize', 16, 'HorizontalAlignment', 'center');
end
end

%% ========================================================================
%  FIGURE 19: SR_EQE Under Bias — Effect of Bias Voltage (single device)
%  Shows how EQE changes at different reverse-bias points for the best SF device
%% ========================================================================
if ismember(19, figures_to_plot)
figure('Name', 'SR_EQE_vs_Bias', 'Position', [50 50 800 550], 'Color', 'w');

% Use best SF device (eta=2, best SRV, optimal EF)
mask = (T_success.SF_flag == 1) & ...
       (abs(T_success.eta_SF - 2) < tol) & ...
       (T_success.SRV_cm_s == optimal_SRV) & ...
       (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
       (T_success.junction_depth_nm == optimal_jd) & ...
       (abs(T_success.Phi_R_offset_eV - 0) < tol);

idx = find(mask, 1);
if ~isempty(idx)
    run_idx = T_success.run_number(idx);
    r = runs(run_idx);

    wl = r.EQE.wavelengths;
    hold on;

    % Check if SR data exists in this run
    % Postprocess stores: r.SR.bias_vals, r.SR.EQE, r.SR.J
    if isfield(r, 'SR') && isfield(r.SR, 'bias_vals') && ~isempty(r.SR.bias_vals)
        bias_vals = r.SR.bias_vals;
        SR_EQE = r.SR.EQE;        % [n_wavelengths x n_bias]
    else
        fprintf('  Warning: No SR data found in run structure. Plotting standard EQE only.\n');
        bias_vals = 0;
        SR_EQE = r.EQE.EQE(:);
    end

    n_bias = length(bias_vals);

    % Colors: most negative bias = dark blue, V=0 = green
    bias_colors = [0.00 0.00 0.70;   % -1.0 V  - dark blue
                   0.40 0.40 0.80;   % -0.5 V  - medium blue
                   0.00 0.60 0.00];  % 0.0 V   - green
    if n_bias ~= 3
        bias_colors = lines(n_bias);
    end

    h_bias = gobjects(n_bias, 1);
    for ib = 1:n_bias
        if size(SR_EQE, 2) >= ib
            eqe_at_bias = SR_EQE(:, ib);
        else
            eqe_at_bias = SR_EQE(:);
        end
        h_bias(ib) = plot(wl, eqe_at_bias, '-', 'LineWidth', 2.5, ...
                          'Color', bias_colors(min(ib, size(bias_colors,1)), :));
    end

    plot([400 800], [100 100], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

    xlabel('Wavelength (nm)', 'FontSize', 20);
    ylabel('EQE (%)', 'FontSize', 20);
    title('Spectral Response Under Reverse Bias (SF, \eta=2)', 'FontSize', 20);
    xlim([400 800]);
    set(gca, 'FontSize', 16, 'LineWidth', 1.2, 'Box', 'on');
    grid on;

    leg_labels = arrayfun(@(v) sprintf('V_{bias} = %.1f V', v), bias_vals, 'UniformOutput', false);
    valid = isgraphics(h_bias);
    legend(h_bias(valid), leg_labels(valid), 'Location', 'northeast', 'FontSize', 14);

    hold off;
else
    text(0.5, 0.5, 'No matching SF run found', 'FontSize', 16, 'HorizontalAlignment', 'center');
end
end

%% ========================================================================
%  FIGURE 20: SR_EQE Under Bias — SF vs Baseline at Each Bias Point
%  Multi-panel: one panel per bias voltage comparing SF (eta=2) vs baseline
%% ========================================================================
if ismember(20, figures_to_plot)

% Determine number of bias points from first successful run with SR data
sample_mask = (T_success.SF_flag == 1) & (abs(T_success.eta_SF - 2) < tol);
sample_idx = find(sample_mask, 1);
if ~isempty(sample_idx)
    sample_run = runs(T_success.run_number(sample_idx));
    if isfield(sample_run, 'SR') && isfield(sample_run.SR, 'bias_vals') && ~isempty(sample_run.SR.bias_vals)
        bias_vals = sample_run.SR.bias_vals;
    else
        bias_vals = 0;
    end
else
    bias_vals = 0;
end

n_bias = length(bias_vals);
figure('Name', 'SR_EQE_SF_vs_Baseline', ...
       'Position', [50 50 420*n_bias 400], 'Color', 'w');

for ib = 1:n_bias
    subplot(1, n_bias, ib);
    hold on;

    % SF device (eta=2)
    mask_sf = (T_success.SF_flag == 1) & ...
              (abs(T_success.eta_SF - 2) < tol) & ...
              (T_success.SRV_cm_s == optimal_SRV) & ...
              (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
              (T_success.junction_depth_nm == optimal_jd) & ...
              (abs(T_success.Phi_R_offset_eV - 0) < tol);

    idx_sf = find(mask_sf, 1);
    if ~isempty(idx_sf)
        r_sf = runs(T_success.run_number(idx_sf));
        wl = r_sf.EQE.wavelengths;

        if isfield(r_sf, 'SR') && isfield(r_sf.SR, 'EQE') && ~isempty(r_sf.SR.EQE)
            sr_eqe_sf = r_sf.SR.EQE(:, ib);
        else
            sr_eqe_sf = r_sf.EQE.EQE(:);
        end

        h_sf = plot(wl, sr_eqe_sf, '-', 'Color', [0 0.6 0], 'LineWidth', 2.5);
    end

    % Baseline
    mask_base = (T_success.SF_flag == 0) & ...
                (T_success.SRV_cm_s == optimal_SRV) & ...
                (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
                (T_success.junction_depth_nm == optimal_jd) & ...
                (abs(T_success.Phi_R_offset_eV - 0) < tol);

    idx_base = find(mask_base, 1);
    if ~isempty(idx_base)
        r_base = runs(T_success.run_number(idx_base));

        if isfield(r_base, 'SR') && isfield(r_base.SR, 'EQE') && ~isempty(r_base.SR.EQE)
            sr_eqe_base = r_base.SR.EQE(:, ib);
        else
            sr_eqe_base = r_base.EQE.EQE(:);
        end

        h_bl = plot(wl, sr_eqe_base, 'b-', 'LineWidth', 2.5);
    end

    plot([400 800], [100 100], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

    xlabel('Wavelength (nm)', 'FontSize', 18);
    if ib == 1
        ylabel('EQE (%)', 'FontSize', 18);
    end
    title(sprintf('V_{bias} = %.1f V', bias_vals(ib)), 'FontSize', 18);
    xlim([400 800]); ylim([50 130]);
    set(gca, 'FontSize', 14, 'Box', 'on');
    grid on;

    if ib == n_bias && exist('h_sf', 'var') && exist('h_bl', 'var')
        legend([h_bl, h_sf], {'Baseline', 'SF (\eta=2)'}, ...
               'Location', 'southeast', 'FontSize', 12);
    end

    hold off;
end
end

%% ========================================================================
%  FIGURE 21: SR_EQE Ratio — Reverse-Bias EQE / Short-Circuit EQE
%  Shows how much additional collection occurs under reverse bias
%  (useful for identifying field-dependent collection losses)
%% ========================================================================
if ismember(21, figures_to_plot)
figure('Name', 'SR_EQE_Ratio', 'Position', [100 100 800 500], 'Color', 'w');

% Compare SF (eta=2) and baseline at optimal params
configs = {'SF (\eta=2)', 'Baseline'};
config_colors = {[0 0.6 0], [0.27 0.33 0.80]};
config_styles = {'-', '--'};

for ic = 1:2
    if ic == 1
        mask = (T_success.SF_flag == 1) & (abs(T_success.eta_SF - 2) < tol);
    else
        mask = (T_success.SF_flag == 0);
    end
    mask = mask & (T_success.SRV_cm_s == optimal_SRV) & ...
           (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
           (T_success.junction_depth_nm == optimal_jd) & ...
           (abs(T_success.Phi_R_offset_eV - 0) < tol);

    idx = find(mask, 1);
    if ~isempty(idx)
        r = runs(T_success.run_number(idx));
        wl = r.EQE.wavelengths;

        if isfield(r, 'SR') && isfield(r.SR, 'EQE') && ~isempty(r.SR.EQE)
            sr_eqe = r.SR.EQE;
            bias_vals_run = r.SR.bias_vals;
        else
            continue;
        end

        % Find V=0 column (reference) and most negative bias column
        [~, idx_0V] = min(abs(bias_vals_run - 0));
        [~, idx_neg] = min(bias_vals_run);  % most negative

        eqe_0V  = sr_eqe(:, idx_0V);
        eqe_neg = sr_eqe(:, idx_neg);

        % Ratio: reverse-bias / short-circuit
        ratio = eqe_neg ./ max(eqe_0V, 0.1);   % avoid div by zero

        hold on;
        plot(wl, ratio, config_styles{ic}, 'LineWidth', 2.5, ...
             'Color', config_colors{ic}, ...
             'DisplayName', sprintf('%s (V=%.1fV / V=0)', configs{ic}, bias_vals_run(idx_neg)));
    end
end

plot([400 800], [1 1], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

xlabel('Wavelength (nm)', 'FontSize', 20);
ylabel('EQE(V_{rev}) / EQE(V=0)', 'FontSize', 20);
title('Reverse-Bias Collection Enhancement', 'FontSize', 20);
xlim([400 800]);
set(gca, 'FontSize', 16, 'LineWidth', 1.2, 'Box', 'on');
grid on;
legend('Location', 'best', 'FontSize', 14);
hold off;
end

%% ========================================================================
%  FIGURE 22: Heatmap — Integrated Blue SR_EQE vs (SRV, Bias)
%  Shows how blue-region EQE changes with SRV and applied bias
%% ========================================================================
if ismember(22, figures_to_plot)
figure('Name', 'Blue_SREQE_Heatmap', 'Position', [150 150 600 400], 'Color', 'w');

% Get wavelength info from a sample run
sample_mask = strcmp({runs.status}, 'success');
sample_idx = find(sample_mask, 1);
wl = runs(sample_idx).EQE.wavelengths;
blue_mask_wl = (wl >= 405) & (wl <= 500);

% Determine bias values
r_sample = runs(sample_idx);
if isfield(r_sample, 'SR') && isfield(r_sample.SR, 'bias_vals') && ~isempty(r_sample.SR.bias_vals)
    bias_vals_plot = r_sample.SR.bias_vals;
else
    bias_vals_plot = 0;
end
n_bias = length(bias_vals_plot);

% Build matrix: rows = SRV, cols = bias
blue_SR_EQE = NaN(length(SRV_list), n_bias);

for k = 1:length(SRV_list)
    mask = (T_success.SF_flag == 1) & ...
           (abs(T_success.eta_SF - 2) < tol) & ...
           (T_success.SRV_cm_s == SRV_list(k)) & ...
           (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
           (T_success.junction_depth_nm == optimal_jd) & ...
           (abs(T_success.Phi_R_offset_eV - 0) < tol);
    idx = find(mask, 1);
    if ~isempty(idx)
        r = runs(T_success.run_number(idx));

        if isfield(r, 'SR') && isfield(r.SR, 'EQE') && ~isempty(r.SR.EQE)
            sr_eqe = r.SR.EQE;
        else
            sr_eqe = repmat(r.EQE.EQE(:), 1, n_bias);
        end

        for ib = 1:n_bias
            blue_SR_EQE(k, ib) = mean(sr_eqe(blue_mask_wl, ib), 'omitnan');
        end
    end
end

h = imagesc(1:n_bias, 1:length(SRV_list), blue_SR_EQE);
colormap(redgreen(256));
set(h, 'AlphaData', ~isnan(blue_SR_EQE));
set(gca, 'Color', 'w');

cb = colorbar;
ylabel(cb, 'Mean Blue EQE (405-500nm) (%)', 'FontSize', 14);
set(gca, 'YDir', 'normal');

set(gca, 'XTick', 1:n_bias, ...
         'XTickLabel', arrayfun(@(v) sprintf('%.1f V', v), bias_vals_plot, 'UniformOutput', false));
set(gca, 'YTick', 1:length(SRV_list), ...
         'YTickLabel', arrayfun(@(x) sprintf('%.0e', x), SRV_list, 'UniformOutput', false));

xlabel('Applied Bias (V)', 'FontSize', 18);
ylabel('Surface recombination velocity s_p (cm/s)', 'FontSize', 18);
title('Blue-Region SR-EQE vs Bias and SRV (SF, \eta=2)', 'FontSize', 18);
set(gca, 'FontSize', 14, 'LineWidth', 1.2);

for k = 1:length(SRV_list)
    for ib = 1:n_bias
        if ~isnan(blue_SR_EQE(k, ib))
            text(ib, k, sprintf('%.1f', blue_SR_EQE(k, ib)), ...
                 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
        end
    end
end
end

%% ========================================================================
%  FIGURE 23: SR_EQE Under Bias — Effect of SRV at Fixed Reverse Bias
%  Multi-panel: one per bias point, lines for each SRV (SF eta=2)
%% ========================================================================
if ismember(23, figures_to_plot)

% Get bias values from a sample run
sample_mask = (T_success.SF_flag == 1) & (abs(T_success.eta_SF - 2) < tol);
sample_idx = find(sample_mask, 1);
if ~isempty(sample_idx)
    r_sample = runs(T_success.run_number(sample_idx));
    if isfield(r_sample, 'SR') && isfield(r_sample.SR, 'bias_vals') && ~isempty(r_sample.SR.bias_vals)
        bias_vals_fig = r_sample.SR.bias_vals;
    else
        bias_vals_fig = [];
    end
else
    bias_vals_fig = [];
end

if ~isempty(bias_vals_fig)
n_bias = length(bias_vals_fig);
figure('Name', 'SR_EQE_SRV_per_bias', ...
       'Position', [50 50 420*n_bias 400], 'Color', 'w');

SRV_colors_23 = [0.00 0.45 0.00; 0.90 0.80 0.20; 0.80 0.00 0.00];

for ib = 1:n_bias
    subplot(1, n_bias, ib);
    hold on;

    h_srv23 = gobjects(length(SRV_list), 1);
    for s = 1:length(SRV_list)
        mask = (T_success.SF_flag == 1) & ...
               (abs(T_success.eta_SF - 2) < tol) & ...
               (T_success.SRV_cm_s == SRV_list(s)) & ...
               (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
               (T_success.junction_depth_nm == optimal_jd) & ...
               (abs(T_success.Phi_R_offset_eV - 0) < tol);
        idx = find(mask, 1);
        if ~isempty(idx)
            r = runs(T_success.run_number(idx));
            wl = r.EQE.wavelengths;
            if isfield(r, 'SR') && isfield(r.SR, 'EQE') && ~isempty(r.SR.EQE)
                sr_eqe = r.SR.EQE(:, ib);
            else
                sr_eqe = r.EQE.EQE(:);
            end
            h_srv23(s) = plot(wl, sr_eqe, '-', 'LineWidth', 2, ...
                              'Color', SRV_colors_23(min(s, size(SRV_colors_23,1)), :), ...
                              'DisplayName', sprintf('s_p = %.0e', SRV_list(s)));
        end
    end

    plot([400 800], [100 100], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlabel('Wavelength (nm)', 'FontSize', 18);
    if ib == 1, ylabel('EQE (%)', 'FontSize', 18); end
    title(sprintf('V_{bias} = %.1f V', bias_vals_fig(ib)), 'FontSize', 18);
    xlim([400 800]); ylim([50 130]);
    set(gca, 'FontSize', 14, 'Box', 'on'); grid on;

    if ib == n_bias
        valid = isgraphics(h_srv23);
        legend(h_srv23(valid), 'Location', 'southeast', 'FontSize', 12);
    end
    hold off;
end
else
    fprintf('Figure 23 skipped: no SR data found.\n');
end
end

%% ========================================================================
%  FIGURE 24: SR_EQE Under Bias — Effect of eta_SF at Fixed Reverse Bias
%  Multi-panel: one per bias point, lines for each eta (at optimal SRV)
%% ========================================================================
if ismember(24, figures_to_plot)

sample_mask = (T_success.SF_flag == 1) & (abs(T_success.eta_SF - 2) < tol);
sample_idx = find(sample_mask, 1);
if ~isempty(sample_idx)
    r_sample = runs(T_success.run_number(sample_idx));
    if isfield(r_sample, 'SR') && isfield(r_sample.SR, 'bias_vals') && ~isempty(r_sample.SR.bias_vals)
        bias_vals_fig = r_sample.SR.bias_vals;
    else
        bias_vals_fig = [];
    end
else
    bias_vals_fig = [];
end

if ~isempty(bias_vals_fig)
n_bias = length(bias_vals_fig);
figure('Name', 'SR_EQE_eta_per_bias', ...
       'Position', [50 100 420*n_bias 400], 'Color', 'w');

eta_vals_24 = [2, 1, 0];
eta_colors_24 = [0.00 0.60 0.00; 0.90 0.80 0.20; 0.80 0.00 0.00];

for ib = 1:n_bias
    subplot(1, n_bias, ib);
    hold on;

    h_eta24 = gobjects(numel(eta_vals_24), 1);
    for ie = 1:numel(eta_vals_24)
        eta_val = eta_vals_24(ie);
        mask = (T_success.SF_flag == 1) & ...
               (abs(T_success.eta_SF - eta_val) < tol) & ...
               (T_success.SRV_cm_s == optimal_SRV) & ...
               (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
               (T_success.junction_depth_nm == optimal_jd) & ...
               (abs(T_success.Phi_R_offset_eV - 0) < tol);
        idx = find(mask, 1);
        if ~isempty(idx)
            r = runs(T_success.run_number(idx));
            wl = r.EQE.wavelengths;
            if isfield(r, 'SR') && isfield(r.SR, 'EQE') && ~isempty(r.SR.EQE)
                sr_eqe = r.SR.EQE(:, ib);
            else
                sr_eqe = r.EQE.EQE(:);
            end
            h_eta24(ie) = plot(wl, sr_eqe, '-', 'LineWidth', 2, ...
                               'Color', eta_colors_24(ie, :));
        end
    end

    % Also plot baseline
    mask_base = (T_success.SF_flag == 0) & ...
                (T_success.SRV_cm_s == optimal_SRV) & ...
                (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
                (T_success.junction_depth_nm == optimal_jd) & ...
                (abs(T_success.Phi_R_offset_eV - 0) < tol);
    idx_base = find(mask_base, 1);
    if ~isempty(idx_base)
        r_base = runs(T_success.run_number(idx_base));
        if isfield(r_base, 'SR') && isfield(r_base.SR, 'EQE') && ~isempty(r_base.SR.EQE)
            sr_eqe_base = r_base.SR.EQE(:, ib);
        else
            sr_eqe_base = r_base.EQE.EQE(:);
        end
        h_base24 = plot(r_base.EQE.wavelengths, sr_eqe_base, 'b--', 'LineWidth', 2);
    end

    plot([400 800], [100 100], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlabel('Wavelength (nm)', 'FontSize', 18);
    if ib == 1, ylabel('EQE (%)', 'FontSize', 18); end
    title(sprintf('V_{bias} = %.1f V', bias_vals_fig(ib)), 'FontSize', 18);
    xlim([400 800]); ylim([50 130]);
    set(gca, 'FontSize', 14, 'Box', 'on'); grid on;

    if ib == n_bias
        leg_labels = arrayfun(@(v) sprintf('\\eta_{SF} = %d', v), eta_vals_24, 'UniformOutput', false);
        valid = isgraphics(h_eta24);
        all_h = h_eta24(valid);
        all_l = leg_labels(valid);
        if exist('h_base24', 'var') && isgraphics(h_base24)
            all_h = [all_h; h_base24];
            all_l = [all_l, {'Baseline'}];
        end
        legend(all_h, all_l, 'Location', 'southeast', 'FontSize', 12);
    end
    hold off;
end
else
    fprintf('Figure 24 skipped: no SR data found.\n');
end
end

%% ========================================================================
%  FIGURE 25: Delta EQE (reverse-bias minus short-circuit) — SF Gain Under Bias
%  Shows how much extra current is collected under reverse bias, decomposed
%  into baseline collection gain vs SF-specific gain
%% ========================================================================
if ismember(25, figures_to_plot)
figure('Name', 'DeltaEQE_ReverseBias', 'Position', [100 100 900 500], 'Color', 'w');

hold on;

% Find the most negative bias index
configs_25 = {
    struct('label', 'SF (\eta=2)', 'SF_flag', 1, 'eta', 2, 'color', [0 0.6 0], 'style', '-');
    struct('label', 'SF (\eta=1)', 'SF_flag', 1, 'eta', 1, 'color', [0.9 0.8 0.2], 'style', '-');
    struct('label', 'Baseline',    'SF_flag', 0, 'eta', 0, 'color', [0.27 0.33 0.8], 'style', '--');
};

h_25 = gobjects(length(configs_25), 1);

for ic = 1:length(configs_25)
    cfg = configs_25{ic};

    if cfg.SF_flag == 0
        mask = (T_success.SF_flag == 0);
    else
        mask = (T_success.SF_flag == 1) & (abs(T_success.eta_SF - cfg.eta) < tol);
    end
    mask = mask & (T_success.SRV_cm_s == optimal_SRV) & ...
           (abs(T_success.EF_emitter_eV - optimal_EF) < tol) & ...
           (T_success.junction_depth_nm == optimal_jd) & ...
           (abs(T_success.Phi_R_offset_eV - 0) < tol);

    idx = find(mask, 1);
    if ~isempty(idx)
        r = runs(T_success.run_number(idx));
        wl = r.EQE.wavelengths;

        if isfield(r, 'SR') && isfield(r.SR, 'EQE') && ~isempty(r.SR.EQE)
            bias_vals_run = r.SR.bias_vals;
            [~, idx_0V]  = min(abs(bias_vals_run - 0));
            [~, idx_neg] = min(bias_vals_run);

            eqe_0V  = r.SR.EQE(:, idx_0V);
            eqe_neg = r.SR.EQE(:, idx_neg);
            delta_eqe = eqe_neg - eqe_0V;

            h_25(ic) = plot(wl, delta_eqe, cfg.style, 'LineWidth', 2.5, ...
                            'Color', cfg.color, 'DisplayName', ...
                            sprintf('%s (V=%.1fV - V=0)', cfg.label, bias_vals_run(idx_neg)));
        end
    end
end

plot([400 800], [0 0], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

xlabel('Wavelength (nm)', 'FontSize', 20);
ylabel('\Delta EQE (%) [V_{rev} - V_{SC}]', 'FontSize', 20);
title('Additional Collection Under Reverse Bias', 'FontSize', 20);
xlim([400 800]);
set(gca, 'FontSize', 16, 'LineWidth', 1.2, 'Box', 'on');
grid on;
valid = isgraphics(h_25);
if any(valid)
    legend(h_25(valid), 'Location', 'best', 'FontSize', 14);
end
hold off;
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
    if nargin < 1, m = 256; end
    r = [linspace(0, 1, m/2), ones(1, m/2)];
    g = [linspace(0, 1, m/2), linspace(1, 0, m/2)];
    b = [ones(1, m/2), linspace(1, 0, m/2)];
    c = [r(:), g(:), b(:)];
end

function c = redgreen(m)
    if nargin < 1, m = 256; end
    n1 = floor(m/2); n2 = m - n1;
    r1 = linspace(0.65, 1.0, n1); g1 = linspace(0.00, 1.0, n1); b1 = linspace(0.00, 1.0, n1);
    r2 = linspace(1.0, 0.00, n2); g2 = linspace(1.0, 0.50, n2); b2 = linspace(1.0, 0.00, n2);
    c = [[r1, r2]', [g1, g2]', [b1, b2]'];
end

function mEQE = local_mean_eqe_at_lambda(runs, run_idx_list, lambda_target)
    eqe_vals = NaN(numel(run_idx_list), 1);
    for ii = 1:numel(run_idx_list)
        r = runs(run_idx_list(ii));
        wl = r.EQE.wavelengths;
        [~, idx] = min(abs(wl - lambda_target));
        eqe_vals(ii) = r.EQE.EQE(idx);
    end
    mEQE = mean(eqe_vals, 'omitnan');
end