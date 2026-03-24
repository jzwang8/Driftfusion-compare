function generate_poster_figures_profiles(summaryFile)
%GENERATE_POSTER_FIGURES_PROFILES Create figures for doping profile sweep analysis
%
%   generate_poster_figures_profiles('sweep_summary_profiles_20260227.mat')
%
%   Profile-comparison figures:
%     1  - Grouped bar: PCE by profile type at each doping level
%     4  - EQE comparison: 3 profiles overlaid (best-case params)
%     5  - Multi-panel EQE: profile effect at different SRV
%     6  - Line plot: PCE vs SRV for each profile (with and without SF)
%     7  - Scatter: PCE(exponential) vs PCE(uniform) parity plot
%     10 - Box plot: PCE distribution by profile type
%     11 - Line plot: PCE vs N_D for each profile (SF=0 and SF=1,eta=2)
%     12 - Multi-panel EQE: profile effect at different N_D
%     13 - Multi-panel EQE: Qf effect by profile type (SF-off vs SF-on)
%
%   Uniform-profile figures (adapted from main/main_baseline):
%     14 - EQE comparison: SF signature (EQE > 100%)
%     15 - Heatmap: PCE vs (SRV, N_D)
%     16 - Line plot: SF benefit vs SRV
%     17 - Blue-region EQE vs eta_SF
%     18 - Multi-panel EQE: effect of eta, SRV, j_d, N_D (with baselines)
%     19 - EQE vs Phi_R offset (with baseline)
%     20 - Box plot: PCE distribution by eta_SF
%     21 - Correlation matrix: parameters vs metrics
%     22 - Scatter matrix: PCE vs key parameters
%     23 - Delta PCE heatmap (SF - baseline)
%     24 - Radar chart: best vs worst config
%     25 - Stacked area: EQE decomposition
%     26 - Summary table
%     27 - Best device EQE
%     28 - Best device J-V
%     29 - Best device vs baseline comparison (EQE + J-V)
%     30 - Line plot: PCE vs junction depth (SF vs baseline)
%     31 - Line plot: EQE(520nm) vs junction depth (SF vs baseline)

close all;

%% ========================================================================
%  CONFIGURATION
%% ========================================================================

FIGURES_TO_PLOT = 'all';
% FIGURES_TO_PLOT = [1:5];

if ischar(FIGURES_TO_PLOT) && strcmp(FIGURES_TO_PLOT, 'all')
    figures_to_plot = [1:13, 14:31];
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

success_mask = strcmp(T.status, 'success');
T_success = T(success_mask, :);

%% Shared settings
tol = 0.01;
n_prof = 3;
profile_types = {'uniform', 'gaussian', 'exponential'};
profile_labels = {'Uniform', 'Gaussian', 'Exponential'};
profile_colors = [0.27 0.33 0.80;    % blue  - uniform
                  0.90 0.60 0.00;    % orange - gaussian
                  0.13 0.55 0.13];   % green  - exponential

% EF -> N_D mapping for axis labels
EF_map = [-4.05, -4.104, -4.223, -4.342];
ND_map = [8e19,   1e19,   1e17,   1e15];

EF_vals_all  = unique(T_success.EF_peak_eV);
SRV_vals_all = unique(T_success.SRV_cm_s);

% Map EF values to N_D labels
ND_axis = NaN(size(EF_vals_all));
for ii = 1:numel(EF_vals_all)
    [~, idx_map] = min(abs(EF_map - EF_vals_all(ii)));
    ND_axis(ii) = ND_map(idx_map);
end

ND_labels = arrayfun(@(x) sprintf('%.0e', x), ND_axis, 'UniformOutput', false);
SRV_labels = arrayfun(@(x) sprintf('%.0e', x), SRV_vals_all, 'UniformOutput', false);

%% ========================================================================
%  COLORBLIND-FRIENDLY PALETTE (used by uniform-profile figures 14-29)
%% ========================================================================
CB_BLUE   = [0.12 0.47 0.71];
CB_ORANGE = [0.85 0.37 0.01];

CB_RAMP5 = [ ...
    0.12 0.47 0.71;   % level 1 (best / lowest)
    0.50 0.70 0.86;   % level 2
    0.60 0.60 0.60;   % level 3 (neutral)
    0.93 0.60 0.25;   % level 4
    0.85 0.37 0.01];  % level 5 (worst / highest)

CB_CAT4 = [ ...
    0.50 0.50 0.50;   % baseline (grey)
    0.85 0.37 0.01;   % SF eta=0 (orange)
    0.55 0.34 0.64;   % SF eta=1 (purple)
    0.12 0.47 0.71];  % SF eta=2 (blue)

% Shared optimal params for uniform-profile figures
optimal_SRV    = 10000;
optimal_EF     = -4.104;    % ~1e19 doping
optimal_PhiOff = 0;
optimal_jd     = 200;       % shallowest junction

% Filter for uniform profile only (used by figures 14-29)
T_uni_success = T_success(strcmp(T_success.profile_type_col, 'uniform'), :);

% SRV levels (may be a subset of SRV_vals_all)
SRV_list = sort(unique(T_uni_success.SRV_cm_s))';
if isempty(SRV_list)
    SRV_list = SRV_vals_all';
end

% Junction depth levels
jd_vals_all = sort(unique(T_uni_success.junction_depth_nm))';
if isempty(jd_vals_all)
    jd_vals_all = [200, 250, 300, 350, 400];
end


%% ========================================================================
%  FIGURE 4: EQE Comparison - 3 Profiles Overlaid (best-case params)
%% ========================================================================
if ismember(4, figures_to_plot)
figure('Name', 'EQE_3Profiles', 'Position', [200 200 700 500], 'Color', 'w');
hold on;

% Use SF=1, eta=2, EF=-4.104, SRV=1000, PhiOff=0 as "best case"
EF_fig4  = -4.104;
SRV_fig4 = 1000;
eta_fig4 = 2;

h_prof = gobjects(n_prof, 1);
for j = 1:n_prof
    mask = (T_success.SF_flag == 1) & ...
           (abs(T_success.eta_SF - eta_fig4) < tol) & ...
           (abs(T_success.EF_peak_eV - EF_fig4) < tol) & ...
           (T_success.SRV_cm_s == SRV_fig4) & ...
           (abs(T_success.Phi_R_offset_eV) < tol) & ...
           strcmp(T_success.profile_type_col, profile_types{j});
    idx = find(mask, 1);
    if ~isempty(idx)
        run_idx = T_success.run_number(idx);
        r = runs(run_idx);
        h_prof(j) = plot(r.EQE.wavelengths, r.EQE.EQE, '-', ...
                         'LineWidth', 2.5, 'Color', profile_colors(j, :));
    end
end

plot([400 800], [100 100], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

xlabel('Wavelength (nm)', 'FontSize', 20);
ylabel('EQE (%)', 'FontSize', 20);
xlim([400 800]);
set(gca, 'FontSize', 16, 'LineWidth', 1.2, 'Box', 'on');
grid on;

valid_h = h_prof(isgraphics(h_prof));
valid_labels = profile_labels(isgraphics(h_prof));
legend(valid_h, valid_labels, 'Location', 'southeast', 'FontSize', 16);

title(sprintf('EQE: SF \\eta=%.0f, N_0\\approx1e19, s_p=%.0e, \\Phi_{off}=0', ...
      eta_fig4, SRV_fig4), 'FontSize', 16);
hold off;

end

%% ========================================================================
%  FIGURE 5: Multi-panel EQE - Profile Effect at Different SRV
%% ========================================================================
if ismember(5, figures_to_plot)
SRV_panels = [1000, 10000, 100000];
n_panels = numel(SRV_panels);

figure('Name', 'EQE_Profile_vs_SRV', ...
       'Position', [50 50, 400*n_panels, 400], 'Color', 'w');

EF_fig5 = -4.104;
eta_fig5 = 2;

for sp = 1:n_panels
    subplot(1, n_panels, sp);
    hold on;

    for j = 1:n_prof
        mask = (T_success.SF_flag == 1) & ...
               (abs(T_success.eta_SF - eta_fig5) < tol) & ...
               (abs(T_success.EF_peak_eV - EF_fig5) < tol) & ...
               (T_success.SRV_cm_s == SRV_panels(sp)) & ...
               (abs(T_success.Phi_R_offset_eV) < tol) & ...
               strcmp(T_success.profile_type_col, profile_types{j});
        idx = find(mask, 1);
        if ~isempty(idx)
            r = runs(T_success.run_number(idx));
            plot(r.EQE.wavelengths, r.EQE.EQE, '-', ...
                 'LineWidth', 2.5, 'Color', profile_colors(j, :));
        end
    end

    plot([400 800], [100 100], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlabel('Wavelength (nm)', 'FontSize', 18);
    if sp == 1
        ylabel('EQE (%)', 'FontSize', 18);
    end
    title(sprintf('s_p = %.0e cm/s', SRV_panels(sp)), 'FontSize', 16);
    xlim([400 800]);
    ylim([40 130]);
    set(gca, 'FontSize', 14, 'Box', 'on');
    grid on;

    if sp == n_panels
        legend(profile_labels, 'Location', 'southeast', 'FontSize', 12);
    end
    hold off;
end

end

%% ========================================================================
%  FIGURE 6: Line Plot - PCE vs SRV for Each Profile
%% ========================================================================
if ismember(6, figures_to_plot)
figure('Name', 'PCE_vs_SRV_byProfile', 'Position', [250 250 700 500], 'Color', 'w');
hold on;

EF_fig6 = -4.104;
line_styles = {'-o', '-s', '-^'};

% --- SF=0 ---
for j = 1:n_prof
    PCE_vs_SRV = NaN(numel(SRV_vals_all), 1);
    for s = 1:numel(SRV_vals_all)
        mask = (T_success.SF_flag == 0) & ...
               (abs(T_success.EF_peak_eV - EF_fig6) < tol) & ...
               (T_success.SRV_cm_s == SRV_vals_all(s)) & ...
               (abs(T_success.Phi_R_offset_eV) < tol) & ...
               strcmp(T_success.profile_type_col, profile_types{j});
        if any(mask)
            PCE_vs_SRV(s) = mean(T_success.PCE_pct(mask), 'omitnan');
        end
    end
    semilogx(SRV_vals_all, PCE_vs_SRV, line_styles{j}, ...
             'LineWidth', 2, 'MarkerSize', 8, ...
             'Color', profile_colors(j, :), ...
             'MarkerFaceColor', profile_colors(j, :), ...
             'DisplayName', sprintf('%s (no SF)', profile_labels{j}));
end

% --- SF=1, eta=2 (dashed) ---
for j = 1:n_prof
    PCE_vs_SRV = NaN(numel(SRV_vals_all), 1);
    for s = 1:numel(SRV_vals_all)
        mask = (T_success.SF_flag == 1) & ...
               (abs(T_success.eta_SF - 2) < tol) & ...
               (abs(T_success.EF_peak_eV - EF_fig6) < tol) & ...
               (T_success.SRV_cm_s == SRV_vals_all(s)) & ...
               (abs(T_success.Phi_R_offset_eV) < tol) & ...
               strcmp(T_success.profile_type_col, profile_types{j});
        if any(mask)
            PCE_vs_SRV(s) = mean(T_success.PCE_pct(mask), 'omitnan');
        end
    end
    semilogx(SRV_vals_all, PCE_vs_SRV, ['--' line_styles{j}(2)], ...
             'LineWidth', 2, 'MarkerSize', 8, ...
             'Color', profile_colors(j, :), ...
             'MarkerFaceColor', 'none', ...
             'DisplayName', sprintf('%s (SF \\eta=2)', profile_labels{j}));
end

xlabel('Surface recombination velocity s_p (cm/s)', 'FontSize', 20);
ylabel('PCE (%)', 'FontSize', 20);
set(gca, 'XScale', 'log', 'FontSize', 16, 'LineWidth', 1.2, 'Box', 'on');
grid on;
legend('Location', 'southwest', 'FontSize', 11, 'NumColumns', 2);
title(sprintf('PCE vs SRV by Profile (N_0 \\approx 1e19)'), 'FontSize', 16);
hold off;

end

%% ========================================================================
%  FIGURE 7: Parity Plot - PCE(exponential) vs PCE(uniform)
%% ========================================================================
if ismember(7, figures_to_plot)
figure('Name', 'Parity_Exp_vs_Uni', 'Position', [300 300 550 500], 'Color', 'w');
hold on;

% Match runs with same (SF, eta, EF, SRV, PhiOff) but different profiles
T_uni = T_success(strcmp(T_success.profile_type_col, 'uniform'), :);
T_exp = T_success(strcmp(T_success.profile_type_col, 'exponential'), :);
T_gau = T_success(strcmp(T_success.profile_type_col, 'gaussian'), :);

% Match exponential to uniform
[PCE_uni_e, PCE_exp] = match_profiles(T_uni, T_exp, tol);
% Match gaussian to uniform
[PCE_uni_g, PCE_gau] = match_profiles(T_uni, T_gau, tol);

% Plot
scatter(PCE_uni_e, PCE_exp, 50, profile_colors(3,:), 'filled', ...
        'MarkerFaceAlpha', 0.7, 'DisplayName', 'Exponential');
scatter(PCE_uni_g, PCE_gau, 50, profile_colors(2,:), 'filled', ...
        'MarkerFaceAlpha', 0.7, 'DisplayName', 'Gaussian');

% Parity line
all_pce = [PCE_uni_e; PCE_uni_g; PCE_exp; PCE_gau];
lims = [min(all_pce)*0.95, max(all_pce)*1.02];
plot(lims, lims, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

xlabel('PCE Uniform (%)', 'FontSize', 20);
ylabel('PCE Graded (%)', 'FontSize', 20);
xlim(lims);
ylim(lims);
set(gca, 'FontSize', 16, 'LineWidth', 1.2, 'Box', 'on');
axis square;
grid on;
legend('Location', 'northwest', 'FontSize', 14);
title('Parity: Graded vs Uniform Doping', 'FontSize', 16);
hold off;

end

%% ========================================================================
%  FIGURE 10: Box Plot - PCE Distribution by Profile Type
%% ========================================================================
if ismember(10, figures_to_plot)
figure('Name', 'PCE_BoxPlot_Profile', 'Position', [450 150 500 450], 'Color', 'w');

PCE_by_prof = cell(n_prof, 1);
for j = 1:n_prof
    mask = strcmp(T_success.profile_type_col, profile_types{j});
    PCE_by_prof{j} = T_success.PCE_pct(mask);
end

all_PCE    = [];
all_groups = [];
for j = 1:n_prof
    all_PCE    = [all_PCE;    PCE_by_prof{j}]; %#ok<AGROW>
    all_groups = [all_groups; j * ones(numel(PCE_by_prof{j}), 1)]; %#ok<AGROW>
end

boxplot(all_PCE, all_groups, 'Colors', 'k', 'Widths', 0.6);
hold on;

set(gca, 'XTick', 1:n_prof, 'XTickLabel', profile_labels);

for j = 1:n_prof
    x_jitter = j + 0.2 * (rand(numel(PCE_by_prof{j}), 1) - 0.5);
    scatter(x_jitter, PCE_by_prof{j}, 25, profile_colors(j, :), ...
            'filled', 'MarkerFaceAlpha', 0.4);
end

ylabel('PCE (%)', 'FontSize', 20);
set(gca, 'FontSize', 18, 'LineWidth', 1.2, 'Box', 'on');
grid on;

for j = 1:n_prof
    med_val = median(PCE_by_prof{j}, 'omitnan');
    n_runs  = numel(PCE_by_prof{j});
    yl = ylim;
    text(j, yl(2) - 0.05*(yl(2)-yl(1)), ...
         sprintf('%.2f%%\n(n=%d)', med_val, n_runs), ...
         'HorizontalAlignment', 'center', 'FontSize', 14);
end

title('PCE Distribution by Doping Profile', 'FontSize', 16);
hold off;

end

%% ========================================================================
%  FIGURE 12: Multi-panel EQE - Profile Effect at Different N_D
%% ========================================================================
if ismember(12, figures_to_plot)
n_EF_panels = numel(EF_vals_all) - 1;  % skip first ND
figure('Name', 'EQE_Profile_vs_ND', ...
       'Position', [50 50, 350*n_EF_panels, 400], 'Color', 'w');

SRV_fig12 = 10000;
eta_fig12 = 2;

for ie = 2:numel(EF_vals_all)
    subplot(1, n_EF_panels, ie - 1);  % panel index starts at 1
    hold on;

    for j = 1:n_prof
        mask = (T_success.SF_flag == 1) & ...
               (abs(T_success.eta_SF - eta_fig12) < tol) & ...
               (abs(T_success.EF_peak_eV - EF_vals_all(ie)) < tol) & ...
               (T_success.SRV_cm_s == SRV_fig12) & ...
               (abs(T_success.Phi_R_offset_eV) < tol) & ...
               strcmp(T_success.profile_type_col, profile_types{j});
        idx = find(mask, 1);
        if ~isempty(idx)
            r = runs(T_success.run_number(idx));
            plot(r.EQE.wavelengths, r.EQE.EQE, '-', ...
                 'LineWidth', 2.5, 'Color', profile_colors(j, :));
        end
    end

    plot([400 800], [100 100], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlabel('Wavelength (nm)', 'FontSize', 16);
    if ie == 2
        ylabel('EQE (%)', 'FontSize', 16);
    end
    title(sprintf('N_0 \\approx %s', ND_labels{ie}), 'FontSize', 14);
    xlim([400 800]);
    ylim([40 130]);
    set(gca, 'FontSize', 13, 'Box', 'on');
    grid on;

    if ie == numel(EF_vals_all)
        legend(profile_labels, 'Location', 'southeast', 'FontSize', 11);
    end
    hold off;
end

sgtitle(sprintf('EQE by Profile at Different Doping (SF \\eta=2, s_p=%.0e)', ...
        SRV_fig12), 'FontSize', 16);

end

%% ========================================================================
%  FIGURE 13: Multi-panel EQE - Qf Effect by Profile Type
%
%  Layout: 2 rows x 3 columns
%    Rows:    SRV = 1000 cm/s (top), SRV = 10000 cm/s (bottom)
%    Columns: Uniform | Gaussian | Exponential
%  Colors match Figure 6 Panel C convention:
%    Qf > 0  -> green  [0   0.5 0  ]
%    Qf = 0  -> grey   [0.5 0.5 0.5]
%    Qf < 0  -> red    [0.8 0   0  ]
%  SF: on only (eta_SF = 2)
%  Fixed: EF = -4.104 (N0 ~ 1e19 cm^-3)
%% ========================================================================
if ismember(13, figures_to_plot)

EF_fig13  = -4.104;
SRV_rows  = [1000];
eta_fig13 = 2;

% Collect and sort Phi_R_offset values (ascending: most negative first)
PhiOff_vals = sort(unique(T_success.Phi_R_offset_eV));
n_Qf = numel(PhiOff_vals);

% Assign colors matching Fig 6 Panel C convention
Qf_colors = zeros(n_Qf, 3);
Qf_lbls   = cell(n_Qf, 1);
for iq = 1:n_Qf
    phi = PhiOff_vals(iq);
    if phi < -tol
        Qf_colors(iq, :) = [0.8, 0.0, 0.0];   % red  - Qf < 0
        Qf_lbls{iq}      = 'Q_f < 0';
    elseif phi > tol
        Qf_colors(iq, :) = [0.0, 0.5, 0.0];   % green - Qf > 0
        Qf_lbls{iq}      = 'Q_f > 0';
    else
        Qf_colors(iq, :) = [0.5, 0.5, 0.5];   % grey  - Qf = 0
        Qf_lbls{iq}      = 'Q_f = 0';
    end
end

n_rows = numel(SRV_rows);
figure('Name', 'EQE_Qf_byProfile_2SRV', ...
       'Position', [50 50, 430*n_prof, 370*n_rows], 'Color', 'w');

for row = 1:n_rows
    SRV_val = SRV_rows(row);

    for col = 1:n_prof
        sp_idx = (row - 1) * n_prof + col;
        subplot(n_rows, n_prof, sp_idx);
        hold on;

        h_Qf = gobjects(n_Qf, 1);

        for iq = 1:n_Qf
            mask = (T_success.SF_flag == 1) & ...
                   (abs(T_success.eta_SF - eta_fig13) < tol) & ...
                   (abs(T_success.EF_peak_eV - EF_fig13) < tol) & ...
                   (T_success.SRV_cm_s == SRV_val) & ...
                   (abs(T_success.Phi_R_offset_eV - PhiOff_vals(iq)) < tol) & ...
                   strcmp(T_success.profile_type_col, profile_types{col});
            idx = find(mask, 1);
            if ~isempty(idx)
                r = runs(T_success.run_number(idx));
                h_Qf(iq) = plot(r.EQE.wavelengths, r.EQE.EQE, '-', ...
                                'LineWidth', 2.5, 'Color', Qf_colors(iq, :));
            end
        end

        plot([400 800], [100 100], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

        if col == 1
            ylabel('EQE (%)', 'FontSize', 15);
        end
        if row == n_rows
            xlabel('Wavelength (nm)', 'FontSize', 15);
        end
        if row == 1
            title(profile_labels{col}, 'FontSize', 15);
        end

        xlim([400 800]);
        ylim([40 130]);
        set(gca, 'FontSize', 13, 'LineWidth', 1.2, 'Box', 'on');
        grid on;

        % Legend on bottom-right panel only
        if row == n_rows && col == n_prof
            valid = isgraphics(h_Qf);
            if any(valid)
                legend(h_Qf(valid), Qf_lbls(valid), ...
                       'Location', 'southeast', 'FontSize', 12);
            end
        end

        hold off;
    end
end

sgtitle(sprintf('EQE vs Q_f by Doping Profile  (\\eta_{SF}=%.0f, N_0\\approx1e19, s_p=%.0e cm/s)', ...
        eta_fig13, SRV_rows(1)), 'FontSize', 15);

end

%% ########################################################################
%  UNIFORM-PROFILE FIGURES (adapted from main/main_baseline)
%  All figures below filter to profile_type = 'uniform' and use
%  EF_peak_eV instead of EF_emitter_eV. No junction_depth_nm.
%% ########################################################################

%% ========================================================================
%  FIGURE 14: EQE Comparison - SF Signature (EQE > 100%)
%% ========================================================================
if ismember(14, figures_to_plot)
figure('Name', 'EQE_SF_Comparison', 'Position', [50 50 900 600], 'Color', 'w');

for eta_val = [0, 2]
    mask = (T_uni_success.SF_flag == 1) & ...
           (abs(T_uni_success.eta_SF - eta_val) < tol) & ...
           (T_uni_success.SRV_cm_s == optimal_SRV) & ...
           (abs(T_uni_success.EF_peak_eV - optimal_EF) < tol) & ...
           (abs(T_uni_success.Phi_R_offset_eV - optimal_PhiOff) < tol) & ...
           (T_uni_success.junction_depth_nm == optimal_jd);

    idx = find(mask);
    if ~isempty(idx)
        run_idx = T_uni_success.run_number(idx(1));
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

% Baseline (SF=0)
mask = (T_uni_success.SF_flag == 0) & ...
       (T_uni_success.SRV_cm_s == optimal_SRV) & ...
       (abs(T_uni_success.EF_peak_eV - optimal_EF) < tol) & ...
       (abs(T_uni_success.Phi_R_offset_eV - optimal_PhiOff) < tol) & ...
       (T_uni_success.junction_depth_nm == optimal_jd);
idx = find(mask);
if ~isempty(idx)
    run_idx = T_uni_success.run_number(idx(1));
    r = runs(run_idx);
    plot(r.EQE.wavelengths, r.EQE.EQE, '-', 'LineWidth', 2.5, ...
         'Color', [0.50 0.50 0.50], ...
         'DisplayName', 'No SF layer (baseline)');
end

xlims = [400 650];
plot(xlims, [100 100], 'k--', 'LineWidth', 1.5, 'DisplayName', 'EQE = 100%');
fill([xlims(1) xlims(2) xlims(2) xlims(1)], [100 100 150 150], ...
     [0.85 0.92 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.3, ...
     'HandleVisibility', 'off');

xlabel('Wavelength (nm)', 'FontSize', 18);
ylabel('EQE (%)', 'FontSize', 18);
title('Singlet Fission Enhancement of EQE (Uniform Profile)', 'FontSize', 20);
xlim(xlims);
ylim([0 150]);
set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'on');
grid on;
legend('Location', 'northeast', 'FontSize', 12);
text(420, 130, 'SF benefit region', 'FontSize', 12, 'Color', CB_BLUE);
hold off;

end

%% ========================================================================
%  FIGURE 15: Heatmap - PCE vs (SRV, N_D) for ideal SF
%% ========================================================================
if ismember(15, figures_to_plot)
figure('Name', 'Heatmap_PCE_SRV_ND', 'Position', [100 100 600 450], 'Color', 'w');

mask = (T_uni_success.SF_flag == 1) & ...
       (abs(T_uni_success.eta_SF - 2) < tol) & ...
       (T_uni_success.junction_depth_nm == 300) & ...
       (abs(T_uni_success.Phi_R_offset_eV - 0) < tol);

T_sub = T_uni_success(mask, :);

SRV_vals = unique(T_sub.SRV_cm_s);
EF_vals  = unique(T_sub.EF_peak_eV);

ND_axis_15 = NaN(size(EF_vals));
for ii = 1:numel(EF_vals)
    [~, idx_map] = min(abs(EF_map - EF_vals(ii)));
    ND_axis_15(ii) = ND_map(idx_map);
end

PCE_matrix = NaN(length(EF_vals), length(SRV_vals));
for i = 1:length(EF_vals)
    for j = 1:length(SRV_vals)
        idx = find(abs(T_sub.EF_peak_eV - EF_vals(i)) < tol & ...
                   abs(T_sub.SRV_cm_s - SRV_vals(j))/SRV_vals(j) < tol, 1);
        if ~isempty(idx)
            PCE_matrix(i, j) = T_sub.PCE_pct(idx);
        end
    end
end

h = imagesc(1:length(SRV_vals), 1:length(EF_vals), PCE_matrix);
colormap(blueorange(256));
set(h, 'AlphaData', ~isnan(PCE_matrix));
set(gca, 'Color', 'w');

cb = colorbar;
ylabel(cb, 'PCE (%)', 'FontSize', 18);
set(gca, 'YDir', 'normal');

set(gca, 'XTick', 1:length(SRV_vals), ...
         'XTickLabel', arrayfun(@(x) sprintf('%.0e', x), SRV_vals, 'UniformOutput', false));
set(gca, 'YTick', 1:length(EF_vals), ...
         'YTickLabel', arrayfun(@(x) sprintf('%.1e', x), ND_axis_15, 'UniformOutput', false));

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
%  FIGURE 16: Line Plot - SF Benefit vs SRV
%% ========================================================================
if ismember(16, figures_to_plot)
figure('Name', 'SF_Benefit_vs_SRV', 'Position', [150 150 450 500], 'Color', 'w');

nSRV = numel(SRV_list);
PCE_baseline_16 = NaN(nSRV, 1);
PCE_SF_16       = NaN(nSRV, 1);

for k = 1:nSRV
    mask = (T_uni_success.SF_flag == 0) & ...
           (T_uni_success.SRV_cm_s == SRV_list(k)) & ...
           (abs(T_uni_success.EF_peak_eV - optimal_EF) < tol) & ...
           (T_uni_success.junction_depth_nm == 300) & ...
           (abs(T_uni_success.Phi_R_offset_eV - 0) < tol);
    if any(mask)
        PCE_baseline_16(k) = mean(T_uni_success.PCE_pct(mask), 'omitnan');
    end

    mask = (T_uni_success.SF_flag == 1) & ...
           (abs(T_uni_success.eta_SF - 2) < tol) & ...
           (T_uni_success.SRV_cm_s == SRV_list(k)) & ...
           (abs(T_uni_success.EF_peak_eV - optimal_EF) < tol) & ...
           (T_uni_success.junction_depth_nm == 300) & ...
           (abs(T_uni_success.Phi_R_offset_eV - 0) < tol);
    if any(mask)
        PCE_SF_16(k) = mean(T_uni_success.PCE_pct(mask), 'omitnan');
    end
end

hold on;

h1 = semilogx(SRV_list, PCE_baseline_16, '-o', ...
              'LineWidth', 2.5, 'MarkerSize', 10, ...
              'MarkerFaceColor', [0.50 0.50 0.50], ...
              'Color', [0.50 0.50 0.50], ...
              'DisplayName', 'Baseline');

h2 = semilogx(SRV_list, PCE_SF_16, '-s', ...
              'LineWidth', 2.5, 'MarkerSize', 10, ...
              'MarkerFaceColor', CB_BLUE, ...
              'Color', CB_BLUE, ...
              'DisplayName', 'SF (\eta=2)');

ylabel('PCE (%)', 'FontSize', 20);
xlabel('Surface recombination velocity s_p (cm/s)', 'FontSize', 20);

legend([h1 h2], {'Baseline', 'SF (\eta=2)'}, ...
       'Location', 'northeast', 'FontSize', 16);

set(gca, 'XScale', 'log', 'XTick', SRV_list);
set(gca, 'TickLabelInterpreter', 'tex');
xticklabels(arrayfun(@(x) sprintf('10^{%d}', round(log10(x))), ...
                     SRV_list, 'UniformOutput', false));
ax = gca;
ax.XMinorTick = 'off';
ax.XMinorGrid = 'off';
grid on;
xlim([min(SRV_list)*0.8, max(SRV_list)*1.2]);
set(gca, 'Box', 'on', 'LineWidth', 1.5, 'FontSize', 18);

ymin16 = nanmin([PCE_baseline_16(:); PCE_SF_16(:)]);
ymax16 = nanmax([PCE_baseline_16(:); PCE_SF_16(:)]);
dy16   = 0.02*(ymax16 - ymin16);
for k = 1:nSRV
    if ~isnan(PCE_baseline_16(k))
        text(SRV_list(k), PCE_baseline_16(k) + dy16, ...
             sprintf('%.2f%%', PCE_baseline_16(k)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontSize', 16, 'Color', [0.30 0.30 0.30]);
    end
    if ~isnan(PCE_SF_16(k))
        text(SRV_list(k), PCE_SF_16(k) + 2*dy16, ...
             sprintf('%.2f%%', PCE_SF_16(k)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontSize', 16, 'Color', CB_BLUE);
    end
end

hold off;
end

%% ========================================================================
%  FIGURE 17: Integrated EQE in Blue Region (405-500nm) vs eta_SF
%% ========================================================================
if ismember(17, figures_to_plot)
figure('Name', 'Blue_EQE_Line', 'Position', [250 250 700 500], 'Color', 'w');

eta_vals_17 = [0, 0.5, 1, 1.5, 2];
n_eta_17    = numel(eta_vals_17);
blue_EQE_17 = NaN(1 + n_eta_17, 1);

% Get wavelength vector from a successful uniform run
r_sample_17 = [];
for rr = 1:numel(runs)
    if strcmp(runs(rr).status, 'success')
        r_sample_17 = runs(rr);
        break;
    end
end

if ~isempty(r_sample_17)
    wl17 = r_sample_17.EQE.wavelengths;
    blue_mask_wl = (wl17 >= 405) & (wl17 <= 500);

    % Baseline
    mask = (T_uni_success.SF_flag == 0) & ...
           (T_uni_success.SRV_cm_s == optimal_SRV) & ...
           (abs(T_uni_success.EF_peak_eV - optimal_EF) < tol) & ...
           (T_uni_success.junction_depth_nm == optimal_jd) & ...
           (abs(T_uni_success.Phi_R_offset_eV - 0) < tol);
    if any(mask)
        run_idx = T_uni_success.run_number(find(mask, 1));
        blue_EQE_17(1) = mean(runs(run_idx).EQE.EQE(blue_mask_wl), 'omitnan');
    end

    % SF runs at each eta
    for e = 1:n_eta_17
        mask = (T_uni_success.SF_flag == 1) & ...
               (abs(T_uni_success.eta_SF - eta_vals_17(e)) < tol) & ...
               (T_uni_success.SRV_cm_s == optimal_SRV) & ...
               (abs(T_uni_success.EF_peak_eV - optimal_EF) < tol) & ...
               (T_uni_success.junction_depth_nm == optimal_jd) & ...
               (abs(T_uni_success.Phi_R_offset_eV - 0) < tol);
        if any(mask)
            run_idx = T_uni_success.run_number(find(mask, 1));
            blue_EQE_17(e+1) = mean(runs(run_idx).EQE.EQE(blue_mask_wl), 'omitnan');
        end
    end

    hold on;
    if ~isnan(blue_EQE_17(1))
        h_base17 = yline(blue_EQE_17(1), '--', 'Baseline', ...
                         'LineWidth', 2, 'Color', [0.50 0.50 0.50]);
    else
        h_base17 = [];
    end

    sf_vals17 = blue_EQE_17(2:end);
    valid_sf17 = ~isnan(sf_vals17);
    h_sf17 = plot(eta_vals_17(valid_sf17), sf_vals17(valid_sf17), '-o', ...
                  'LineWidth', 2.5, 'MarkerSize', 8, ...
                  'Color', CB_BLUE, 'MarkerFaceColor', [0.50 0.70 0.86], ...
                  'DisplayName', 'SF device');

    h_100_17 = yline(100, 'k--', 'EQE = 100%', 'LineWidth', 1.5);

    xlabel('\eta_{SF}', 'FontSize', 16);
    ylabel('Mean EQE in Blue Region (405-500nm) (%)', 'FontSize', 16);
    title('Blue-Region EQE Enhancement vs \eta_{SF} (Uniform)', 'FontSize', 18);
    xlim([min(eta_vals_17) - 0.1, max(eta_vals_17) + 0.1]);

    valid_vals17 = blue_EQE_17(~isnan(blue_EQE_17));
    if ~isempty(valid_vals17)
        ylim([0, max(valid_vals17)*1.15]);
    end
    set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'on');
    grid on;

    if ~isempty(h_base17)
        legend([h_sf17, h_base17, h_100_17], ...
               {'SF device', 'Baseline', 'EQE = 100%'}, ...
               'Location', 'best', 'FontSize', 12);
    else
        legend([h_sf17, h_100_17], {'SF device', 'EQE = 100%'}, ...
               'Location', 'best', 'FontSize', 12);
    end

    for k = 1:n_eta_17
        if ~isnan(sf_vals17(k))
            text(eta_vals_17(k), sf_vals17(k) + 3, ...
                 sprintf('%.1f%%', sf_vals17(k)), ...
                 'HorizontalAlignment', 'center', ...
                 'FontSize', 12, 'FontWeight', 'bold');
        end
    end
    hold off;
end

end

%% ========================================================================
%  FIGURE 18: Multi-panel EQE - effect of eta, SRV, N_D (with baselines)
%  Layout: 2x2 but only 3 panels (skip junction depth)
%% ========================================================================
if ismember(18, figures_to_plot)
figure('Name', 'EQE_Parameter_Effects', 'Position', [50 50 1000 800], 'Color', 'w');

%% Panel A: Effect of eta_SF
subplot(2, 2, 1);
hold on;

eta_vals_18 = [2, 1.5, 1, 0.5, 0];
eta_colors_18 = [ ...
    0.12 0.47 0.71;
    0.50 0.70 0.86;
    0.60 0.60 0.60;
    0.93 0.60 0.25;
    0.85 0.37 0.01];

h_eta = gobjects(numel(eta_vals_18), 1);
for i = 1:numel(eta_vals_18)
    mask = (T_uni_success.SF_flag == 1) & ...
           (abs(T_uni_success.eta_SF - eta_vals_18(i)) < tol) & ...
           (T_uni_success.SRV_cm_s == optimal_SRV) & ...
           (abs(T_uni_success.EF_peak_eV - optimal_EF) < tol) & ...
           (T_uni_success.junction_depth_nm == optimal_jd) & ...
           (abs(T_uni_success.Phi_R_offset_eV - 0) < tol);
    if any(mask)
        run_idx = T_uni_success.run_number(find(mask, 1));
        r = runs(run_idx);
        h_eta(i) = plot(r.EQE.wavelengths, r.EQE.EQE, '-', ...
                        'LineWidth', 2, 'Color', eta_colors_18(i, :));
    end
end

% Baseline
mask_base_A = (T_uni_success.SF_flag == 0) & ...
              (T_uni_success.SRV_cm_s == optimal_SRV) & ...
              (abs(T_uni_success.EF_peak_eV - optimal_EF) < tol) & ...
              (T_uni_success.junction_depth_nm == optimal_jd) & ...
              (abs(T_uni_success.Phi_R_offset_eV - 0) < tol);
h_base_A = gobjects(0);
if any(mask_base_A)
    run_idx = T_uni_success.run_number(find(mask_base_A, 1));
    r = runs(run_idx);
    h_base_A = plot(r.EQE.wavelengths, r.EQE.EQE, 'k--', 'LineWidth', 1.5);
end

plot([400 800], [100 100], 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Wavelength (nm)', 'FontSize', 20);
ylabel('EQE (%)', 'FontSize', 20);
xlim([400 800]);
% ylim([50 120]);
set(gca, 'FontSize', 13, 'Box', 'on'); grid on;

leg_labels_eta = arrayfun(@(v) sprintf('\\eta_{SF} = %s', num2str(v)), ...
                          eta_vals_18, 'UniformOutput', false);
valid_eta = isgraphics(h_eta);
if isgraphics(h_base_A)
    legend([h_eta(valid_eta); h_base_A], [leg_labels_eta(valid_eta), {'Baseline'}], ...
           'Location', 'southeast', 'FontSize', 13);
elseif any(valid_eta)
    legend(h_eta(valid_eta), leg_labels_eta(valid_eta), 'Location', 'southeast', 'FontSize', 11);
end
hold off;

%% Panel B: Effect of SRV (with baseline dashed)
subplot(2, 2, 2);
hold on;

SRV_colors_18 = CB_RAMP5;
% If SRV_list has different number of entries, pick colors accordingly
n_srv_18 = numel(SRV_list);
if n_srv_18 <= 5
    srv_cmap = CB_RAMP5(round(linspace(1, 5, n_srv_18)), :);
else
    srv_cmap = interp1(linspace(0,1,5), CB_RAMP5, linspace(0,1,n_srv_18));
end

h_srv = gobjects(n_srv_18, 1);
for s = 1:n_srv_18
    mask = (T_uni_success.SF_flag == 1) & ...
           (abs(T_uni_success.eta_SF - 2) < tol) & ...
           (T_uni_success.SRV_cm_s == SRV_list(s)) & ...
           (abs(T_uni_success.EF_peak_eV - optimal_EF) < tol) & ...
           (T_uni_success.junction_depth_nm == optimal_jd) & ...
           (abs(T_uni_success.Phi_R_offset_eV - 0) < tol);
    if any(mask)
        run_idx = T_uni_success.run_number(find(mask, 1));
        r = runs(run_idx);
        h_srv(s) = plot(r.EQE.wavelengths, r.EQE.EQE, '-', 'LineWidth', 2, ...
                        'Color', srv_cmap(s, :), ...
                        'DisplayName', sprintf('s_p = %.0d', SRV_list(s)));
    end
end

% Baselines (dashed, same color)
for s = 1:n_srv_18
    mask_base = (T_uni_success.SF_flag == 0) & ...
                (T_uni_success.SRV_cm_s == SRV_list(s)) & ...
                (abs(T_uni_success.EF_peak_eV - optimal_EF) < tol) & ...
                (T_uni_success.junction_depth_nm == optimal_jd) & ...
                (abs(T_uni_success.Phi_R_offset_eV - 0) < tol);
    if any(mask_base)
        run_idx = T_uni_success.run_number(find(mask_base, 1));
        r = runs(run_idx);
        plot(r.EQE.wavelengths, r.EQE.EQE, '--', 'LineWidth', 1.5, ...
             'Color', srv_cmap(s, :), 'HandleVisibility', 'off');
    end
end

h_leg_solid  = plot(NaN, NaN, 'k-',  'LineWidth', 2);
h_leg_dashed = plot(NaN, NaN, 'k--', 'LineWidth', 1.5);

plot([400 800], [100 100], 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Wavelength (nm)', 'FontSize', 20);
ylabel('EQE (%)', 'FontSize', 20);
xlim([400 800]); 
% ylim([50 120]);
legend([h_srv(isgraphics(h_srv)); h_leg_solid; h_leg_dashed], ...
       [arrayfun(@(x) sprintf('s_p = %.0d', x), SRV_list(isgraphics(h_srv)), 'UniformOutput', false), ...
        {'SF (\eta=2)', 'Baseline'}], ...
       'Location', 'southeast', 'FontSize', 11);
set(gca, 'FontSize', 13, 'Box', 'on'); grid on;
hold off;

%% Panel C: Effect of junction depth j_d (with baseline dashed)
subplot(2, 2, 3);
hold on;

jd_vals_panel = jd_vals_all;
n_jd_18 = numel(jd_vals_panel);
if n_jd_18 <= 5
    jd_cmap = CB_RAMP5(round(linspace(1, 5, n_jd_18)), :);
else
    jd_cmap = interp1(linspace(0,1,5), CB_RAMP5, linspace(0,1,n_jd_18));
end

h_jd = gobjects(n_jd_18, 1);
for jj = 1:n_jd_18
    mask = (T_uni_success.SF_flag == 1) & ...
           (abs(T_uni_success.eta_SF - 2) < tol) & ...
           (T_uni_success.SRV_cm_s == 10000) & ...
           (abs(T_uni_success.EF_peak_eV - optimal_EF) < tol) & ...
           (T_uni_success.junction_depth_nm == jd_vals_panel(jj)) & ...
           (abs(T_uni_success.Phi_R_offset_eV - 0) < tol);
    if any(mask)
        run_idx = T_uni_success.run_number(find(mask, 1));
        r = runs(run_idx);
        h_jd(jj) = plot(r.EQE.wavelengths, r.EQE.EQE, '-', ...
                        'LineWidth', 2, 'Color', jd_cmap(jj, :));
    end
end

% Baselines (dashed, same color)
for jj = 1:n_jd_18
    mask_base = (T_uni_success.SF_flag == 0) & ...
                (T_uni_success.SRV_cm_s == 10000) & ...
                (abs(T_uni_success.EF_peak_eV - optimal_EF) < tol) & ...
                (T_uni_success.junction_depth_nm == jd_vals_panel(jj)) & ...
                (abs(T_uni_success.Phi_R_offset_eV - 0) < tol);
    if any(mask_base)
        run_idx = T_uni_success.run_number(find(mask_base, 1));
        r = runs(run_idx);
        plot(r.EQE.wavelengths, r.EQE.EQE, '--', 'LineWidth', 1.5, ...
             'Color', jd_cmap(jj, :), 'HandleVisibility', 'off');
    end
end

h_leg_solid_C  = plot(NaN, NaN, 'k-',  'LineWidth', 2);
h_leg_dashed_C = plot(NaN, NaN, 'k--', 'LineWidth', 1.5);

plot([400 800], [100 100], 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Wavelength (nm)', 'FontSize', 20);
ylabel('EQE (%)', 'FontSize', 20);
xlim([400 800]); 
% ylim([50 120]);
set(gca, 'FontSize', 13, 'Box', 'on'); grid on;

leg_labels_jd = arrayfun(@(d) sprintf('j_d = %d nm', d), ...
                         jd_vals_panel, 'UniformOutput', false);
valid_jd = isgraphics(h_jd);
legend([h_jd(valid_jd); h_leg_solid_C; h_leg_dashed_C], ...
       [leg_labels_jd(valid_jd), {'SF (\eta=2)', 'Baseline'}], ...
       'Location', 'southeast', 'FontSize', 11);
hold off;

%% Panel D: Effect of doping concentration N_D (with baseline dashed)
subplot(2, 2, 4);
hold on;

ND_vals_panel = [8e19, 1e19, 1e17, 1e15];
EF_vals_panel = [-4.05, -4.10, -4.22, -4.34];
ND_colors_18  = [ ...
    0.85 0.37 0.01;
    0.93 0.60 0.25;
    0.50 0.70 0.86;
    0.12 0.47 0.71];

h_nd = gobjects(numel(ND_vals_panel), 1);
for nn = 1:numel(ND_vals_panel)
    mask = (T_uni_success.SF_flag == 1) & ...
           (abs(T_uni_success.eta_SF - 2) < tol) & ...
           (T_uni_success.SRV_cm_s == optimal_SRV) & ...
           (abs(T_uni_success.EF_peak_eV - EF_vals_panel(nn)) < tol) & ...
           (T_uni_success.junction_depth_nm == optimal_jd) & ...
           (abs(T_uni_success.Phi_R_offset_eV - 0) < tol);
    if any(mask)
        run_idx = T_uni_success.run_number(find(mask, 1));
        r = runs(run_idx);
        h_nd(nn) = plot(r.EQE.wavelengths, r.EQE.EQE, '-', ...
                        'LineWidth', 2, 'Color', ND_colors_18(nn, :));
    end
end

% Baselines (dashed, same color)
for nn = 1:numel(ND_vals_panel)
    mask_base = (T_uni_success.SF_flag == 0) & ...
                (T_uni_success.SRV_cm_s == optimal_SRV) & ...
                (abs(T_uni_success.EF_peak_eV - EF_vals_panel(nn)) < tol) & ...
                (T_uni_success.junction_depth_nm == optimal_jd) & ...
                (abs(T_uni_success.Phi_R_offset_eV - 0) < tol);
    if any(mask_base)
        run_idx = T_uni_success.run_number(find(mask_base, 1));
        r = runs(run_idx);
        plot(r.EQE.wavelengths, r.EQE.EQE, '--', 'LineWidth', 1.5, ...
             'Color', ND_colors_18(nn, :), 'HandleVisibility', 'off');
    end
end

h_leg_solid_D  = plot(NaN, NaN, 'k-',  'LineWidth', 2);
h_leg_dashed_D = plot(NaN, NaN, 'k--', 'LineWidth', 1.5);

plot([400 800], [100 100], 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Wavelength (nm)', 'FontSize', 20);
ylabel('EQE (%)', 'FontSize', 20);
xlim([400 800]); 
% ylim([50 120]);
set(gca, 'FontSize', 13, 'Box', 'on'); grid on;

leg_labels_nd = arrayfun(@(d) sprintf('N_D = %.0e cm^{-3}', d), ...
                         ND_vals_panel, 'UniformOutput', false);
valid_nd = isgraphics(h_nd);
legend([h_nd(valid_nd); h_leg_solid_D; h_leg_dashed_D], ...
       [leg_labels_nd(valid_nd), {'SF (\eta=2)', 'Baseline'}], ...
       'Location', 'southeast', 'FontSize', 11);
hold off;

end  % end of figure 18

%% ========================================================================
%  FIGURE 19: EQE vs Phi_R offset (with baseline)
%% ========================================================================
if ismember(19, figures_to_plot)
figure('Name', 'EQE_vs_PhiRoffset_sp1e4', ...
       'Position', [200 200 430 400], 'Color', 'w');
hold on;

PhiOff_vals_19   = [0.0539, 0, -0.0539];
PhiOff_colors_19 = [CB_BLUE; 0.50 0.50 0.50; CB_ORANGE];
PhiOff_labels_19 = {'Q_f>0', 'Q_f=0', 'Q_f<0'};

SRV_plot_19 = 1e4;
nPhi_19 = numel(PhiOff_vals_19);
h_phi = gobjects(nPhi_19, 1);

for p = 1:nPhi_19
    mask = (T_uni_success.SF_flag == 1) & ...
           (abs(T_uni_success.eta_SF - 2) < tol) & ...
           (T_uni_success.SRV_cm_s == SRV_plot_19) & ...
           (abs(T_uni_success.EF_peak_eV - optimal_EF) < tol) & ...
           (T_uni_success.junction_depth_nm == optimal_jd) & ...
           (abs(T_uni_success.Phi_R_offset_eV - PhiOff_vals_19(p)) < tol);
    if any(mask)
        run_idx = T_uni_success.run_number(find(mask, 1));
        r = runs(run_idx);
        h_phi(p) = plot(r.EQE.wavelengths, r.EQE.EQE, '-', ...
                        'LineWidth', 2, 'Color', PhiOff_colors_19(p, :), ...
                        'DisplayName', PhiOff_labels_19{p});
    end
end

% Baseline
mask_base_phi = (T_uni_success.SF_flag == 0) & ...
                (T_uni_success.SRV_cm_s == SRV_plot_19) & ...
                (abs(T_uni_success.EF_peak_eV - optimal_EF) < tol) & ...
                (T_uni_success.junction_depth_nm == optimal_jd) & ...
                (abs(T_uni_success.Phi_R_offset_eV - 0) < tol);
h_base_phi = gobjects(0);
if any(mask_base_phi)
    run_idx = T_uni_success.run_number(find(mask_base_phi, 1));
    r = runs(run_idx);
    h_base_phi = plot(r.EQE.wavelengths, r.EQE.EQE, 'k--', ...
                      'LineWidth', 1.5, 'DisplayName', 'Baseline');
end

plot([400 800], [100 100], 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Wavelength (nm)', 'FontSize', 20);
ylabel('EQE (%)', 'FontSize', 20);
xlim([400 800]); 
% ylim([40 120]);
set(gca, 'FontSize', 15, 'Box', 'on'); grid on;

valid_phi = isgraphics(h_phi);
if isgraphics(h_base_phi)
    legend([h_phi(valid_phi); h_base_phi], [PhiOff_labels_19(valid_phi), {'Baseline'}], ...
           'Location', 'southeast', 'FontSize', 15);
elseif any(valid_phi)
    legend(h_phi(valid_phi), PhiOff_labels_19(valid_phi), 'Location', 'southeast', 'FontSize', 15);
end
hold off;

end

%% ========================================================================
%  FIGURE 20: Box Plot - PCE Distribution by eta_SF (uniform only)
%% ========================================================================
if ismember(20, figures_to_plot)
figure('Name', 'PCE_Distribution_eta', 'Position', [350 150 500 400], 'Color', 'w');

PCE_by_cat_20 = cell(4, 1);
mask = (T_uni_success.SF_flag == 0);
PCE_by_cat_20{1} = T_uni_success.PCE_pct(mask);
for e = 0:2
    mask = (T_uni_success.SF_flag == 1) & (abs(T_uni_success.eta_SF - e) < tol);
    PCE_by_cat_20{e+2} = T_uni_success.PCE_pct(mask);
end

all_PCE_20    = [];
all_groups_20 = [];
for k = 1:4
    all_PCE_20    = [all_PCE_20;    PCE_by_cat_20{k}]; %#ok<AGROW>
    all_groups_20 = [all_groups_20; k * ones(length(PCE_by_cat_20{k}), 1)]; %#ok<AGROW>
end

boxplot(all_PCE_20, all_groups_20, 'Colors', 'k', 'Widths', 0.6);
hold on;

set(gca, 'XTick', 1:4, ...
         'XTickLabel', {'Baseline', 'SF (\eta=0)', 'SF (\eta=1)', 'SF (\eta=2)'}, ...
         'TickLabelInterpreter', 'tex');

for k = 1:4
    x_jitter = k + 0.15 * (rand(length(PCE_by_cat_20{k}), 1) - 0.5);
    scatter(x_jitter, PCE_by_cat_20{k}, 30, CB_CAT4(k,:), ...
            'filled', 'MarkerFaceAlpha', 0.5);
end

ylabel('PCE (%)', 'FontSize', 20);
set(gca, 'FontSize', 18, 'LineWidth', 1.2, 'Box', 'on');
grid on;

med_vals_20 = zeros(1,4);
for k = 1:4
    med_vals_20(k) = median(PCE_by_cat_20{k}, 'omitnan');
end
yl = ylim;
y_top_vals = yl(2) - 0.08*(yl(2)-yl(1));
y_top_text = yl(2) - 0.01*(yl(2)-yl(1));
for k = 1:4
    if ~isnan(med_vals_20(k))
        text(k, y_top_vals, sprintf('%.2f%%', med_vals_20(k)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
             'FontSize', 16, 'Color', 'k', 'Clipping', 'off');
    end
end
text(2.5, y_top_text, 'Median PCE', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
     'FontSize', 14, 'Clipping', 'off');
title('PCE Distribution by \eta_{SF} (Uniform Profile)', 'FontSize', 16);
hold off;

end

%% ========================================================================
%  FIGURE 21: Correlation Matrix (uniform only)
%% ========================================================================
if ismember(21, figures_to_plot)
figure('Name', 'Correlation_Matrix', 'Position', [400 200 700 600], 'Color', 'w');

vars_21 = {'SF', '\eta_{SF}', 'j_d', 'log(SRV)', 'E_F', '\Phi_R off', 'PCE', 'V_{oc}', 'J_{sc}', 'FF'};
data_matrix_21 = [T_uni_success.SF_flag, T_uni_success.eta_SF, ...
                  T_uni_success.junction_depth_nm, ...
                  log10(T_uni_success.SRV_cm_s), T_uni_success.EF_peak_eV, ...
                  T_uni_success.Phi_R_offset_eV, ...
                  T_uni_success.PCE_pct, T_uni_success.Voc_V, ...
                  T_uni_success.Jsc_mA_cm2, T_uni_success.FF_pct];

R_21 = corrcoef(data_matrix_21, 'Rows', 'complete');
imagesc(R_21);
colormap(blueorange(256));
caxis([-1 1]);
cb = colorbar;
ylabel(cb, 'Correlation', 'FontSize', 12);

set(gca, 'XTick', 1:length(vars_21), 'XTickLabel', vars_21, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:length(vars_21), 'YTickLabel', vars_21);
title('Parameter-Metric Correlation (Uniform Profile)', 'FontSize', 16);
set(gca, 'FontSize', 11);

for i = 1:size(R_21, 1)
    for j = 1:size(R_21, 2)
        if abs(R_21(i,j)) > 0.3
            text(j, i, sprintf('%.2f', R_21(i,j)), 'HorizontalAlignment', 'center', ...
                 'FontSize', 9, 'Color', 'k', 'FontWeight', 'bold');
        end
    end
end

end

%% ========================================================================
%  FIGURE 22: Scatter Matrix - PCE vs Key Parameters (uniform only)
%% ========================================================================
if ismember(22, figures_to_plot)
figure('Name', 'Scatter_Matrix', 'Position', [100 50 1000 800], 'Color', 'w');

subplot(2, 2, 1);
scatter(T_uni_success.eta_SF, T_uni_success.PCE_pct, 40, T_uni_success.SF_flag, ...
        'filled', 'MarkerFaceAlpha', 0.6);
xlabel('\eta_{SF}', 'FontSize', 14);
ylabel('PCE (%)', 'FontSize', 14);
title('PCE vs \eta_{SF}', 'FontSize', 14);
colormap(gca, [CB_ORANGE; CB_BLUE]);
set(gca, 'FontSize', 12, 'Box', 'on'); grid on;

subplot(2, 2, 2);
scatter(log10(T_uni_success.SRV_cm_s), T_uni_success.PCE_pct, 40, ...
        T_uni_success.eta_SF, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('log_{10}(SRV)', 'FontSize', 14);
ylabel('PCE (%)', 'FontSize', 14);
title('PCE vs SRV', 'FontSize', 14);
colormap(gca, blueorange(256));
cb = colorbar; ylabel(cb, '\eta_{SF}');
set(gca, 'FontSize', 12, 'Box', 'on'); grid on;

subplot(2, 2, 3);
scatter(T_uni_success.junction_depth_nm, T_uni_success.PCE_pct, 40, ...
        T_uni_success.eta_SF, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('Junction Depth (nm)', 'FontSize', 14);
ylabel('PCE (%)', 'FontSize', 14);
title('PCE vs Junction Depth', 'FontSize', 14);
colormap(gca, blueorange(256));
cb = colorbar; ylabel(cb, '\eta_{SF}');
set(gca, 'FontSize', 12, 'Box', 'on'); grid on;

subplot(2, 2, 4);
scatter(T_uni_success.EF_peak_eV, T_uni_success.PCE_pct, 40, ...
        log10(T_uni_success.SRV_cm_s), 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('E_F (eV)', 'FontSize', 14);
ylabel('PCE (%)', 'FontSize', 14);
title('PCE vs Fermi Level', 'FontSize', 14);
colormap(gca, blueorange(256));
cb = colorbar; ylabel(cb, 'log_{10}(SRV)');
set(gca, 'FontSize', 12, 'Box', 'on'); grid on;

end

%% ========================================================================
%  FIGURE 23: Delta PCE Heatmap (SF - Baseline, uniform only)
%% ========================================================================
if ismember(23, figures_to_plot)
figure('Name', 'Delta_PCE_Heatmap', 'Position', [150 100 600 450], 'Color', 'w');

EF_vals_23  = unique(T_uni_success.EF_peak_eV);
SRV_vals_23 = unique(T_uni_success.SRV_cm_s);

ND_axis_23 = NaN(size(EF_vals_23));
for ii = 1:numel(EF_vals_23)
    [~, idx_map] = min(abs(EF_map - EF_vals_23(ii)));
    ND_axis_23(ii) = ND_map(idx_map);
end

delta_PCE_23 = NaN(length(EF_vals_23), length(SRV_vals_23));
for i = 1:length(EF_vals_23)
    for j = 1:length(SRV_vals_23)
        mask_base = (T_uni_success.SF_flag == 0) & ...
                    (T_uni_success.SRV_cm_s == SRV_vals_23(j)) & ...
                    (abs(T_uni_success.EF_peak_eV - EF_vals_23(i)) < tol) & ...
                    (T_uni_success.junction_depth_nm == 300) & ...
                    (abs(T_uni_success.Phi_R_offset_eV - 0) < tol);
        mask_sf = (T_uni_success.SF_flag == 1) & ...
                  (abs(T_uni_success.eta_SF - 2) < tol) & ...
                  (T_uni_success.SRV_cm_s == SRV_vals_23(j)) & ...
                  (abs(T_uni_success.EF_peak_eV - EF_vals_23(i)) < tol) & ...
                  (T_uni_success.junction_depth_nm == 300) & ...
                  (abs(T_uni_success.Phi_R_offset_eV - 0) < tol);
        if any(mask_base) && any(mask_sf)
            PCE_base = mean(T_uni_success.PCE_pct(mask_base), 'omitnan');
            PCE_sf   = mean(T_uni_success.PCE_pct(mask_sf),   'omitnan');
            delta_PCE_23(i, j) = PCE_sf - PCE_base;
        end
    end
end

h = imagesc(1:length(SRV_vals_23), 1:length(EF_vals_23), delta_PCE_23);
colormap(blueorange(256));
set(h, 'AlphaData', ~isnan(delta_PCE_23));
set(gca, 'Color', 'w');

cb = colorbar;
ylabel(cb, '\Delta PCE (%) [SF - Baseline]', 'FontSize', 18);
set(gca, 'YDir', 'normal');

set(gca, 'XTick', 1:length(SRV_vals_23), ...
         'XTickLabel', arrayfun(@(x) sprintf('%.0e', x), SRV_vals_23, 'UniformOutput', false));
set(gca, 'YTick', 1:length(EF_vals_23), ...
         'YTickLabel', arrayfun(@(x) sprintf('%.1e', x), ND_axis_23, 'UniformOutput', false));

xlabel('Surface recombination velocity s_p (cm/s)', 'FontSize', 20);
ylabel('n-emitter doping concentration N_D (cm^{-3})', 'FontSize', 20);
set(gca, 'FontSize', 16);

for i = 1:length(EF_vals_23)
    for j = 1:length(SRV_vals_23)
        if ~isnan(delta_PCE_23(i,j))
            if delta_PCE_23(i,j) > 0
                text(j, i, sprintf('+%.2f', delta_PCE_23(i,j)), ...
                     'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', 'k');
            else
                text(j, i, sprintf('%.2f', delta_PCE_23(i,j)), ...
                     'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', 'w');
            end
        end
    end
end

end

%% ========================================================================
%  FIGURE 24: Radar Chart - Best vs Worst (uniform only)
%% ========================================================================
if ismember(24, figures_to_plot)
figure('Name', 'Radar_Chart', 'Position', [200 150 600 550], 'Color', 'w');

[~, best_idx_24]  = max(T_uni_success.PCE_pct);
[~, worst_idx_24] = min(T_uni_success.PCE_pct);

metrics_24 = {'PCE', 'V_{oc}', 'J_{sc}', 'FF', 'Blue EQE'};
n_metrics_24 = length(metrics_24);

best_run_24  = runs(T_uni_success.run_number(best_idx_24));
worst_run_24 = runs(T_uni_success.run_number(worst_idx_24));

wl_24 = best_run_24.EQE.wavelengths;
blue_mask_24 = (wl_24 >= 405) & (wl_24 <= 500);
best_blue_eqe_24  = mean(best_run_24.EQE.EQE(blue_mask_24), 'omitnan');
worst_blue_eqe_24 = mean(worst_run_24.EQE.EQE(blue_mask_24), 'omitnan');

PCE_range_24     = [min(T_uni_success.PCE_pct), max(T_uni_success.PCE_pct)];
Voc_range_24     = [min(T_uni_success.Voc_V),   max(T_uni_success.Voc_V)];
Jsc_range_24     = [min(T_uni_success.Jsc_mA_cm2), max(T_uni_success.Jsc_mA_cm2)];
FF_range_24      = [min(T_uni_success.FF_pct),  max(T_uni_success.FF_pct)];
BlueEQE_range_24 = [50, 130];

best_norm_24 = [(T_uni_success.PCE_pct(best_idx_24) - PCE_range_24(1)) / diff(PCE_range_24), ...
                (T_uni_success.Voc_V(best_idx_24)   - Voc_range_24(1)) / diff(Voc_range_24), ...
                (T_uni_success.Jsc_mA_cm2(best_idx_24) - Jsc_range_24(1)) / diff(Jsc_range_24), ...
                (T_uni_success.FF_pct(best_idx_24) - FF_range_24(1)) / diff(FF_range_24), ...
                (best_blue_eqe_24 - BlueEQE_range_24(1)) / diff(BlueEQE_range_24)];

worst_norm_24 = [(T_uni_success.PCE_pct(worst_idx_24) - PCE_range_24(1)) / diff(PCE_range_24), ...
                 (T_uni_success.Voc_V(worst_idx_24)   - Voc_range_24(1)) / diff(Voc_range_24), ...
                 (T_uni_success.Jsc_mA_cm2(worst_idx_24) - Jsc_range_24(1)) / diff(Jsc_range_24), ...
                 (T_uni_success.FF_pct(worst_idx_24) - FF_range_24(1)) / diff(FF_range_24), ...
                 (worst_blue_eqe_24 - BlueEQE_range_24(1)) / diff(BlueEQE_range_24)];

best_norm_24  = [best_norm_24, best_norm_24(1)];
worst_norm_24 = [worst_norm_24, worst_norm_24(1)];
angles_24 = linspace(0, 2*pi, n_metrics_24 + 1);

polarplot(angles_24, best_norm_24, '-', 'LineWidth', 2.5, 'Color', CB_BLUE, ...
          'DisplayName', 'Best config');
hold on;
polarplot(angles_24, worst_norm_24, '-', 'LineWidth', 2.5, 'Color', CB_ORANGE, ...
          'DisplayName', 'Worst config');

ax24 = gca;
ax24.ThetaTick = rad2deg(angles_24(1:end-1));
ax24.ThetaTickLabel = metrics_24;
ax24.RLim = [0 1.1];
ax24.FontSize = 12;
title('Best vs Worst (Uniform Profile)', 'FontSize', 14);
legend('Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 12);
hold off;

end

%% ========================================================================
%  FIGURE 25: Stacked Area - EQE Decomposition (uniform only)
%% ========================================================================
if ismember(25, figures_to_plot)
figure('Name', 'EQE_Stacked', 'Position', [250 200 800 500], 'Color', 'w');

mask_base_25 = (T_uni_success.SF_flag == 0) & ...
               (T_uni_success.SRV_cm_s == optimal_SRV) & ...
               (abs(T_uni_success.EF_peak_eV - optimal_EF) < tol) & ...
               (T_uni_success.junction_depth_nm == optimal_jd) & ...
               (abs(T_uni_success.Phi_R_offset_eV - 0) < tol);

mask_sf_25 = (T_uni_success.SF_flag == 1) & ...
             (abs(T_uni_success.eta_SF - 2) < tol) & ...
             (T_uni_success.SRV_cm_s == optimal_SRV) & ...
             (abs(T_uni_success.EF_peak_eV - optimal_EF) < tol) & ...
             (T_uni_success.junction_depth_nm == optimal_jd) & ...
             (abs(T_uni_success.Phi_R_offset_eV - 0) < tol);

if any(mask_base_25) && any(mask_sf_25)
    base_run_25 = runs(T_uni_success.run_number(find(mask_base_25, 1)));
    sf_run_25   = runs(T_uni_success.run_number(find(mask_sf_25, 1)));

    wl25 = base_run_25.EQE.wavelengths;
    eqe_base_25 = base_run_25.EQE.EQE;
    eqe_sf_25   = sf_run_25.EQE.EQE;

    area(wl25, eqe_base_25, 'FaceColor', [0.50 0.70 0.86], 'FaceAlpha', 0.6, ...
         'EdgeColor', 'none', 'DisplayName', 'Baseline Si');
    hold on;
    eqe_diff_25 = max(0, eqe_sf_25 - eqe_base_25);
    area(wl25, eqe_base_25 + eqe_diff_25, 'FaceColor', CB_ORANGE, 'FaceAlpha', 0.6, ...
         'EdgeColor', 'none', 'DisplayName', 'SF Enhancement');
    plot(wl25, eqe_sf_25, '-', 'LineWidth', 2, 'Color', CB_ORANGE, ...
         'DisplayName', 'Total (SF, \eta=2)');
    plot([400 800], [100 100], 'k--', 'LineWidth', 1.5, 'DisplayName', '100% EQE');

    xlabel('Wavelength (nm)', 'FontSize', 16);
    ylabel('EQE (%)', 'FontSize', 16);
    title('EQE Decomposition (Uniform Profile)', 'FontSize', 18);
    xlim([400 700]); ylim([0 150]);
    legend('Location', 'northeast', 'FontSize', 12);
    set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'on');
    grid on;
    hold off;
end

end

%% ========================================================================
%  FIGURE 26: Summary Table (uniform only)
%% ========================================================================
if ismember(26, figures_to_plot)
figure('Name', 'Summary_Table', 'Position', [300 250 700 400], 'Color', 'w');
axis off;

[max_PCE_26, best_idx_26] = max(T_uni_success.PCE_pct);
best_26 = T_uni_success(best_idx_26, :);

n_uni = height(T_uni_success);
summary_text_26 = {
    'BEST CONFIGURATION (Uniform Profile):'
    sprintf('  PCE = %.3f %%', best_26.PCE_pct)
    sprintf('  Voc = %.4f V', best_26.Voc_V)
    sprintf('  Jsc = %.4f mA/cm^2', best_26.Jsc_mA_cm2)
    sprintf('  FF = %.2f %%', best_26.FF_pct)
    ''
    'PARAMETERS:'
    sprintf('  SF flag = %d', best_26.SF_flag)
    sprintf('  eta_SF = %.1f', best_26.eta_SF)
    sprintf('  Junction depth = %d nm', best_26.junction_depth_nm)
    sprintf('  SRV = %.0e cm/s', best_26.SRV_cm_s)
    sprintf('  E_F = %.3f eV', best_26.EF_peak_eV)
    sprintf('  Phi_R offset = %.4f eV', best_26.Phi_R_offset_eV)
    ''
    'SWEEP STATISTICS (uniform profile):'
    sprintf('  Uniform runs: %d', n_uni)
    sprintf('  Mean PCE (all): %.3f %%', mean(T_uni_success.PCE_pct, 'omitnan'))
    sprintf('  Mean PCE (SF=1, eta=2): %.3f %%', ...
            mean(T_uni_success.PCE_pct(T_uni_success.SF_flag==1 & T_uni_success.eta_SF==2), 'omitnan'))
};

text(0.1, 0.95, summary_text_26, 'VerticalAlignment', 'top', 'FontSize', 11, ...
     'FontName', 'FixedWidth');
title('Parameter Sweep Summary (Uniform Profile)', 'FontSize', 18);

end

%% ========================================================================
%  FIGURE 27: Best Device EQE (uniform only)
%% ========================================================================
if ismember(27, figures_to_plot)
figure('Name', 'Best_Device_EQE', 'Position', [350 300 800 550], 'Color', 'w');

[~, best_idx_27] = max(T_uni_success.PCE_pct);
best_27 = T_uni_success(best_idx_27, :);
best_run_27 = runs(best_27.run_number);

wl27  = best_run_27.EQE.wavelengths;
eqe27 = best_run_27.EQE.EQE;
iqe27 = best_run_27.EQE.IQE;

plot(wl27, eqe27, '-', 'LineWidth', 2.5, 'Color', CB_BLUE, 'DisplayName', 'EQE');
hold on;
plot(wl27, iqe27, '--', 'LineWidth', 2, 'Color', CB_ORANGE, 'DisplayName', 'IQE');
plot([400 800], [100 100], 'k--', 'LineWidth', 1.5, 'DisplayName', '100%');

fill([400 550 550 400], [0 0 160 160], [0.85 0.92 1.0], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off');

idx_above_27 = eqe27 > 100;
if any(idx_above_27)
    wl_above = wl27(idx_above_27);
    eqe_above = eqe27(idx_above_27);
    fill([wl_above, fliplr(wl_above)], [100*ones(size(wl_above)), fliplr(eqe_above)], ...
         [0.85 0.92 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
end

xlabel('Wavelength (nm)', 'FontSize', 18);
ylabel('Quantum Efficiency (%)', 'FontSize', 18);
title(sprintf('Best Device EQE (PCE = %.3f%%, Uniform)', best_27.PCE_pct), 'FontSize', 20);
xlim([400 800]); ylim([0 160]);
legend('Location', 'northeast', 'FontSize', 14);
set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'on');
grid on;

[max_eqe_27, max_idx_27] = max(eqe27);
text(wl27(max_idx_27), max_eqe_27 + 5, sprintf('Peak: %.1f%%', max_eqe_27), ...
     'FontSize', 12, 'HorizontalAlignment', 'center', 'Color', CB_BLUE);

param_str_27 = sprintf(['SF = %d, \\eta_{SF} = %.0f\n' ...
                         'j_d = %d nm\n' ...
                         'SRV = %.0e cm/s\n' ...
                         'E_F = %.3f eV'], ...
                        best_27.SF_flag, best_27.eta_SF, ...
                        best_27.junction_depth_nm, ...
                        best_27.SRV_cm_s, best_27.EF_peak_eV);
annotation('textbox', [0.15, 0.65, 0.25, 0.2], 'String', param_str_27, ...
           'FontSize', 11, 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
           'FitBoxToText', 'on', 'Interpreter', 'tex');
hold off;

end

%% ========================================================================
%  FIGURE 28: Best Device J-V (uniform only)
%% ========================================================================
if ismember(28, figures_to_plot)
figure('Name', 'Best_Device_JV', 'Position', [400 350 800 550], 'Color', 'w');

[~, best_idx_28] = max(T_uni_success.PCE_pct);
best_28     = T_uni_success(best_idx_28, :);
best_run_28 = runs(best_28.run_number);

V28 = best_run_28.JV.Vapp;
J28 = best_run_28.JV.Jtot;
[~, idx_max_28] = max(V28);
V_fwd_28 = V28(1:idx_max_28);
J_fwd_28 = J28(1:idx_max_28);

plot(V_fwd_28, J_fwd_28, '-', 'LineWidth', 2.5, 'Color', CB_BLUE);
hold on;
plot(xlim, [0 0], 'k--', 'LineWidth', 1);

% Jsc
[~, idx_v0_28] = min(abs(V_fwd_28));
Jsc_plot_28 = J_fwd_28(idx_v0_28);
plot(0, Jsc_plot_28, 'o', 'MarkerSize', 10, 'MarkerFaceColor', CB_ORANGE, ...
     'MarkerEdgeColor', CB_ORANGE);
text(0.02, Jsc_plot_28, sprintf('J_{sc} = %.4f mA/cm^2', abs(Jsc_plot_28)), ...
     'FontSize', 11, 'VerticalAlignment', 'bottom');

% Voc
Voc_plot_28 = best_28.Voc_V;
plot(Voc_plot_28, 0, 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0.55 0.34 0.64], ...
     'MarkerEdgeColor', [0.55 0.34 0.64]);
text(Voc_plot_28, 0.02, sprintf('V_{oc} = %.4f V', Voc_plot_28), ...
     'FontSize', 11, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% MPP
power_28 = V_fwd_28 .* J_fwd_28;
[MPP_28, idx_mpp_28] = min(power_28);
V_mpp_28 = V_fwd_28(idx_mpp_28(1));
J_mpp_28 = J_fwd_28(idx_mpp_28(1));
plot(V_mpp_28, J_mpp_28, 's', 'MarkerSize', 12, ...
     'MarkerFaceColor', [0.55 0.34 0.64], 'MarkerEdgeColor', [0.55 0.34 0.64]);
text(V_mpp_28 + 0.02, J_mpp_28, sprintf('MPP\nP = %.4f mW/cm^2', abs(MPP_28(1))), ...
     'FontSize', 10, 'VerticalAlignment', 'middle');

fill([0 V_mpp_28 V_mpp_28 0], [0 0 J_mpp_28 J_mpp_28], [0.85 0.92 1.0], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.3);

xlabel('Voltage (V)', 'FontSize', 18);
ylabel('Current Density (mA/cm^2)', 'FontSize', 18);
title(sprintf('Best Device J-V (PCE = %.3f%%, FF = %.2f%%, Uniform)', ...
              best_28.PCE_pct, best_28.FF_pct), 'FontSize', 20);
set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'on');
grid on;

metrics_str_28 = sprintf(['PCE = %.3f %%\n' ...
                           'V_{oc} = %.4f V\n' ...
                           'J_{sc} = %.4f mA/cm^2\n' ...
                           'FF = %.2f %%'], ...
                          best_28.PCE_pct, best_28.Voc_V, ...
                          best_28.Jsc_mA_cm2, best_28.FF_pct);
annotation('textbox', [0.15, 0.15, 0.25, 0.2], 'String', metrics_str_28, ...
           'FontSize', 12, 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
           'FitBoxToText', 'on', 'Interpreter', 'tex');
hold off;

end

%% ========================================================================
%  FIGURE 29: Best Device vs Baseline Comparison (uniform only)
%% ========================================================================
if ismember(29, figures_to_plot)
figure('Name', 'Best_vs_Baseline', ...
       'Position', [100 100 1400 500], 'Color', 'w');

[~, best_idx_29] = max(T_uni_success.PCE_pct);
best_29     = T_uni_success(best_idx_29, :);
best_run_29 = runs(best_29.run_number);

mask_baseline_29 = (T_uni_success.SF_flag == 0) & ...
                   (T_uni_success.junction_depth_nm == best_29.junction_depth_nm) & ...
                   (T_uni_success.SRV_cm_s == best_29.SRV_cm_s) & ...
                   (abs(T_uni_success.EF_peak_eV - best_29.EF_peak_eV) < 0.01) & ...
                   (abs(T_uni_success.Phi_R_offset_eV - best_29.Phi_R_offset_eV) < 0.01);

if any(mask_baseline_29)
    baseline_29     = T_uni_success(find(mask_baseline_29, 1), :);
    baseline_run_29 = runs(baseline_29.run_number);

    %% Panel A: EQE comparison
    subplot(1, 2, 1);
    hold on;

    h_sf_eqe_29 = plot(best_run_29.EQE.wavelengths, best_run_29.EQE.EQE, ...
                       '-', 'Color', CB_ORANGE, 'LineWidth', 2.5);
    h_bl_eqe_29 = plot(baseline_run_29.EQE.wavelengths, baseline_run_29.EQE.EQE, ...
                       '-', 'Color', CB_BLUE, 'LineWidth', 2.5);

    plot([400 800], [100 100], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlabel('Wavelength (nm)', 'FontSize', 24);
    ylabel('EQE (%)', 'FontSize', 24);
    xlim([400 800]); 
    %ylim([80 110]);
    set(gca, 'FontSize', 18, 'LineWidth', 1.2, 'Box', 'on');
    grid on;

    [eqe_max_29, idx_max_29] = max(best_run_29.EQE.EQE);
    lambda_max_29 = best_run_29.EQE.wavelengths(idx_max_29);
    plot(lambda_max_29, eqe_max_29, 'o', 'MarkerSize', 8, ...
         'MarkerFaceColor', CB_ORANGE, 'MarkerEdgeColor', [0.65 0.25 0.00], ...
         'HandleVisibility', 'off');
    yl29 = ylim;
    text(lambda_max_29 + 5, eqe_max_29 + 0.02*(yl29(2)-yl29(1)), ...
         sprintf('%.1f%% @ %dnm', eqe_max_29, round(lambda_max_29)), ...
         'FontSize', 18, 'Color', [0.65 0.25 0.00], ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

    legend([h_bl_eqe_29, h_sf_eqe_29], ...
           {'Baseline', sprintf('SF (\\eta=%.0f)', best_29.eta_SF)}, ...
           'Location', 'northeast', 'FontSize', 16);
    hold off;

    %% Panel B: J-V comparison
    ax2_29 = subplot(1, 2, 2); hold on;

    V_best_29 = best_run_29.JV.Vapp;
    J_best_29 = best_run_29.JV.Jtot;
    [~, idx_max_best_29] = max(V_best_29);
    h_sf_29 = plot(V_best_29(1:idx_max_best_29), J_best_29(1:idx_max_best_29), ...
                   '-', 'Color', CB_ORANGE, 'LineWidth', 2.5);

    V_base_29 = baseline_run_29.JV.Vapp;
    J_base_29 = baseline_run_29.JV.Jtot;
    [~, idx_max_base_29] = max(V_base_29);
    h_bl_29 = plot(V_base_29(1:idx_max_base_29), J_base_29(1:idx_max_base_29), ...
                   '-', 'Color', CB_BLUE, 'LineWidth', 2.5);

    plot(xlim, [0 0], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    yl29b = ylim;
    plot([0 0], yl29b, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

    xlabel('Voltage (V)', 'FontSize', 24);
    ylabel('Current Density (mA/cm^2)', 'FontSize', 24);
    xlim([-0.3, 0.8]); ylim([-0.05, 0.1]);
    set(gca, 'FontSize', 18, 'LineWidth', 1.2, 'Box', 'on');
    grid on;

    base_line1_29 = sprintf('Baseline: PCE = %.3f%%            V_{oc} = %.4f V', ...
                             baseline_29.PCE_pct, baseline_29.Voc_V);
    base_line2_29 = sprintf('                J_{sc} = %.4f mA/cm^2   FF = %.2f%%', ...
                             baseline_29.Jsc_mA_cm2, baseline_29.FF_pct);
    sf_line1_29   = sprintf('SF (\\eta=%.0f): PCE = %.3f%%            V_{oc} = %.4f V', ...
                             best_29.eta_SF, best_29.PCE_pct, best_29.Voc_V);
    sf_line2_29   = sprintf('                J_{sc} = %.4f mA/cm^2   FF = %.2f%%', ...
                             best_29.Jsc_mA_cm2, best_29.FF_pct);

    h_dummy1_29 = plot(NaN, NaN, 'w', 'LineStyle', 'none', 'Marker', 'none');
    h_dummy2_29 = plot(NaN, NaN, 'w', 'LineStyle', 'none', 'Marker', 'none');
    legend([h_bl_29, h_dummy1_29, h_sf_29, h_dummy2_29], ...
           {base_line1_29, base_line2_29, sf_line1_29, sf_line2_29}, ...
           'Location', 'northwest', 'FontSize', 16, 'Interpreter', 'tex');
    hold off;
else
    text(0.5, 0.5, 'No matching baseline found', 'FontSize', 16, ...
         'HorizontalAlignment', 'center');
end

end

%% ========================================================================
%  FIGURE 30: PCE vs Junction Depth (SF vs Baseline, uniform only)
%% ========================================================================
if ismember(30, figures_to_plot)
figure('Name', 'PCE_vs_JD', 'Position', [200 200 600 400], 'Color', 'w');
hold on;

eta_val_30 = 2;

% SF (eta=2)
PCE_vs_jd_30 = NaN(length(jd_vals_all), 1);
for j = 1:length(jd_vals_all)
    mask = (T_uni_success.SF_flag == 1) & ...
           (abs(T_uni_success.eta_SF - eta_val_30) < tol) & ...
           (T_uni_success.junction_depth_nm == jd_vals_all(j)) & ...
           (T_uni_success.SRV_cm_s == optimal_SRV) & ...
           (abs(T_uni_success.EF_peak_eV - optimal_EF) < tol) & ...
           (abs(T_uni_success.Phi_R_offset_eV - 0) < tol);
    if any(mask)
        PCE_vs_jd_30(j) = mean(T_uni_success.PCE_pct(mask), 'omitnan');
    end
end

h_sf_30 = plot(jd_vals_all, PCE_vs_jd_30, '-o', ...
        'Color', CB_BLUE, 'LineWidth', 2, 'MarkerSize', 8, ...
        'MarkerFaceColor', CB_BLUE, 'DisplayName', '\eta_{SF} = 2');

% Baseline (SF=0)
PCE_base_jd_30 = NaN(length(jd_vals_all), 1);
for j = 1:length(jd_vals_all)
    mask = (T_uni_success.SF_flag == 0) & ...
           (T_uni_success.junction_depth_nm == jd_vals_all(j)) & ...
           (T_uni_success.SRV_cm_s == optimal_SRV) & ...
           (abs(T_uni_success.EF_peak_eV - optimal_EF) < tol) & ...
           (abs(T_uni_success.Phi_R_offset_eV - 0) < tol);
    if any(mask)
        PCE_base_jd_30(j) = mean(T_uni_success.PCE_pct(mask), 'omitnan');
    end
end

h_base_30 = plot(jd_vals_all, PCE_base_jd_30, '--^', ...
           'Color', [0.50 0.50 0.50], 'LineWidth', 2, 'MarkerSize', 8, ...
           'MarkerFaceColor', [0.50 0.50 0.50], 'DisplayName', 'Baseline (no SF)');

% Linear fits & slope comparison
valid_sf   = ~isnan(PCE_vs_jd_30);
valid_base = ~isnan(PCE_base_jd_30);
if sum(valid_sf) >= 2 && sum(valid_base) >= 2
    p_sf   = polyfit(jd_vals_all(valid_sf),   PCE_vs_jd_30(valid_sf),       1);
    p_base = polyfit(jd_vals_all(valid_base), PCE_base_jd_30(valid_base), 1);

    fprintf('FIG 30 slopes: baseline = %.3g %%/nm,  eta=2 = %.3g %%/nm\n', ...
            p_base(1), p_sf(1));

    all_pce = [PCE_vs_jd_30; PCE_base_jd_30];
    y_txt = nanmin(all_pce) + 0.5 * (nanmax(all_pce) - nanmin(all_pce));
    txt_str = sprintf(['Slope (per 100nm):\\newline', ...
                       'Baseline: %.4f%% / 100nm\\newline', ...
                       '\\eta_{SF}=2: %.4f%% / 100nm'], ...
                       100*p_base(1), 100*p_sf(1));
    text(min(jd_vals_all)+5, y_txt, txt_str, ...
         'FontSize', 14, 'BackgroundColor', 'w', ...
         'EdgeColor', [0.5 0.5 0.5], 'Margin', 6);
end

xlabel('Junction depth (nm)', 'FontSize', 20);
ylabel('PCE (%)', 'FontSize', 20);
legend([h_sf_30, h_base_30], 'Location', 'best', 'FontSize', 16);
set(gca, 'FontSize', 16, 'LineWidth', 1.2, 'Box', 'on');
xlim([min(jd_vals_all)-20, max(jd_vals_all)+20]);
grid on;
title('PCE vs Junction Depth (Uniform Profile)', 'FontSize', 16);
hold off;

end

%% ========================================================================
%  FIGURE 31: EQE(520nm) vs Junction Depth (SF vs Baseline, uniform only)
%% ========================================================================
if ismember(31, figures_to_plot)
figure('Name', 'EQE520_vs_JD', 'Position', [200 200 600 400], 'Color', 'w');
hold on;

eta_val_31  = 2;
lambda_pk31 = 520;

EQE520_sf_31   = NaN(length(jd_vals_all), 1);
EQE520_base_31 = NaN(length(jd_vals_all), 1);

for j = 1:length(jd_vals_all)
    jd = jd_vals_all(j);

    mask_base = (T_uni_success.SF_flag == 0) & ...
                (T_uni_success.junction_depth_nm == jd) & ...
                (T_uni_success.SRV_cm_s == optimal_SRV) & ...
                (abs(T_uni_success.EF_peak_eV - optimal_EF) < tol) & ...
                (abs(T_uni_success.Phi_R_offset_eV - 0) < tol);
    if any(mask_base)
        run_idx_list = T_uni_success.run_number(mask_base);
        EQE520_base_31(j) = local_mean_eqe_at_lambda(runs, run_idx_list, lambda_pk31);
    end

    mask_sf = (T_uni_success.SF_flag == 1) & ...
              (abs(T_uni_success.eta_SF - eta_val_31) < tol) & ...
              (T_uni_success.junction_depth_nm == jd) & ...
              (T_uni_success.SRV_cm_s == optimal_SRV) & ...
              (abs(T_uni_success.EF_peak_eV - optimal_EF) < tol) & ...
              (abs(T_uni_success.Phi_R_offset_eV - 0) < tol);
    if any(mask_sf)
        run_idx_list = T_uni_success.run_number(mask_sf);
        EQE520_sf_31(j) = local_mean_eqe_at_lambda(runs, run_idx_list, lambda_pk31);
    end
end

h_sf_31 = plot(jd_vals_all, EQE520_sf_31, '-o', ...
    'Color', CB_BLUE, 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', CB_BLUE, 'DisplayName', '\eta_{SF} = 2');

h_base_31 = plot(jd_vals_all, EQE520_base_31, '--^', ...
    'Color', [0.50 0.50 0.50], 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', [0.50 0.50 0.50], 'DisplayName', 'Baseline (no SF)');

xlabel('Junction depth (nm)', 'FontSize', 20);
ylabel('EQE at 520nm (%)', 'FontSize', 20);
legend([h_sf_31, h_base_31], 'Location', 'best', 'FontSize', 16);
set(gca, 'FontSize', 16, 'LineWidth', 1.2, 'Box', 'on');
xlim([min(jd_vals_all)-20, max(jd_vals_all)+20]);
grid on;

all_vals_31 = [EQE520_sf_31(:); EQE520_base_31(:)];
all_vals_31 = all_vals_31(~isnan(all_vals_31));
if ~isempty(all_vals_31)
    yrng = max(all_vals_31) - min(all_vals_31);
    if yrng == 0, yrng = max(1, 0.1*max(all_vals_31)); end
    ylim([min(all_vals_31) - 0.15*yrng, max(all_vals_31) + 0.15*yrng]);
end

valid_sf   = ~isnan(EQE520_sf_31);
valid_base = ~isnan(EQE520_base_31);
if sum(valid_sf) >= 2 && sum(valid_base) >= 2
    p_sf   = polyfit(jd_vals_all(valid_sf),   EQE520_sf_31(valid_sf),      1);
    p_base = polyfit(jd_vals_all(valid_base), EQE520_base_31(valid_base),  1);

    fprintf('FIG 31 (EQE 520 nm) slopes: baseline = %.3g %%/nm,  eta=2 = %.3g %%/nm\n', ...
            p_base(1), p_sf(1));

    y_txt = min(all_vals_31) + 0.5 * (max(all_vals_31) - min(all_vals_31));
    txt_str = sprintf(['Slope:\\newline', ...
                       '\\eta_{SF}=2: %.2f%% / 100nm\\newline',...
                        'Baseline: %.2f%% / 100nm'], ...
                       100*p_sf(1), 100*p_base(1));
    text(min(jd_vals_all)+5, y_txt, txt_str, ...
         'FontSize', 14, 'BackgroundColor', 'w', ...
         'EdgeColor', [0.5 0.5 0.5], 'Margin', 6);
end

title('EQE(520nm) vs Junction Depth (Uniform Profile)', 'FontSize', 16);
hold off;

end

fprintf('\n========================================\n');
fprintf('Generated %d figures!\n', length(figures_to_plot));
fprintf('========================================\n');

end  % end main function


%% ========================================================================
%  HELPER FUNCTIONS
%% ========================================================================
function [PCE_ref, PCE_test] = match_profiles(T_ref, T_test, tol)
%MATCH_PROFILES  Find matched (SF, eta, EF, SRV, PhiOff) pairs

    PCE_ref  = [];
    PCE_test = [];

    for i = 1:height(T_ref)
        mask = (T_test.SF_flag == T_ref.SF_flag(i)) & ...
               (abs(T_test.eta_SF - T_ref.eta_SF(i)) < tol) & ...
               (abs(T_test.EF_peak_eV - T_ref.EF_peak_eV(i)) < tol) & ...
               (T_test.SRV_cm_s == T_ref.SRV_cm_s(i)) & ...
               (abs(T_test.Phi_R_offset_eV - T_ref.Phi_R_offset_eV(i)) < tol);
        idx = find(mask, 1);
        if ~isempty(idx)
            PCE_ref(end+1)  = T_ref.PCE_pct(i); %#ok<AGROW>
            PCE_test(end+1) = T_test.PCE_pct(idx); %#ok<AGROW>
        end
    end

    PCE_ref  = PCE_ref(:);
    PCE_test = PCE_test(:);
end


function c = redgreen(m)
%REDGREEN Red-white-green diverging colormap

    if nargin < 1, m = 256; end

    n1 = floor(m/2);
    n2 = m - n1;

    r1 = linspace(0.65, 1.0, n1);
    g1 = linspace(0.00, 1.0, n1);
    b1 = linspace(0.00, 1.0, n1);

    r2 = linspace(1.0, 0.00, n2);
    g2 = linspace(1.0, 0.50, n2);
    b2 = linspace(1.0, 0.00, n2);

    c = [r1(:), g1(:), b1(:); r2(:), g2(:), b2(:)];
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


function c = blueorange(m)
%BLUEORANGE  Colorblind-friendly diverging colormap (blue-white-orange)

    if nargin < 1, m = 256; end

    n1 = floor(m/2);
    n2 = m - n1;

    r1 = linspace(0.12, 1.0, n1);
    g1 = linspace(0.47, 1.0, n1);
    b1 = linspace(0.71, 1.0, n1);

    r2 = linspace(1.0, 0.85, n2);
    g2 = linspace(1.0, 0.37, n2);
    b2 = linspace(1.0, 0.01, n2);

    c = [r1(:), g1(:), b1(:); r2(:), g2(:), b2(:)];
end