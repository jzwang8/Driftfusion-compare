function T = postprocess_sweep_qsspc()
%POSTPROCESS_SWEEP_QSSPC Load per-run .mat files and assemble results table
%
%   Reads all QSSPC_*.mat files from the sweep output directory,
%   assembles a results table, and generates analysis figures.
%
%   Usage:
%     T = postprocess_sweep_qsspc();

%% Load config
config = sweep_config_qsspc();
outputDir = config.outputDir;

%% Find all result files
files = dir(fullfile(outputDir, 'QSSPC_*.mat'));
n_files = length(files);

if n_files == 0
    error('No QSSPC result files found in %s', outputDir);
end

fprintf('Found %d result files in %s\n', n_files, outputDir);

%% Load and assemble
sp_r_arr        = NaN(n_files, 1);
Qf_arr          = NaN(n_files, 1);
Vi_arr          = NaN(n_files, 1);
Phi_right_arr   = NaN(n_files, 1);
tau_eff_arr     = NaN(n_files, 1);
Delta_p_avg_arr = NaN(n_files, 1);
G_avg_arr       = NaN(n_files, 1);
p_surface_arr   = NaN(n_files, 1);
S_front_eff_arr = NaN(n_files, 1);
delta_QFL_arr   = NaN(n_files, 1);

for ii = 1:n_files
    data = load(fullfile(files(ii).folder, files(ii).name), 'thisRun');
    r = data.thisRun;

    sp_r_arr(ii)        = r.sp_r;
    Qf_arr(ii)          = r.Qf_cm2;
    Vi_arr(ii)          = r.Vi_V;
    Phi_right_arr(ii)   = r.Phi_right;
    tau_eff_arr(ii)     = r.tau_eff_s;
    Delta_p_avg_arr(ii) = r.Delta_p_avg;
    G_avg_arr(ii)       = r.G_avg;
    p_surface_arr(ii)   = r.p_surface;
    S_front_eff_arr(ii) = r.S_front_eff;
    delta_QFL_arr(ii)   = r.delta_QFL;
end

T = table(sp_r_arr, Qf_arr, Vi_arr, Phi_right_arr, ...
          tau_eff_arr, Delta_p_avg_arr, G_avg_arr, p_surface_arr, ...
          S_front_eff_arr, delta_QFL_arr, ...
    'VariableNames', {'sp_r', 'Qf_cm2', 'Vi_V', 'Phi_right_eV', ...
    'tau_eff_s', 'Delta_p_avg', 'G_avg', 'p_surface', ...
    'S_front_eff', 'delta_QFL_eV'});

% Sort by sp_r then Qf
T = sortrows(T, {'sp_r', 'Qf_cm2'});

fprintf('Assembled %d results\n', height(T));

%% Save assembled table
save(fullfile(outputDir, 'T_qsspc.mat'), 'T', 'config');
writetable(T, fullfile(outputDir, 'T_qsspc.csv'));
fprintf('Saved T_qsspc.mat and T_qsspc.csv\n');

%% ==================== FIGURES ====================

sp_r_values  = config.sp_r_values;
Qf_values    = config.Qf_values_cm2;
tau_bulk     = config.tau_bulk_s;
W            = config.wafer_thickness_cm;
sp_l         = config.sp_l;
n_sp = length(sp_r_values);
n_Qf = length(Qf_values);
colors = lines(n_sp);

%% Figure 1: tau_eff vs Qf for each sp_r
figure('Position', [100, 100, 700, 500]);
hold on;
for jj = 1:n_sp
    mask = T.sp_r == sp_r_values(jj);
    Tsub = sortrows(T(mask, :), 'Qf_cm2');
    plot(Tsub.Qf_cm2, Tsub.tau_eff_s * 1e6, '-o', ...
        'Color', colors(jj,:), 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'DisplayName', sprintf('s_p = %g cm/s', sp_r_values(jj)));
end
yline(tau_bulk * 1e6, '--k', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('\\tau_{bulk} = %.0f ms', tau_bulk * 1e3));
hold off;
set(gca, 'YScale', 'log');
xlabel('Q_f [cm^{-2}]');
ylabel('\tau_{eff} [\mus]');
title('Effective lifetime vs fixed charge density');
legend('Location', 'southeast', 'FontSize', 8);
grid on;
saveas(gcf, fullfile(outputDir, 'fig1_tau_vs_Qf.png'));

%% Figure 2: S_front_eff vs Qf
figure('Position', [100, 100, 700, 500]);
hold on;
for jj = 1:n_sp
    mask = T.sp_r == sp_r_values(jj);
    Tsub = sortrows(T(mask, :), 'Qf_cm2');
    S_plot = Tsub.S_front_eff;
    S_plot(S_plot <= 0) = NaN;  % bulk-limited => S undefined
    plot(Tsub.Qf_cm2, S_plot, '-o', ...
        'Color', colors(jj,:), 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'DisplayName', sprintf('s_p = %g cm/s', sp_r_values(jj)));
end
hold off;
set(gca, 'YScale', 'log');
xlabel('Q_f [cm^{-2}]');
ylabel('S_{eff,front} [cm/s]');
title('Effective front SRV vs Q_f');
legend('Location', 'best', 'FontSize', 8);
grid on;
saveas(gcf, fullfile(outputDir, 'fig2_Seff_vs_Qf.png'));

%% Figure 3: Heatmap tau_eff(sp_r, Qf)
figure('Position', [100, 100, 700, 550]);

tau_grid = NaN(n_Qf, n_sp);
for jj = 1:n_sp
    for kk = 1:n_Qf
        mask = T.sp_r == sp_r_values(jj) & T.Qf_cm2 == Qf_values(kk);
        if any(mask)
            tau_grid(kk, jj) = T.tau_eff_s(mask) * 1e6;
        end
    end
end

imagesc(1:n_sp, 1:n_Qf, tau_grid);
set(gca, 'YDir', 'normal');
cb = colorbar;
cb.Label.String = '\tau_{eff} [\mus]';
colormap(flipud(hot));

xticks(1:n_sp);
xticklabels(arrayfun(@(x) sprintf('%.0e', x), sp_r_values, 'UniformOutput', false));
xlabel('s_p [cm/s]');

yticks(1:n_Qf);
yticklabels(arrayfun(@(x) sprintf('%.1e', x), Qf_values, 'UniformOutput', false));
ylabel('Q_f [cm^{-2}]');

title(sprintf('\\tau_{eff} [\\mus] | n-Si %.0f \\Omega\\cdotcm', config.resistivity_ohmcm));

% Annotate
for kk = 1:n_Qf
    for jj = 1:n_sp
        if ~isnan(tau_grid(kk, jj))
            val = tau_grid(kk, jj);
            if val > 100
                txt = sprintf('%.0f', val);
            else
                txt = sprintf('%.1f', val);
            end
            text(jj, kk, txt, 'HorizontalAlignment', 'center', ...
                'FontSize', 7, 'Color', 'w');
        end
    end
end
saveas(gcf, fullfile(outputDir, 'fig3_heatmap.png'));

%% Figure 4: Light soaking overlay (placeholder experimental data)
% Replace these arrays with your actual QSSPC measurements
t_soak_min  = [0, 1, 2, 5, 10, 20, 30];
tau_exp_us  = [50, 80, 150, 250, 350, 450, 500];
tau_exp_s   = tau_exp_us * 1e-6;

% Choose sp_r for mapping (your estimate of chemical passivation quality)
sp_r_for_mapping = 1000;  % [cm/s]

mask_sp = T.sp_r == sp_r_for_mapping;
if ~any(mask_sp)
    fprintf('Warning: sp_r = %g not found in results. Skipping Fig 4.\n', sp_r_for_mapping);
else
    Tsub = sortrows(T(mask_sp, :), 'Qf_cm2');

    % Interpolate: experimental tau -> Qf
    Qf_mapped = NaN(size(tau_exp_s));
    for kk = 1:length(tau_exp_s)
        if tau_exp_s(kk) >= min(Tsub.tau_eff_s) && tau_exp_s(kk) <= max(Tsub.tau_eff_s)
            Qf_mapped(kk) = interp1(Tsub.tau_eff_s, Tsub.Qf_cm2, tau_exp_s(kk), 'pchip');
        end
    end

    figure('Position', [100, 100, 1000, 450]);

    subplot(1,2,1);
    hold on;
    plot(Tsub.Qf_cm2, Tsub.tau_eff_s * 1e6, '-k', 'LineWidth', 2, ...
        'DisplayName', sprintf('Simulation (s_p = %g)', sp_r_for_mapping));
    scatter(Qf_mapped, tau_exp_us, 80, t_soak_min, 'filled', 'MarkerEdgeColor', 'k');
    cb2 = colorbar; cb2.Label.String = 'Soak time [min]';
    hold off;
    xlabel('Q_f [cm^{-2}]');
    ylabel('\tau_{eff} [\mus]');
    title('Experimental \tau mapped to Q_f');
    legend('Location', 'northwest');
    grid on;

    subplot(1,2,2);
    yyaxis left;
    plot(t_soak_min, tau_exp_us, 'o-b', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    ylabel('\tau_{eff} [\mus]');
    set(gca, 'YColor', 'b');

    yyaxis right;
    valid = ~isnan(Qf_mapped);
    plot(t_soak_min(valid), Qf_mapped(valid), 's-r', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
    ylabel('Inferred Q_f [cm^{-2}]');
    set(gca, 'YColor', 'r');

    xlabel('Light soaking time [min]');
    title('Light soaking trajectory');
    grid on;

    saveas(gcf, fullfile(outputDir, 'fig4_light_soaking.png'));
end

%% Print summary
fprintf('\n========== SUMMARY ==========\n');
fprintf('Vi <-> Qf conversion (Ci = %.3e F/cm^2):\n', config.Ci_F_per_cm2);
for kk = 1:length(Qf_values)
    fprintf('  Qf = %+.2e  =>  Vi = %+.4f V\n', ...
        Qf_values(kk), config.Vi_values_V(kk));
end

fprintf('\nResults at Qf = 0:\n');
mask_Qf0 = T.Qf_cm2 == 0;
T_Qf0 = T(mask_Qf0, :);
for kk = 1:height(T_Qf0)
    fprintf('  sp_r = %6.0f => tau_eff = %8.2f us | S_front_eff = %8.1f cm/s\n', ...
        T_Qf0.sp_r(kk), T_Qf0.tau_eff_s(kk) * 1e6, T_Qf0.S_front_eff(kk));
end

fprintf('\n=== Done ===\n');

end