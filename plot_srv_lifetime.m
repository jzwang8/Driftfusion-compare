function plot_srv_lifetime()
%PLOT_SRV_LIFETIME  Two-panel figure showing the relationship between
%   surface recombination velocity (SRV) and minority carrier lifetime.
%
%   Left panel:  tau (µs) vs SRV (cm/s)  — x: SRV, y: lifetime
%   Right panel: SRV (cm/s) vs tau (µs)  — x: lifetime, y: SRV (inverted)
%
%   Two lifetime models shown:
%     Emitter curves:  tau = jd / S      (surface-limited, 200nm emitter)
%     Wafer curve:     tau = W / (2*S)   (both surfaces active, QSSPC/µ-PCD)
%
%   The wafer formula applies to your measured ~60 µs values.
%   The emitter formula applies to the drift-diffusion model parameter s_p.

close all;

%% ========================================================================
%  PARAMETERS
%% ========================================================================

% Emitter junction depths (nm) — drift-diffusion model
jd_vals_nm = [50, 100, 200, 300, 500];
jd_colors  = [0.85 0.33 0.10;   % brick red   - 50 nm
              0.93 0.69 0.13;   % amber       - 100 nm
              0.47 0.67 0.19;   % olive green - 200 nm
              0.30 0.60 0.74;   % steel blue  - 300 nm
              0.49 0.18 0.56];  % purple      - 500 nm

% Wafer parameters (QSSPC / µ-PCD measurement context)
W_wafer_um  = 280;              % wafer thickness (µm)
W_wafer_cm  = W_wafer_um * 1e-4;
wafer_color = [0.1 0.1 0.1];   % near-black

SRV_min = 1e0;
SRV_max = 1e6;
SRV_vec = logspace(log10(SRV_min), log10(SRV_max), 500);

% Passivation regime boundaries
regimes = {
    struct('lo', 1e0, 'hi', 1e2, 'label', 'Excellent', 'color', [0.80 0.95 0.80]);
    struct('lo', 1e2, 'hi', 1e4, 'label', 'Good',      'color', [0.95 0.95 0.80]);
    struct('lo', 1e4, 'hi', 1e6, 'label', 'Poor',      'color', [0.98 0.88 0.80]);
};

%% ========================================================================
%  COMPUTE LIFETIMES
%% ========================================================================

n_jd   = numel(jd_vals_nm);
tau_emitter = zeros(n_jd, numel(SRV_vec));
for k = 1:n_jd
    jd_cm            = jd_vals_nm(k) * 1e-7;
    tau_emitter(k,:) = (jd_cm ./ SRV_vec) * 1e6;   % µs
end

% Wafer: tau = W / (2S), both surfaces active
tau_wafer = (W_wafer_cm ./ (2 * SRV_vec)) * 1e6;   % µs

% Combined y-axis limits
tau_all_min = min([tau_emitter(:); tau_wafer(:)]);
tau_all_max = max([tau_emitter(:); tau_wafer(:)]);

%% ========================================================================
%  FIGURE
%% ========================================================================

figure('Name', 'SRV_Lifetime_TwoPanel', ...
       'Position', [80 100 1200 540], 'Color', 'w');

%% ----------------------------------------------------------------
%  LEFT PANEL: tau (y) vs SRV (x)
%% ----------------------------------------------------------------
ax1 = subplot(1, 2, 1);
hold(ax1, 'on');

% Shaded regimes
for ir = 1:numel(regimes)
    patch(ax1, ...
          [regimes{ir}.lo regimes{ir}.hi regimes{ir}.hi regimes{ir}.lo], ...
          [tau_all_min*0.1 tau_all_min*0.1 tau_all_max*10 tau_all_max*10], ...
          regimes{ir}.color, 'EdgeColor', 'none', 'FaceAlpha', 0.55, ...
          'HandleVisibility', 'off');
    text(ax1, sqrt(regimes{ir}.lo * regimes{ir}.hi), tau_all_max * 3, ...
         regimes{ir}.label, 'FontSize', 10, 'HorizontalAlignment', 'center', ...
         'Color', [0.25 0.25 0.25], 'FontWeight', 'bold');
end

% Emitter curves
h_emitter = gobjects(n_jd, 1);
for k = 1:n_jd
    h_emitter(k) = loglog(ax1, SRV_vec, tau_emitter(k,:), '-', ...
                          'LineWidth', 1.8, 'Color', jd_colors(k,:));
end

% Wafer curve
h_wafer = loglog(ax1, SRV_vec, tau_wafer, 'k-', ...
                 'LineWidth', 2.5, 'Color', wafer_color);

xlabel(ax1, 'Surface recombination velocity s_p (cm/s)', 'FontSize', 15);
ylabel(ax1, 'Minority carrier lifetime \tau_p (\mus)', 'FontSize', 15);
title(ax1, '\tau_p vs s_p', 'FontSize', 15);
set(ax1, 'XScale', 'log', 'YScale', 'log', ...
         'FontSize', 13, 'LineWidth', 1.2, 'Box', 'on');
xlim(ax1, [SRV_min, SRV_max]);
ylim(ax1, [tau_all_min * 0.3, tau_all_max * 5]);
grid(ax1, 'on');
hold(ax1, 'off');

%% ----------------------------------------------------------------
%  RIGHT PANEL: SRV (y) vs tau (x)
%% ----------------------------------------------------------------
ax2 = subplot(1, 2, 2);
hold(ax2, 'on');

% Shaded regimes as horizontal bands
for ir = 1:numel(regimes)
    patch(ax2, ...
          [tau_all_min*0.1 tau_all_max*10 tau_all_max*10 tau_all_min*0.1], ...
          [regimes{ir}.lo regimes{ir}.lo regimes{ir}.hi regimes{ir}.hi], ...
          regimes{ir}.color, 'EdgeColor', 'none', 'FaceAlpha', 0.55, ...
          'HandleVisibility', 'off');
    text(ax2, tau_all_max * 2.5, sqrt(regimes{ir}.lo * regimes{ir}.hi), ...
         regimes{ir}.label, 'FontSize', 10, 'HorizontalAlignment', 'right', ...
         'Color', [0.25 0.25 0.25], 'FontWeight', 'bold');
end

% Emitter curves (axes swapped)
for k = 1:n_jd
    loglog(ax2, tau_emitter(k,:), SRV_vec, '-', ...
           'LineWidth', 1.8, 'Color', jd_colors(k,:), ...
           'HandleVisibility', 'off');
end

% Wafer curve (axes swapped)
loglog(ax2, tau_wafer, SRV_vec, 'k-', ...
       'LineWidth', 2.5, 'Color', wafer_color, ...
       'HandleVisibility', 'off');

xlabel(ax2, 'Minority carrier lifetime \tau_p (\mus)', 'FontSize', 15);
ylabel(ax2, 'Surface recombination velocity s_p (cm/s)', 'FontSize', 15);
title(ax2, 's_p vs \tau_p', 'FontSize', 15);
set(ax2, 'XScale', 'log', 'YScale', 'log', ...
         'FontSize', 13, 'LineWidth', 1.2, 'Box', 'on');
xlim(ax2, [tau_all_min * 0.3, tau_all_max * 5]);
ylim(ax2, [SRV_min, SRV_max]);
grid(ax2, 'on');
hold(ax2, 'off');

%% ----------------------------------------------------------------
%  Shared legend
%% ----------------------------------------------------------------
jd_labels = arrayfun(@(x) sprintf('Emitter: j_d = %d nm  [\\tau = j_d/s_p]', x), ...
                     jd_vals_nm, 'UniformOutput', false);
all_h      = [h_emitter; h_wafer];
all_labels = [jd_labels, {sprintf('Wafer: W = %d \\mum  [\\tau = W/(2s_p)]', W_wafer_um)}];

lg = legend(ax1, all_h, all_labels, 'Location', 'southwest', 'FontSize', 11);

sgtitle('SRV \leftrightarrow Minority Carrier Lifetime', 'FontSize', 16);

end