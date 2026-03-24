function compare_sweep_profiles()
%COMPARE_SWEEP_PROFILES
% Load multiple Driftfusion sweep .mat files and create publication-quality
% overlay plots for carrier density, charge density, electric field,
% potential, band energies, and recombination rates.
%
% Designed for research poster/presentation aesthetics.
%
% Usage:
%   1. Set the sweep_name and sweep_vals to define which variable to compare
%   2. Set fixed parameters for all other variables
%   3. Run the script
%
% Files are loaded from sweep_results/ using the naming convention:
%   SF{0/1}_eta{X}_jd{X}_sp{X}_ND{X}_PhiOff{X}.mat

close all;

%% ========================================================================
%  USER CONFIGURATION
%% ========================================================================

results_dir = "sweep_results/sweep_results_poster";   % folder containing .mat files

% ---- Sweep configuration ----
% Choose which variable to sweep for comparison:
%   "junction_depth" | "sp_r" | "eta_sf" | "EF" | "Phi_R_offset"
sweep_name = "Phi_R_offset";

% Values to overlay (typically 3-5 runs for clarity)
% sweep_vals = [-.0539,0,.0539];   % example: junction depth (nm)
sweep_vals = [-.0539,0,.0539];

% Special config for eta_sf sweep (SF_flag and eta pairs)
% Used only when sweep_name = "eta_sf"
sweep_vals_SF  = [0, 1, 1, 1];   % SF_flag for each run
sweep_vals_eta = [0, 0, 1, 2];   % eta_sf for each run

% ---- Fixed parameters (held constant while sweeping one variable) ----
SF_flag_fixed        = 1;
eta_sf_fixed         = 2.0;      % if SF_flag_fixed==0, eta is forced to 0
junction_depth_fixed = 300;      % nm (ignored if sweeping junction_depth)
sp_r_fixed           = 1e4;      % cm/s (ignored if sweeping sp_r)
EF_fixed             = -4.104;   % eV (ignored if sweeping EF)
Phi_R_offset_fixed   = 0;  % eV (ignored if sweeping Phi_R_offset)

% ---- Plot settings ----
tarr_mode = "end";               % "end" or numeric time (s)
xrange    = [1.795e5, 1.804e5];  % x-axis range in nm (emitter/SF region)
% xrange    = [1.802e5, 1.8032e5];

% Which plots to generate:
do_carrier_density = true;       % n, p (log y-axis)
do_charge_density  = true;       % rho
do_electric_field  = true;       % F
do_potential       = true;       % V
do_band_energies   = true;       % Ecb, Evb, Efn, Efp
do_recombination   = true;       % btb, SRH, vsr

% Figure size [width, height] in pixels
fig_size = [600, 400];

% Physical constraint for validation
E_C = -4.05;  % Si conduction band edge (eV)

%% ========================================================================
%  COLOR SCHEME - Colorblind-friendly (blue / grey / orange)
%% ========================================================================

% Qualitative color palette — matches poster figures script
% Designed for deuteranopia/protanopia accessibility
colors = [
    0.12, 0.47, 0.71;   % blue
    0.50, 0.50, 0.50;   % grey
    0.85, 0.37, 0.01;   % orange
    0.55, 0.34, 0.64;   % purple
    0.00, 0.60, 0.50;   % teal
    0.90, 0.60, 0.00;   % amber
    0.64, 0.08, 0.18;   % crimson
    0.50, 0.70, 0.86;   % light blue
    0.93, 0.60, 0.25;   % light orange
    0.40, 0.76, 0.65;   % mint
    0.60, 0.60, 0.00;   % olive
    0.80, 0.40, 0.70;   % pink-purple
];

%% ========================================================================
%  BUILD RUN LIST
%% ========================================================================

runs = struct( ...
    'SF_flag', {}, 'eta_sf', {}, 'junction_depth', {}, 'sp_r', {}, ...
    'EF_emitter', {}, 'N_D', {}, 'N_D_token', {}, 'Phi_R_offset', {}, 'Phi_R', {}, ...
    'filename', {}, 'label', {}, 'color', {} );

% Determine number of runs based on sweep type
if sweep_name == "eta_sf"
    n_runs = numel(sweep_vals_SF);
else
    n_runs = numel(sweep_vals);
end

for k = 1:n_runs
    p = struct();
    
    p.SF_flag = SF_flag_fixed;
    
    % eta handling
    if p.SF_flag == 0
        p.eta_sf = 0;
    else
        p.eta_sf = eta_sf_fixed;
    end
    
    % start with fixed values
    p.junction_depth = junction_depth_fixed;
    p.sp_r           = sp_r_fixed;
    p.EF_emitter     = EF_fixed;
    p.Phi_R_offset   = Phi_R_offset_fixed;
    
    % override swept variable
    switch sweep_name
        case "junction_depth"
            p.junction_depth = sweep_vals(k);
        case "sp_r"
            p.sp_r = sweep_vals(k);
        case "eta_sf"
            % Special handling: sweep_vals encodes [SF_flag, eta_sf] pairs
            % Use sweep_vals_eta and sweep_vals_SF defined in config
            p.SF_flag = sweep_vals_SF(k);
            p.eta_sf = sweep_vals_eta(k);
        case "EF"
            p.EF_emitter = sweep_vals(k);
        case "Phi_R_offset"
            p.Phi_R_offset = sweep_vals(k);
        otherwise
            error("Unknown sweep_name: %s", sweep_name);
    end
    
    % compute N_D from EF
    [N_D, ~] = getDopingConcFromEf(p.EF_emitter);
    p.N_D = N_D;
    
    % quantize ND to match file token convention ND%.0e
    p.N_D_token = str2double(sprintf('%.0e', p.N_D));
    
    % compute Phi_R and validate
    p.Phi_R = p.EF_emitter + p.Phi_R_offset;
    if p.Phi_R > E_C
        fprintf("Skipping unphysical combo: Phi_R=%.4f > E_C=%.4f\n", p.Phi_R, E_C);
        continue;
    end
    
    % filename from generator
    p.filename = generate_filename(p.SF_flag, p.eta_sf, p.junction_depth, ...
                                   p.sp_r, p.N_D_token, p.Phi_R_offset);
    
    % legend label
    p.label = make_label(p, sweep_name);
    
    % assign color (cycle if more runs than colors)
    color_idx = mod(k-1, size(colors, 1)) + 1;
    p.color = colors(color_idx, :);
    
    runs(end+1) = p; %#ok<AGROW>
end

if isempty(runs)
    error("No valid runs constructed. Check sweep/fixed parameters.");
end

fprintf('\n=== Comparing %d runs ===\n', numel(runs));
for i = 1:numel(runs)
    fprintf('  %d. %s\n', i, runs(i).label);
end
fprintf('\n');

%% ========================================================================
%  LOAD ALL RUNS
%% ========================================================================

data = struct([]);

for i = 1:numel(runs)
    fpath = fullfile(results_dir, runs(i).filename);
    
    if ~isfile(fpath)
        warning("Missing file: %s", fpath);
        continue;
    end
    
    S = load(fpath);
    
    if ~isfield(S, "thisRun") || ~isfield(S.thisRun, "solEq") || ~isfield(S.thisRun.solEq, "el")
        warning("File %s missing expected variable thisRun.solEq.el", fpath);
        continue;
    end
    
    sol = S.thisRun.solEq.el;
    
    fprintf("Loaded: %s\n", runs(i).filename);
    
    % Extract data using dfana
    [~, ~, x, par, ~, n, p_carr, ~, ~, V] = dfana.splitsol(sol);
    
    rho = dfana.calcrho(sol, "whole");
    F   = dfana.calcF(sol, "whole");
    
    [Ecb, Evb, Efn, Efp] = dfana.calcEnergies(sol);
    
    rsub  = dfana.calcr(sol, "sub");
    x_sub = par.x_sub;
    
    if tarr_mode == "end"
        tarr = sol.t(end);
    else
        tarr = tarr_mode;
    end
    
    % Find time index
    t_idx = find(sol.t <= tarr, 1, 'last');
    if isempty(t_idx)
        t_idx = 1;
    end
    
    % Store in data struct
    idx = numel(data) + 1;
    data(idx).label = runs(i).label;
    data(idx).color = runs(i).color;
    data(idx).x     = x * 1e7;      % cm -> nm
    data(idx).x_sub = x_sub * 1e7;  % cm -> nm
    data(idx).t_idx = t_idx;
    
    % Extract profiles at time index
    data(idx).n   = n(t_idx, :);
    data(idx).p   = p_carr(t_idx, :);
    data(idx).rho = rho(t_idx, :);
    data(idx).F   = F(t_idx, :);
    data(idx).V   = V(t_idx, :);
    
    data(idx).Efn = Efn(t_idx, :);
    data(idx).Efp = Efp(t_idx, :);
    data(idx).Ecb = Ecb(t_idx, :);
    data(idx).Evb = Evb(t_idx, :);
    
    data(idx).r_btb = rsub.btb(t_idx, :);
    data(idx).r_srh = rsub.srh(t_idx, :);
    data(idx).r_vsr = rsub.vsr(t_idx, :);
    
    data(idx).sol = sol;  % keep for colourblocks
end

if isempty(data)
    error("All selected files were missing or invalid.");
end

fprintf('\nSuccessfully loaded %d runs\n\n', numel(data));

%% ========================================================================
%  GENERATE PLOTS
%% ========================================================================

% Common axis/font settings
ax_fontsize    = 20;
label_fontsize = 20;
legend_fontsize = 16;
line_width     = 3;
line_alpha     = 0.7;  % transparency (0 = invisible, 1 = opaque)

% Get sweep axis label for titles
sweep_label = get_sweep_axis_label(sweep_name);

%% --- Figure 1: Carrier Density (n, p) - LOG Y-AXIS ---
if do_carrier_density
    fig1 = figure('Name', 'Carrier Density', ...
                  'Position', [50, 50, fig_size(1), fig_size(2)], ...
                  'Color', 'w');
    ax1 = axes(fig1);
    hold(ax1, 'on');
    
    h_lines = gobjects(numel(data), 2);  % handles for legend
    
    for i = 1:numel(data)
        h_lines(i,1) = semilogy(ax1, data(i).x, data(i).n, '-', ...
                                'Color', [data(i).color, line_alpha], ...
                                'LineWidth', line_width);
        h_lines(i,2) = semilogy(ax1, data(i).x, data(i).p, ':', ...
                                'Color', [data(i).color, line_alpha], ...
                                'LineWidth', line_width);
    end
    
    % Build legend: show run labels + one entry explaining line styles
    leg_handles = [h_lines(:,1); ...
                   plot(ax1, NaN, NaN, 'k-', 'LineWidth', line_width); ...
                   plot(ax1, NaN, NaN, 'k:', 'LineWidth', line_width)];
    leg_entries = [{data.label}, {'n (solid)'}, {'p (dotted)'}];
    
    xlabel(ax1, 'Position (nm)', 'FontSize', label_fontsize);
    ylabel(ax1, 'Carrier Density (cm^{-3})', 'FontSize', label_fontsize);
    xlim(ax1, xrange);
    % ylim([1e18,1e19])
    set(ax1, 'YScale', 'log', 'FontSize', ax_fontsize, ...
             'LineWidth', 1.2, 'Box', 'on');
    grid(ax1, 'on');
    legend(ax1, leg_handles, leg_entries, 'Location', 'best', ...
           'FontSize', legend_fontsize);
    
    hold(ax1, 'off');
end

%% --- Figure 2: Charge Density (rho) ---
if do_charge_density
    fig2 = figure('Name', 'Charge Density', ...
                  'Position', [100, 100, fig_size(1), fig_size(2)], ...
                  'Color', 'w');
    ax2 = axes(fig2);
    hold(ax2, 'on');
    
    h_lines = gobjects(numel(data), 1);
    
    for i = 1:numel(data)
        h_lines(i) = plot(ax2, data(i).x, data(i).rho, '-', ...
                          'Color', [data(i).color, line_alpha], ...
                          'LineWidth', line_width);
    end
    
    xlabel(ax2, 'Position (nm)', 'FontSize', label_fontsize);
    ylabel(ax2, 'Charge Density (cm^{-3})', 'FontSize', label_fontsize);
    xlim(ax2, xrange);
    set(ax2, 'FontSize', ax_fontsize, 'LineWidth', 1.2, 'Box', 'on');
    grid(ax2, 'on');
    legend(ax2, h_lines, {data.label}, 'Location', 'best', ...
           'FontSize', legend_fontsize);
    
    hold(ax2, 'off');
end

%% --- Figure 3: Electric Field (F) ---
if do_electric_field
    fig3 = figure('Name', 'Electric Field', ...
                  'Position', [150, 150, fig_size(1), fig_size(2)], ...
                  'Color', 'w');
    ax3 = axes(fig3);
    hold(ax3, 'on');
    
    h_lines = gobjects(numel(data), 1);
    
    for i = 1:numel(data)
        h_lines(i) = plot(ax3, data(i).x, data(i).F, '-', ...
                          'Color', [data(i).color, line_alpha], ...
                          'LineWidth', line_width);
    end
    
    xlabel(ax3, 'Position (nm)', 'FontSize', label_fontsize);
    ylabel(ax3, 'Electric Field (V cm^{-1})', 'FontSize', label_fontsize);
    xlim(ax3, xrange);
    set(ax3, 'FontSize', ax_fontsize, 'LineWidth', 1.2, 'Box', 'on');
    grid(ax3, 'on');
    legend(ax3, h_lines, {data.label}, 'Location', 'best', ...
           'FontSize', legend_fontsize);
    
    hold(ax3, 'off');
end

%% --- Figure 4: Electrostatic Potential (V) ---
if do_potential
    fig4 = figure('Name', 'Electrostatic Potential', ...
                  'Position', [200, 200, fig_size(1), fig_size(2)], ...
                  'Color', 'w');
    ax4 = axes(fig4);
    hold(ax4, 'on');
    
    h_lines = gobjects(numel(data), 1);
    
    for i = 1:numel(data)
        h_lines(i) = plot(ax4, data(i).x, data(i).V, '-', ...
                          'Color', [data(i).color, line_alpha], ...
                          'LineWidth', line_width);
    end
    
    xlabel(ax4, 'Position (nm)', 'FontSize', label_fontsize);
    ylabel(ax4, 'Electrostatic Potential (V)', 'FontSize', label_fontsize);
    xlim(ax4, xrange);
    set(ax4, 'FontSize', ax_fontsize, 'LineWidth', 1.2, 'Box', 'on');
    grid(ax4, 'on');
    legend(ax4, h_lines, {data.label}, 'Location', 'best', ...
           'FontSize', legend_fontsize);
    
    hold(ax4, 'off');
end

%% --- Figure 5: Band Energies (Ecb, Evb, Ef) ---
if do_band_energies
    fig5 = figure('Name', 'Band Energies', ...
                  'Position', [250, 250, fig_size(1), fig_size(2)], ...
                  'Color', 'w');
    ax5 = axes(fig5);
    hold(ax5, 'on');
    
    % For band diagram: solid for bands, dotted for quasi-Fermi level
    % Each run gets its own color; line style distinguishes band type
    h_lines = gobjects(numel(data), 3);  % Ecb, Evb, Ef
    
    for i = 1:numel(data)
        h_lines(i,1) = plot(ax5, data(i).x, data(i).Ecb, '-', ...
                            'Color', [data(i).color, line_alpha], ...
                            'LineWidth', line_width);
        h_lines(i,2) = plot(ax5, data(i).x, data(i).Evb, '-', ...
                            'Color', [data(i).color, line_alpha], ...
                            'LineWidth', line_width * 0.7);
        h_lines(i,3) = plot(ax5, data(i).x, data(i).Efn, ':', ...
                            'Color', [data(i).color, line_alpha], ...
                            'LineWidth', line_width);
    end
    
    % Build legend: show run labels + one entry explaining line styles
    leg_handles = [h_lines(:,1); ...
                   plot(ax5, NaN, NaN, 'k-', 'LineWidth', line_width); ...
                   plot(ax5, NaN, NaN, 'k:', 'LineWidth', line_width)];
    leg_entries = [{data.label}, {'E_{CB}, E_{VB} (solid)'}, ...
                   {'E_F (dotted)'}];
    
    xlabel(ax5, 'Position (nm)', 'FontSize', label_fontsize);
    ylabel(ax5, 'Energy (eV)', 'FontSize', label_fontsize);
    xlim(ax5, xrange);
    set(ax5, 'FontSize', ax_fontsize, 'LineWidth', 1.2, 'Box', 'on');
    grid(ax5, 'on');
    legend(ax5, leg_handles, leg_entries, 'Location', 'best', ...
           'FontSize', legend_fontsize);
    
    hold(ax5, 'off');
end

%% --- Figure 6: Recombination Rates (btb, SRH, vsr) ---
if do_recombination
    fig6 = figure('Name', 'Recombination Rates', ...
                  'Position', [300, 300, fig_size(1), fig_size(2)], ...
                  'Color', 'w');
    ax6 = axes(fig6);
    hold(ax6, 'on');
    
    h_lines = gobjects(numel(data), 3);  % btb, srh, vsr
    
    for i = 1:numel(data)
        h_lines(i,1) = semilogy(ax6, data(i).x_sub, abs(data(i).r_btb), '-', ...
                                'Color', [data(i).color, line_alpha], ...
                                'LineWidth', line_width);
        h_lines(i,2) = semilogy(ax6, data(i).x_sub, abs(data(i).r_srh), ':', ...
                                'Color', [data(i).color, line_alpha], ...
                                'LineWidth', line_width);
        h_lines(i,3) = semilogy(ax6, data(i).x_sub, abs(data(i).r_vsr), '--', ...
                                'Color', [data(i).color, line_alpha], ...
                                'LineWidth', line_width);
    end
    
    % Build legend
    leg_handles = [h_lines(:,1); ...
                   plot(ax6, NaN, NaN, 'k-', 'LineWidth', line_width); ...
                   plot(ax6, NaN, NaN, 'k:', 'LineWidth', line_width); ...
                   plot(ax6, NaN, NaN, 'k--', 'LineWidth', line_width)];
    leg_entries = [{data.label}, {'r_{btb} (solid)'}, ...
                   {'r_{SRH} (dotted)'}, {'r_{vsr} (dashed)'}];
    
    xlabel(ax6, 'Position (nm)', 'FontSize', label_fontsize);
    ylabel(ax6, 'Recombination Rate (cm^{-3} s^{-1})', 'FontSize', label_fontsize);
    xlim(ax6, xrange);
    set(ax6, 'YScale', 'log', 'FontSize', ax_fontsize, ...
             'LineWidth', 1.2, 'Box', 'on');
    grid(ax6, 'on');
    legend(ax6, leg_handles, leg_entries, 'Location', 'best', ...
           'FontSize', legend_fontsize);
    
    hold(ax6, 'off');
end

fprintf('=== Generated %d figures ===\n', ...
        do_carrier_density + do_charge_density + do_electric_field + ...
        do_potential + do_band_energies + do_recombination);

end  % end main function


%% ========================================================================
%  HELPER FUNCTIONS
%% ========================================================================

function filename = generate_filename(SF_flag, SF_efficiency, junction_depth, sp_r, N_D, Phi_R_offset)
%GENERATE_FILENAME Create deterministic filename from parameters
%   Format: SF{0/1}_eta{X}_jd{X}_sp{X}_ND{X}_PhiOff{X}.mat

    sf_str  = sprintf('SF%d', SF_flag);
    eta_str = sprintf('eta%.1f', SF_efficiency);
    jd_str  = sprintf('jd%d', junction_depth);
    
    sp_str = sprintf('sp%.0e', sp_r);
    nd_str = sprintf('ND%.0e', N_D);
    
    % Normalize exponent formatting (e.g., 1e05 not 1e+05)
    sp_str = strrep(sp_str, '+', '');
    nd_str = strrep(nd_str, '+', '');
    
    % Ensure two-digit exponent (e.g., e5 -> e05)
    sp_str = regexprep(sp_str, 'e(-?)(\d)$', 'e$10$2');
    nd_str = regexprep(nd_str, 'e(-?)(\d)$', 'e$10$2');
    
    phi_str = sprintf('PhiOff%.4f', Phi_R_offset);
    phi_str = strrep(phi_str, '-', 'n');  % negative sign -> 'n'
    phi_str = strrep(phi_str, '.', 'p');  % decimal point -> 'p'
    
    filename = sprintf('%s_%s_%s_%s_%s_%s.mat', sf_str, eta_str, jd_str, sp_str, nd_str, phi_str);
end


function lbl = make_label(p, sweep_name)
%MAKE_LABEL Create legend label based on swept parameter

    switch sweep_name
        case "junction_depth"
            lbl = sprintf("j_d = %d nm", round(p.junction_depth));
        case "sp_r"
            lbl = sprintf("s_p = %.0e cm/s", p.sp_r);
        case "eta_sf"
            % Special labels for SF_flag and eta combinations
            if p.SF_flag == 0
                lbl = "No SF layer";
            else
                lbl = sprintf("\\eta_{SF} = %.0f", p.eta_sf);
            end
        case "EF"
            lbl = sprintf("E_F = %.3f eV", p.EF_emitter);
        case "Phi_R_offset"
            % Use Qf notation instead of Phi_R offset value
            if p.Phi_R_offset > 0
                lbl = "Q_f > 0";
            elseif p.Phi_R_offset < 0
                lbl = "Q_f < 0";
            else
                lbl = "Q_f = 0";
            end
        otherwise
            lbl = "run";
    end
    lbl = string(lbl);
end


function lbl = get_sweep_axis_label(sweep_name)
%GET_SWEEP_AXIS_LABEL Return axis label string for sweep variable

    switch sweep_name
        case "junction_depth"
            lbl = "Junction Depth (nm)";
        case "sp_r"
            lbl = "Surface Recombination Velocity (cm/s)";
        case "eta_sf"
            lbl = "\eta_{SF}";
        case "EF"
            lbl = "Fermi Level E_F (eV)";
        case "Phi_R_offset"
            lbl = "\Phi_R Offset (eV)";
        otherwise
            lbl = sweep_name;
    end
end