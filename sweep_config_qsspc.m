function config = sweep_config_qsspc()
%SWEEP_CONFIG_QSSPC Single source of truth for QSSPC lifetime sweep parameters
%
%   The CSV file is a TEMPLATE. Sweep parameters override template
%   values per-run in run_sweep_chunk_qsspc.m.
%
%   Swept: sp_r, Vi (Qf via Phi_right offset), ND (wafer doping)
%   For each ND, the per-run override recalculates EF0, mu_n, mu_p, taup.

%% ==================== TEMPLATE CSV ====================
config.file_name = 'Input_files/ntype.csv';

%% ==================== READ TEMPLATE VALUES FROM CSV ====================
csv_data = readtable(config.file_name, 'Delimiter', ',');

active_idx = find(strcmp(csv_data.layer_type, 'active'));
elec_idx   = find(strcmp(csv_data.layer_type, 'electrode'));
left_elec  = elec_idx(1);
right_elec = elec_idx(end);

config.template.thickness = csv_data.thickness(active_idx);
config.template.Phi_EA    = csv_data.Phi_EA(active_idx);
config.template.Phi_IP    = csv_data.Phi_IP(active_idx);
config.template.EF0       = csv_data.EF0(active_idx);
config.template.Et        = csv_data.Et(active_idx);
config.template.Nc        = csv_data.Nc(active_idx);
config.template.Nv        = csv_data.Nv(active_idx);
config.template.mu_n      = csv_data.mu_n(active_idx);
config.template.mu_p      = csv_data.mu_p(active_idx);
config.template.taun      = csv_data.taun(active_idx);
config.template.taup      = csv_data.taup(active_idx);
config.template.epp       = csv_data.epp(active_idx);
config.template.B         = csv_data.B(active_idx);
config.template.sn_l      = csv_data.sn(left_elec);
config.template.sp_l      = csv_data.sp(left_elec);
config.template.sn_r      = csv_data.sn(right_elec);
config.template.sp_r      = csv_data.sp(right_elec);

config.ni_Si = 1.04e10;

%% ==================== SWEEP RANGES ====================

% ND: wafer doping [cm^-3] (from resistivity)
%   For each ND, the per-run override recalculates:
%     EF0  = getEfFromDopingConc(ND)
%     [mu_n, mu_p] = getMobilitiesFromDopingConc(ND)
%     taup = getTaupFromDopingConc(ND)
%     Phi_left  = EF0
%     Phi_right = EF0 + Vi
% config.resistivity_values = [3];                                % [ohm.cm] -- add more later
% config.ND_values = getDopingConcFromResistivity(config.resistivity_values);    % [cm^-3]
config.ND_values = [1e17];    % [cm^-3], 1.65e15 is 3ohm.cm

% sp_r: hole SRV at Tc/ZnPc surface [cm/s]
config.sp_r_values = [10, 100, 1000, 10000];

% Vi: electrode work function offset [V]
% config.Vi_values_V = [-0.05, -0.02, -0.01, 0, 0.01, 0.02, 0.05, 0.1, 0.2];
config.Vi_values_V = [0];

% Corresponding Qf from Vi (for labeling only)
config.Qf_values_cm2 = calcQfFromVi(config.Vi_values_V);

%% ==================== SIMULATION SETTINGS ====================
config.light_intensity = 1;      % [Suns]
config.Rs_ohm_cm2     = 1e6;    % [ohm.cm^2] high Rs for open circuit
config.t_lighton_s    = 1e-1;   % [s] time to reach steady state
config.tpoints        = 200;

%% ==================== INFRASTRUCTURE ====================
config.outputDir      = 'qsspc_sweep_results';
config.runs_per_chunk = 20;

%% ==================== PRINT SUMMARY ====================
n_ND = length(config.ND_values);
n_sp = length(config.sp_r_values);
n_Vi = length(config.Vi_values_V);
n_total = n_ND * n_sp * n_Vi;

fprintf('=== QSSPC Sweep Config ===\n');
fprintf('  Template CSV: %s\n', config.file_name);
fprintf('  Sweep: %d ND x %d sp_r x %d Vi = %d combinations\n', n_ND, n_sp, n_Vi, n_total);
for ii = 1:n_ND
    EF0_ii = getEfFromDopingConc(config.ND_values(ii));
    fprintf('  ND=%.2e  -> EF0=%.4f eV\n', ...
        config.ND_values(ii), EF0_ii);
end
fprintf('  sp_r: [%s] cm/s\n', num2str(config.sp_r_values, '%.0f '));
fprintf('  Vi:   [%s] V\n', num2str(config.Vi_values_V, '%.3f '));
fprintf('  Runs per chunk: %d\n', config.runs_per_chunk);
fprintf('==========================\n');

end