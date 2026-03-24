%% SRV from lifetime (planar wafer)
% Formulas:
%   1/τ_eff = 1/τ_bulk + (S_front + S_back)/W
%   → S_total = S_front + S_back = W * (1/τ_eff - 1/τ_bulk)
%   For identical surfaces: S_front = S_back = S_total/2
%
% Units: W in cm, τ in s, S in cm/s

clear; clc;

%% ---- USER INPUTS -------------------------------------------------------
W_um         = 280;                        % wafer thickness [µm]
tau_eff_ms   = [11.03*1e-3];                % effective lifetimes [ms] (vector OK)
tau_bulk_ms  = 18.283;                         % bulk lifetime [ms]; set [] if unknown
S_back_cmps  = [];                         % known back SRV [cm/s]; [] → assume symmetric

%% ---- UNIT CONVERSIONS --------------------------------------------------
W_cm        = W_um * 1e-4;                 % µm → cm
tau_eff_s   = tau_eff_ms * 1e-3;           % ms → s
if ~isempty(tau_bulk_ms)
    tau_bulk_s = tau_bulk_ms * 1e-3;       % ms → s
else
    tau_bulk_s = Inf;                      % unknown bulk → treat as very large
end

%% ---- COMPUTE SRV -------------------------------------------------------
delta = 1./tau_eff_s - 1./tau_bulk_s;      % s^-1
% Guard against small negative round-off:
delta(delta < 0 & delta > -1e-12) = 0;

S_total_cmps = W_cm .* delta;              % S_front + S_back  [cm/s]

if isempty(S_back_cmps)
    % Symmetric assumption
    S_front_cmps = S_total_cmps/2;
    S_back_used  = S_front_cmps;
else
    S_front_cmps = S_total_cmps - S_back_cmps;
    S_back_used  = S_back_cmps + zeros(size(S_front_cmps)); % vectorize for table
end

%% ---- PRINT RESULTS -----------------------------------------------------
T = table( tau_eff_ms(:), repmat(W_um,numel(tau_eff_ms),1), ...
           S_total_cmps(:), S_front_cmps(:), S_back_used(:), ...
           'VariableNames', {'tau_eff_ms','W_um','S_total_cmps','S_front_cmps','S_back_cmps'} );

disp('SRV results (cm/s):');
disp(T);

if isinf(tau_bulk_s)
    disp('Note: τ_bulk unknown → used upper-bound estimate S_total = W / τ_eff (symmetric: S = W/(2 τ_eff)).');
end

%% ---- OPTIONAL: quick plot (SRV vs τ_eff) ------------------------------
figure('Color','w'); 
plot(tau_eff_ms, S_front_cmps, 'o-','LineWidth',1.8,'MarkerSize',7); grid on; box on;
xlabel('Effective lifetime \tau_{eff} (ms)');
ylabel('Front SRV S_{front} (cm s^{-1})');
title('Surface Recombination Velocity vs. Effective Lifetime');