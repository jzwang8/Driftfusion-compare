%% Combined script: CT mobility and minority lifetime vs doping
%% Common doping range [cm^-3]
N = logspace(14, 20, 300);    % from 1e14 to 1e20 cm^-3

%% 1) Caughey–Thomas mobility vs. doping (log–log)

% Caughey–Thomas model parameters for Si (example values)
% Electrons
mu0_n    = 1350;      % low-doping mobility [cm^2/Vs]
mu_min_n = 65;        % high-doping mobility floor [cm^2/Vs]
Nref_n   = 1e17;      % reference doping [cm^-3]
alpha_n  = 0.72;      % fitting exponent

% Holes
mu0_p    = 480;       % low-doping mobility [cm^2/Vs]
mu_min_p = 50;        % high-doping mobility floor [cm^2/Vs]
Nref_p   = 2e17;      % reference doping [cm^-3]
alpha_p  = 0.70;      % fitting exponent

% Caughey–Thomas mobility function
ct_mu = @(mu0, mu_min, Nref, alpha, N) ...
    mu_min + (mu0 - mu_min) ./ (1 + (N./Nref).^alpha);

mu_n = ct_mu(mu0_n, mu_min_n, Nref_n, alpha_n, N);
mu_p = ct_mu(mu0_p, mu_min_p, Nref_p, alpha_p, N);

% Plot CT mobility
fig1 = figure(1);
ax1  = axes('Parent', fig1);

loglog(ax1, N, mu_n, 'LineWidth', 2); hold(ax1, 'on');
loglog(ax1, N, mu_p, 'LineWidth', 2);

xlabel(ax1, 'Doping concentration N [cm^{-3}]');
ylabel(ax1, 'Mobility \mu [cm^2/Vs]');
title(ax1, 'Caughey–Thomas mobility vs. doping');
legend(ax1, {'Electrons','Holes'}, 'Location', 'northeast');
grid(ax1, 'on');
box(ax1, 'on');

%% 2) Minority carrier lifetime vs. doping (SRH + Auger, log–log)

% Lifetime model parameters
tau_SRH = 5e-3;    % [s] SRH lifetime (good FZ Si)

% Auger coefficients [cm^6/s]
C_n = 1e-31;       % minority electrons in p-type
C_p = 2e-31;       % minority holes in n-type

% Effective minority lifetimes (SRH || Auger)
tau_n_p = 1 ./ (1/tau_SRH + C_n .* N.^2);  % n in p-Si
tau_p_n = 1 ./ (1/tau_SRH + C_p .* N.^2);  % p in n-Si

% Plot lifetimes
fig2 = figure(2);
ax2  = axes('Parent', fig2);

loglog(ax2, N, tau_n_p, 'LineWidth', 2); hold(ax2, 'on');
loglog(ax2, N, tau_p_n, 'LineWidth', 2);

xlabel(ax2, 'Doping concentration N [cm^{-3}]');
ylabel(ax2, 'Minority carrier lifetime \tau [s]');
title(ax2, 'Minority carrier lifetime vs. doping');
legend(ax2, {'electrons in p-Si', 'holes in n-Si'}, 'Location', 'southwest');
grid(ax2, 'on');
box(ax2, 'on');