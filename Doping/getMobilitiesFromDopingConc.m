function [mu_n, mu_p] = getMobilitiesFromDopingConc(N)
%GETMOBILITIESFROMDOPINGCONC  Doping-dependent low-field mobilities in Si at 300K.
%
%   [mu_n, mu_p] = getMobilitiesFromDopingConc(N)
%
%   Input:
%       N  - total ionized impurity concentration (cm^-3), scalar or vector
%
%   Output:
%       mu_n - electron mobility (cm^2/V·s)
%       mu_p - hole mobility (cm^2/V·s)
%
%   Model: Caughey-Thomas empirical fit for Si at 300K.
%          mu(N) = mu_min + (mu_max - mu_min) / (1 + (N/N_ref)^alpha)
%
%   Parameters from Sze, Physics of Semiconductor Devices (also Wikipedia).
%   Original model: D.M. Caughey and R.E. Thomas, Proc. IEEE 55, 2192 (1967).

    % Ensure N is positive
    N = abs(N);

    % ---- Electrons ----
    mumin_n  = 65;       % cm^2/V·s  (high-doping floor)
    mumax_n  = 1330;     % cm^2/V·s  (low-doping plateau, = mumin_n + 1265)
    Nref_n   = 8.5e16;   % cm^-3     (rolloff concentration)
    alpha_n  = 0.72;     % unitless  (transition steepness)

    mu_n = mumin_n + (mumax_n - mumin_n) ./ (1 + (N ./ Nref_n).^alpha_n);

    % ---- Holes ----
    mumin_p  = 48;       % cm^2/V·s  (high-doping floor)
    mumax_p  = 495;      % cm^2/V·s  (low-doping plateau, = mumin_p + 447)
    Nref_p   = 6.3e16;   % cm^-3     (rolloff concentration)
    alpha_p  = 0.76;     % unitless  (transition steepness)

    mu_p = mumin_p + (mumax_p - mumin_p) ./ (1 + (N ./ Nref_p).^alpha_p);

end

% OLD - updated 3/6/26
% function [mu_n, mu_p] = getMobilitiesFromDopingConc(N)
% %   Doping-dependent low-field mobilities in silicon at 300 K.
% %
% %   [mu_n, mu_p] = getMobilitiesFromDopingConc(N)
% %
% %   Input:
% %       N  - total doping concentration (cm^-3), can be scalar or vector
% %            N = NA + ND (absolute value, not signed)
% %
% %   Output:
% %       mu_n - electron mobility (cm^2/V·s)
% %       mu_p - hole mobility (cm^2/V·s)
% %
% %   Model: Caughey–Thomas-type empirical fit for Si at 300 K.
% %          mu(N) = mu_min + (mu_0 - mu_min) / (1 + (N/N_ref)^alpha)
% 
%     % Ensure N is positive (use magnitude of net doping)
%     N = abs(N);
% 
%     % ---- Electrons (n) ----
%     mu0_n   = 1417;       % cm^2/V·s  (low-doping electron mobility)
%     mumin_n = 52.2;       % cm^2/V·s  (high-doping limit)
%     Nref_n  = 9.68e16;    % cm^-3
%     alpha_n = 0.68;       % unitless
% 
%     mu_n = mumin_n + (mu0_n - mumin_n) ./ (1 + (N ./ Nref_n) .^ alpha_n);
% 
%     % ---- Holes (p) ----
%     mu0_p   = 470.5;      % cm^2/V·s  (low-doping hole mobility)
%     mumin_p = 44.9;       % cm^2/V·s  (high-doping limit)
%     Nref_p  = 2.35e17;    % cm^-3
%     alpha_p = 0.7;        % unitless
% 
%     mu_p = mumin_p + (mu0_p - mumin_p) ./ (1 + (N ./ Nref_p) .^ alpha_p);
% end