function [EF_profile, mu_n_profile, mu_p_profile, taup_profile, N_D_profile, x_sublayers] = ...
        generate_emitter_EF_profile(N0_peak, junction_depth_nm, profile_type, N_sublayers)
%GENERATE_EMITTER_EF_PROFILE Generate spatially varying emitter doping profile
%
%   Discretizes the emitter into N_sublayers of equal thickness, each with
%   a Fermi level, mobility, and lifetime derived from its local doping.
%
%   Inputs:
%     N0_peak           - peak/surface doping concentration [cm^-3]
%     junction_depth_nm - total emitter thickness [nm]
%     profile_type      - 'uniform', 'gaussian', or 'exponential'
%     N_sublayers       - number of sublayers (e.g. 50)
%
%   Outputs:
%     EF_profile        - [1 x N_sublayers] Fermi levels [eV, vacuum ref]
%     mu_n_profile      - [1 x N_sublayers] electron mobilities [cm^2/Vs]
%     mu_p_profile      - [1 x N_sublayers] hole mobilities [cm^2/Vs]
%     taup_profile      - [1 x N_sublayers] hole lifetimes [s]
%     N_D_profile       - [1 x N_sublayers] doping concentrations [cm^-3]
%     x_sublayers       - [1 x N_sublayers] sublayer center depths [nm]
%
%   Profile definitions:
%     Gaussian:    N(x) = N0 * exp(-(x - x0)^2 / (2*sigma^2))
%                  x0 = d/2, sigma = d/5  (peak at center, covers emitter)
%     Exponential: N(x) = A * exp(-x/Ld)
%                  A, Ld derived by dose-matching to the Gaussian profile
%                  (same total integrated dopant dose)
%     Uniform:     N(x) = N0_peak everywhere

    d = junction_depth_nm;  % nm

    % Sublayer center positions (depth from surface, nm)
    sublayer_thickness = d / N_sublayers;
    x_sublayers = sublayer_thickness/2 : sublayer_thickness : d - sublayer_thickness/2;

    %% Compute N_D profile
    switch lower(profile_type)

        case 'uniform'
            N_D_profile = N0_peak * ones(1, N_sublayers);

        case 'gaussian'
            x0    = d / 2;      % peak at center of emitter
            sigma = d / 5;      % width ~ covers emitter well
            N_D_profile = N0_peak * exp(-((x_sublayers - x0).^2) / (2 * sigma^2));

            % Clamp minimum to background doping to avoid unphysical values
            N_D_profile = max(N_D_profile, 1e15);

        case 'exponential'
            % Derive Ld and A by dose-matching to the equivalent Gaussian
            % so both profiles have the same peak N0, total dose, junction depth
            x0    = d / 2;
            sigma = d / 5;

            % Full Gaussian profile on fine grid for integration
            x_fine = linspace(0, d, 1000);
            N_gauss = N0_peak * exp(-((x_fine - x0).^2) / (2 * sigma^2));
            total_dose_gauss = trapz(x_fine, N_gauss);  % nm * cm^-3

            % Dose-match: find Ld such that
            % integral_0^d [A * exp(-x/Ld)] dx = total_dose_gauss
            % with A = N0_peak (same surface concentration)
            % integral = A * Ld * (1 - exp(-d/Ld))
            % => Ld * (1 - exp(-d/Ld)) = total_dose_gauss / N0_peak
            target_ratio = total_dose_gauss / N0_peak;  % nm

            dose_residual = @(Ld) Ld * (1 - exp(-d/Ld)) - target_ratio;

            % fzero with bracket [1, d] nm
            try
                Ld = fzero(dose_residual, [0.1, d*10]);
            catch
                warning('generate_emitter_EF_profile: fzero failed, using Ld = d/3');
                Ld = d / 3;
            end

            A = N0_peak;  % surface concentration = peak

            % Check if A would push EF above conduction band (degenerate)
            E_C = -4.05;  % eV
            EF_surface = getEfFromDopingConc(A);
            if EF_surface > E_C
                % Cap at conduction band edge and re-solve for Ld
                A = getDopingConcFromEf(E_C);
                fprintf('  [profile] Surface doping capped at E_C: A=%.2e cm^-3\n', A);
                % Re-solve Ld with capped A
                dose_residual2 = @(Ld) A * Ld * (1 - exp(-d/Ld)) - total_dose_gauss;
                try
                    Ld = fzero(dose_residual2, [0.1, d*10]);
                catch
                    Ld = d / 3;
                end
            end

            N_D_profile = A * exp(-x_sublayers / Ld);

            % Clamp minimum
            N_D_profile = max(N_D_profile, 1e15);

            fprintf('  [profile] Exponential: A=%.2e cm^-3, Ld=%.1f nm\n', A, Ld);
            fprintf('  [profile] Dose check: Gauss=%.3e, Exp=%.3e nm*cm^-3\n', ...
                    total_dose_gauss, trapz(x_sublayers, N_D_profile));

        otherwise
            error('generate_emitter_EF_profile: unknown profile_type "%s". Use "uniform", "gaussian", or "exponential".', profile_type);
    end

    %% Convert N_D -> EF, mobility, lifetime for each sublayer
    EF_profile   = zeros(1, N_sublayers);
    mu_n_profile = zeros(1, N_sublayers);
    mu_p_profile = zeros(1, N_sublayers);
    taup_profile = zeros(1, N_sublayers);

    for i = 1:N_sublayers
        N_D_i = N_D_profile(i);

        % Fermi level
        EF_profile(i) = getEfFromDopingConc(N_D_i);

        % Mobilities (doping-dependent)
        [mu_n_i, mu_p_i] = getMobilitiesFromDopingConc(N_D_i);
        mu_n_profile(i) = mu_n_i;
        mu_p_profile(i) = mu_p_i;

        % Minority carrier (hole) lifetime
        taup_profile(i) = getTaupFromDopingConc(N_D_i);
    end

end