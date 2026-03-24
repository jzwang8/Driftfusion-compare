function ND = getDopingConcFromResistivity(rho_ohmcm)
%RHOTOND_NTYPE Convert resistivity to donor concentration for n-type Si
%
%   ND = rhoToND_ntype(rho_ohmcm)
%
%   Iteratively solves rho = 1/(q * mu_n(ND) * ND) using
%   getMobilitiesFromDopingConc for the Caughey-Thomas mobility.
%
%   Input:
%     rho_ohmcm  - resistivity [ohm.cm], scalar or array
%
%   Output:
%     ND         - donor concentration [cm^-3]

    q = 1.602e-19;  % [C]

    ND = zeros(size(rho_ohmcm));

    for ii = 1:numel(rho_ohmcm)
        rho = rho_ohmcm(ii);

        % Initial guess using low-doping mobility (~1330)
        nd = 1 / (q * 1330 * rho);

        % Iterate
        for iter = 1:50
            mu_n = getMobilitiesFromDopingConc(nd);
            nd_new = 1 / (q * mu_n * rho);
            if abs(nd_new - nd) / nd < 1e-10
                break;
            end
            nd = nd_new;
        end

        ND(ii) = nd;
    end
end