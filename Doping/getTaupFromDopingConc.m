function tau_p = getTaupFromDopingConc(Nd_cm3)
    Cp      = 2e-31;  % cm^6/s
          % Cp follows Fossum (1983) Eq.1: tau_p = 1/(Cp * ND^2)
% This is the eeh Auger coefficient (called Cn in standard R = Cn*n^2*p + Cp*n*p^2)
    tau_SRH = 1e-3;   % s

    if any(Nd_cm3 <= 0, 'all')
        error('Nd_cm3 must be positive.');
    end

    tau_p = 1 ./ ( 1./tau_SRH + Cp .* Nd_cm3.^2 );
end