function Ef = getEfFromDopingConc(Nd, Na)
    % getEfFromDopingConc
    %   Computes Fermi level Ef (in eV relative to vacuum) from doping.
    %
    %   Usage:
    %       Ef = getEfFromDopingConc(NdSigned)
    %           NdSigned > 0 : n-type, Nd = NdSigned, Na = 0
    %           NdSigned < 0 : p-type, Na = |NdSigned|, Nd = 0
    %
    %       Ef = getEfFromDopingConc(Nd, Na)
    %           Nd, Na >= 0 : explicit donors/acceptors
    %
    %   All inputs in cm^-3, output Ef in eV (vacuum level = 0).
    %   Band edges: Ec = -4.05 eV, Ev = -5.17 eV.

    %----------------- constants -----------------
    k_B = 8.617333262e-5;   % eV/K
    T   = 300;              % K
    ni  = 1.5e10;           % cm^-3
    Nv  = 1.04e19;          % cm^-3
    Nc  = 2.8e19;           % cm^-3

    % Absolute band edges (relative to vacuum)
    Ec_abs = -4.05;         % eV, electron affinity
    Ev_abs = -5.17;         % eV, ionization potential

    % Bandgap (should be ~1.12 eV)
    Eg = Ec_abs - Ev_abs;   % eV

    % Intrinsic Fermi level relative to Ev:
    % Ei - Ev = Eg/2 + (3/4) kT ln(Nv/Nc)
    Ei_minus_Ev = Eg/2 + (3/4)*k_B*T*log(Nv/Nc);

    %-------------- handle input forms --------------
    if nargin == 1
        % Nd signed: >0 n-type, <0 p-type
        Nd_signed = Nd;
        Nd = zeros(size(Nd_signed));
        Na = zeros(size(Nd_signed));

        Nd(Nd_signed > 0) = Nd_signed(Nd_signed > 0);   % donors
        Na(Nd_signed < 0) = -Nd_signed(Nd_signed < 0);  % acceptors
    elseif nargin == 2
        % Nd, Na both given explicitly (nothing else to do)
    else
        error('getEfFromDopingConc: expects 1 or 2 input arguments.');
    end

    % ensure same size
    if ~isequal(size(Nd), size(Na))
        error('Nd and Na must have the same size.');
    end

    Ef = zeros(size(Nd));   % absolute Ef [eV]

    % Net doping:
    Nnet = Nd - Na;    % >0 → n-type, <0 → p-type, 0 → intrinsic

    % Indices
    idx_n   = Nnet > 0;
    idx_p   = Nnet < 0;
    idx_int = Nnet == 0;

    % n-type: n ≈ Nd - Na, Ef - Ei = kT ln(n/ni)
    if any(idx_n)
        n = Nnet(idx_n);
        dE = k_B*T*log(n ./ ni);   % Ef - Ei
        % Ef - Ec = (Ef - Ei) + (Ei - Ev) - Eg
        Ef_minus_Ec = dE + (Ei_minus_Ev - Eg);
        % Absolute Ef = Ec_abs + (Ef - Ec)
        Ef(idx_n) = Ec_abs + Ef_minus_Ec;
    end

    % p-type: p ≈ Na - Nd, Ef - Ei = -kT ln(p/ni)
    if any(idx_p)
        p = -Nnet(idx_p);          % positive
        dE = -k_B*T*log(p ./ ni);  % Ef - Ei
        Ef_minus_Ec = dE + (Ei_minus_Ev - Eg);
        Ef(idx_p) = Ec_abs + Ef_minus_Ec;
    end

    % intrinsic: Ef ≈ Ei
    if any(idx_int)
        % For intrinsic, Ef - Ec = (Ei - Ev) - Eg
        Ef_minus_Ec = Ei_minus_Ev - Eg;
        Ef(idx_int) = Ec_abs + Ef_minus_Ec;
    end
end