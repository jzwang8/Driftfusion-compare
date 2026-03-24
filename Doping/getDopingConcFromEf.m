function [Nd, Na] = getDopingConcFromEf(Ef_abs)
    % getDopingConcFromEf
    %   Given an absolute Fermi level Ef_abs (eV, relative to vacuum),
    %   compute the equivalent donor- and acceptor-like concentrations.
    %
    %   Inputs:
    %       Ef_abs  - Fermi level in eV (vacuum level = 0)
    %
    %   Outputs:
    %       Nd      - donor (n-type) concentration [cm^-3]
    %       Na      - acceptor (p-type) concentration [cm^-3]
    %
    %   Model:
    %       n  = ni * exp((Ef - Ei) / (k_B T))
    %       p  = ni * exp((Ei - Ef) / (k_B T))
    %
    %   Here we place the bands at:
    %       Ec_abs = -4.05 eV   (electron affinity)
    %       Ev_abs = -5.17 eV   (ionization potential)
    %   so Eg = Ec_abs - Ev_abs.

    %----------------- constants -----------------
    k_B = 8.617333262e-5;   % eV/K
    T   = 300;              % K
    ni  = 1.5e10;           % cm^-3
    Nv  = 1.04e19;          % cm^-3
    Nc  = 2.8e19;           % cm^-3

    % Absolute band edges (relative to vacuum)
    Ec_abs = -4.05;         % eV, conduction band edge (electron affinity)
    Ev_abs = -5.17;         % eV, valence band edge (ionization potential)

    % Bandgap (should be ~1.12 eV)
    Eg = Ec_abs - Ev_abs;   % eV

    % Intrinsic Fermi level relative to Ev:
    % Ei - Ev = Eg/2 + (3/4) kT ln(Nv/Nc)
    Ei_minus_Ev = Eg/2 + (3/4)*k_B*T*log(Nv/Nc);

    % Absolute intrinsic Fermi level:
    Ei_abs = Ev_abs + Ei_minus_Ev;   % eV

    %----------------- carrier concentrations -----------------
    % Ef - Ei (absolute energies)
    delta = Ef_abs - Ei_abs;         % eV

    % n-type (electron) concentration
    Nd = ni .* exp( delta / (k_B * T) );    % [cm^-3]

    % p-type (hole) concentration
    Na = ni .* exp( -delta / (k_B * T) );   % [cm^-3]
end