function Qf_q_cm2 = calcQfFromVi(Vi)
% calc_Qf_from_Vi  Compute fixed charge density Qf from Vi
%
%   Qf_q_cm2 = calc_Qf_from_Vi(Vi, xc, Ki)
%
%   Inputs:
%       Vi  - interface potential (eV), can be scalar or vector
%       xc  - characteristic length (m) between charge and interface
%       Ki  - relative permittivity of the insulator (dimensionless)
%
%   Output:
%       Qf_q_cm2  - fixed charge density (q/cm^2)

%
%   Make sure units are consistent:
%       Vi in volts, xc in meters, eps0 in F/m.
    
    q=1.6e-19;
    Vi_volts = Vi; % in V

    xc=.75e-9;
    Ki=1.6;

    % Vacuum permittivity (F/m)
    eps0 = 8.854187817e-12;


    % Compute Qf
    Qf_C_m2 = Vi_volts .* Ki .* eps0 ./ xc;
    Qf_C_cm2 = Qf_C_m2*1e-4;
    Qf_q_cm2 = Qf_C_cm2/q;
end