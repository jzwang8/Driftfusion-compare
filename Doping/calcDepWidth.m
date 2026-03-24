eps_0 = 8.85e-14;
eps_si = 11.7;
q = 1.6e-19;
N_A = 2e15;
N_D = 8e17;
Vbi = .87;
V=0;

W = sqrt(2*eps_si*eps_0/q * (N_A+N_D)/(N_A*N_D)*(Vbi-V))
W_p = N_D/(N_A+N_D)*W
W_n = N_A/(N_A+N_D)*W

% N_D=5.11e18;
% n=N_D;
kT=.0259;
% Nc=1.8e19;
% Ef_minus_Ec_fermi = kT*(log(n/Nc)+2^(-3/2)*n/Nc)
% 
% Ec_minus_Ef_boltz = kT*log(Nc/n)

Vs = 0.14; % Surface potential (V)

% W = sqrt(2 * eps_si*eps_0 * Vs / (q * N_D));
% W_nm = W * 1e7

kT_q = (k*T)/q;               % Thermal voltage (V)
epsilon_0 = 8.854e-14;         % Vacuum permittivity (F/cm)
epsilon_si = 11.7 * epsilon_0; % Permittivity of silicon (F/cm)

% Inputs
ND = 1e20;                     % Donor concentration (cm^-3)
Ec_Ef = -0.17;                 % (eV) Ec - Ef (Fermi level is above Ec by 0.17 eV)

% Built-in potential (V)
% Assuming surface is intrinsic or pinned near midgap
Eg = 1.12;                     % Bandgap of silicon at 300K (eV)
Vbi = Eg/2 - Ec_Ef;             % (eV)
Vbi = Vbi;                     % already in eV units

% Calculation of depletion width
W = sqrt(2 * epsilon_si * Vbi / (q * ND)); % Depletion width in cm

% Convert to nm
W_nm = W * 1e7; % (1 cm = 10^7 nm)