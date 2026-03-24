% Define constants
q = 1.6e-19; % Elementary charge (C)
eps0 = 8.85e-14; % Vacuum permittivity (F/cm)
T = 300; % Temperature in Kelvin (Assumed Room Temperature)

% User inputs
Rs = input('Enter measured sheet resistance (ohm/sq): '); % Measured by 4-point probe
xj = input('Enter junction depth (cm): '); % Known doping depth

% Define empirical electron mobility function based on Arora model
% Mobility function for n-type Si (approximation)
mu_n = @(Nd) 88 * (Nd / 1e17).^(-0.68) + 125; % Mobility (cm^2/Vs) for Nd in cm^-3

% Solve for donor concentration N_D iteratively
Nd_guess = 1e15; % Initial guess (cm^-3)
tolerance = 1e-3;
error = 1;

while error > tolerance
    rho = Rs * xj; % Resistivity (ohm-cm)
    Nd_new = 1 / (q * mu_n(Nd_guess) * rho); % New estimate for Nd
    error = abs((Nd_new - Nd_guess) / Nd_guess); % Relative err1or
    Nd_guess = Nd_new; % Update estimate
end

% Display results
fprintf('Estimated n-type emitter doping concentration: %.2e cm^-3\n', Nd_new);