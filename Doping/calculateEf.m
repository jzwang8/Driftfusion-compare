% Gaussian doping parameters
% N0 = gaussFit.N0; % Peak doping concentration in cm^-3
N0 = gaussFit.N0/10; % Peak doping concentration in cm^-3
x0 = gaussFit.x0;       % Mean position (center of the Gaussian) in arbitrary units
sigma = gaussFit.sigma;    % Standard deviation of the Gaussian

% Generate 10 evenly spaced x positions
n=100;
depth = 300;
layerPointsTotal = 5000;
x = linspace(0, depth, n); % n evenly spaced points from 0 to 300

% Get exponential doping parameters
[L,A] = exponentialFromGaussian(N0,x0,sigma);

% Compute doping concentration at each x
Nd_gaussian = N0 * exp(-((x - x0).^2) / (2 * sigma^2));
Nd_exp = A * exp(-x / L);

Ef_Eg_gaussian = getEfFromDopingConc(Nd_gaussian);
Ef_Eg_exp = getEfFromDopingConc(Nd_exp);

thicknesses = repmat(depth/n*1e-7, length(x), 1);
layerPoints = repmat(floor(layerPointsTotal/n), length(x), 1);

% Save results to CSV
% Define column headers
headers = ["x", "thickness", "Nd_gaussian", "Ef_Eg_gaussian", "Nd_exp", "Ef_Eg_exp", "layerPoints"];

% Combine headers and data
data = [flip(x(:)), thicknesses, flip(Nd_gaussian(:)), flip(Ef_Eg_gaussian(:)), flip(Nd_exp(:)), flip(Ef_Eg_exp(:)), layerPoints];

% Write to CSV
writematrix([headers; num2cell(data)], 'fermi_levels.csv');

fprintf('Data saved to fermi_levels.csv\n');
