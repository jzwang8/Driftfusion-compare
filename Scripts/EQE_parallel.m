system('taskkill /F /IM EXCEL.EXE');
delete(gcp('nocreate'))
%%       Input:
file_name = 'input_files/pn_junction_jzw.csv';
scan_start = 400; %nm (must be between 301 - 1200 nm)
scan_end = 680;  %nm (must be between 301 - 1200 nm)
scan_interval = 8; % must be an integer value
Singlet_Fission = true;
SF_efficiency = 2;
charges = 1e19; % [cm-3]
% Save options:
save_TF = true;                             
[~, base_name, ~] = fileparts(file_name);
save_file_name = [datestr(now, 'yyyymmdd_HHMM') '_' base_name '_SF' num2str(Singlet_Fission) '_' num2str(SF_efficiency) '_' num2str(charges)];
% Parallel  calculations on no. cores:
number_of_workers = 12;

%% Initialisation:
initialise_df;
Wavelength_scan = scan_start:scan_interval:scan_end;
number_of_wavelengths = length(Wavelength_scan);
EQE_array = zeros(1,number_of_wavelengths);
IQE_array = zeros(1,number_of_wavelengths);
Jsc_array = zeros(1,number_of_wavelengths);

%% Calculate parameters:
par = pc(file_name);
par.eta = SF_efficiency;
%% Charged layer stuff (make into function)
layer_points_p_layer = par.parr(1);
% Get thickness of charged layer
% charges = 1e19; % [cm-3]
thickness = 1.5;  % [nm]
x = par.xx;
layer_index = find(x - (x(end) - thickness*1e-7) > 0, 1, 'first'); %finds index for the thickness of charged layer
actual_thickness = (x(end)-x(layer_index))*1e7;%compensates for the numerical grid that might not take exactly thickness as thickness
compensation_charges = -1 * actual_thickness/(x(layer_points_p_layer)*1e7) * charges;

par.dev_sub.Extra_charge(layer_index:end) = charges;
par.dev_sub.Extra_charge(1:layer_points_p_layer) = compensation_charges;

%% Caclulate equilibrium solution
soleq = equilibrate(par);

%% Calculate the number of complete chunks and the size of the remaining chunk
numCompleteChunks = floor(number_of_wavelengths / number_of_workers);
remainingIterations = mod(number_of_wavelengths, number_of_workers);

% Loop through each complete chunk
for chunk = 1:numCompleteChunks
    startIdx = (chunk - 1) * number_of_workers + 1;
    endIdx = chunk * number_of_workers;
    disp(['Processing chunk ', num2str(chunk), ': n = ', num2str(startIdx), ' to ', num2str(endIdx)]);

    parfor n = startIdx:endIdx
        % Calculate EQE data for the required wavelengths:
        wavelength = Wavelength_scan(n);
        [Jsc_array(n), EQE_array(n), IQE_array(n)] = EQE_parallel_loop(soleq, wavelength, Singlet_Fission);
    end
    delete(gcp('nocreate'))
    % Close all excel stuff
    system('taskkill /F /IM EXCEL.EXE');
end

delete(gcp('nocreate'))
% Store data in struct:
EQE_sol.EQE = EQE_array;
EQE_sol.Jsc = Jsc_array;
EQE_sol.wavelengths = Wavelength_scan;
EQE_sol.IQE = IQE_array;

% Plot the IQE spectrum:
figure
hold on
% plot(EQE_sol.wavelengths,EQE_sol.IQE, 'DisplayName','IQE')
plot(EQE_sol.wavelengths,EQE_sol.EQE, 'DisplayName','EQE')
xlabel('Wavelength (nm)')
ylabel('Quantum efficiency (%)')
legend show
hold off

% Save the data:
if save_TF == true
    rootFilePath = 'Results/';
    full_file_path = fullfile(rootFilePath,save_file_name);
    save(full_file_path,"EQE_sol")
end


function [Jsc, EQE_calc, IQE_calc] = EQE_parallel_loop(soleq, wavelength, Tetracene_TF)
    %% Constants
    h = 6.62607015e-34;
    e = 1.602176634e-19;
    c = 2.99792458e8;   

    %% Changing generation variables: (no need to recalculate eq. sols for every WL!)
    par_temp = soleq.el.par; % Copy the global parameters to change them locally in par_temp
    par_temp.Tetracene_TF = Tetracene_TF;
    disp(wavelength)
    par_temp.light_source1 = 'laser';
    par_temp.laser_lambda1 = wavelength;
    par_temp.pulsepow = 1; % [mW cm-2]
    par_temp.gx1 = generation(par_temp, 'laser', wavelength); % calculates generation profile based on the used wavelength
    soleq_temp = soleq.el;
    soleq_temp.par = par_temp; % write the changed pars in the temporary eq. solution

    %% Calculate CV for specific laser wavelength:
    Vmin = -0.3;
    Vmax = 0.6;
    scan_rate = 50e-3;
    cycles = 1;
    tpoints = 400;
    CVsol = doCV(soleq_temp, 1, 0, Vmax, Vmin, scan_rate, cycles, tpoints); % doCV uses parameters from soleq_temp.par

    %% Extract Jsc and calculate EQE for specific wavelength:
    V = dfana.calcVapp(CVsol);
    J = dfana.calcJ(CVsol, "sub").tot;
    change_sweep_direction_index = find(V == max(V),1);
    J_f = J(1:change_sweep_direction_index);
    V_f = V(1:change_sweep_direction_index);
    Jsc = abs(interp1(V_f, J_f, 0, 'linear')); % [A cm-2]
    Pin = par_temp.pulsepow * 1e-3;          %dfana.calcPin(CVsol); % [W cm-2]
    I = beerlambertI(par_temp, par_temp.xx, par_temp.light_source1, wavelength, 0);
    Plost = I(1,wavelength-300); % calculates the residue power
    % on the far left of the solar cell for the value of
    % laserlambda. The 301 comes from the beerlambertI
    % funcion where lambda=301:1200
    EQE_calc = 100 * Jsc * h * c / (Pin * wavelength * 1e-9 * e);
    IQE_calc = EQE_calc * Pin /(Pin - Plost);
end