% Calculate EQE w.r.t the wavelength, calculates equilibrium solution
% as well as performing a CV to obtain Jsc. The ratio Jsc/(Pin *
% lambda) gives the EQE for wavelengths between 300 and 1100 nm. The IQE
% is determined by accounting for the light exiting the cell on the
% back side.
% Input:    file_name       = file name of the input file one usually
%           save_TF         = logic if the output data should be saved
%           save_file_name  = file name under which the data will be
%                             saved. Change the rootFilePath in the code if
%                             you want to save it somewhere else.
% Output:   wavelengths     = array of scanned wavelengths [nm]
%           EQE             = array of resulting EQE [%]
%           IQE             = array of resulting IQE [%]
%           Jsc             = array of resulting Jsc [A cm-2]

%% Input:
file_name = 'input_files/No SRH';
save_TF = false;                             % save the output in a file true/false
save_file_name = '20250418 No SRH EQE Tc 0 eh per photon';
Singlet_Fission = false;

scan_start = 400; %nm (must be between 301 - 1200 nm)
scan_end = 500;  %nm (must be between 301 - 1200 nm)
scan_interval = 20;

%% Constants
h = 6.62607015e-34;
e = 1.602176634e-19;
c = 2.99792458e8;


%% initialize system
initialise_df;
par = pc(file_name);

Wavelength_scan = scan_start:scan_interval:scan_end;
number_of_wavelengths = length(Wavelength_scan);
EQE_array = zeros(1,number_of_wavelengths);
IQE_array = zeros(1,number_of_wavelengths);
Jsc_array = zeros(1,number_of_wavelengths);

%% Calculate equilibirum solution:
soleq = equilibrate(par); % calculates equilibrium solution

for n = 1 : number_of_wavelengths
    wavelength = Wavelength_scan(n); % [nm]
    par_temp = soleq.el.par; % Copy the global parameters to change them locally in par_temp
    par_temp.Tetracene_TF = Singlet_Fission;
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
    Jsc_array(n) = Jsc;
    Pin = par_temp.pulsepow * 1e-3;          %dfana.calcPin(CVsol); % [W cm-2]
    I = beerlambertI(par_temp, par_temp.xx, par_temp.light_source1, wavelength, 0);
    Plost = I(1,wavelength-300); % calculates the residue power
    % on the far left of the solar cell for the value of
    % laserlambda. The 301 comes from the beerlambertI
    % funcion where lambda=301:1200
    EQE_calc = 100 * Jsc * h * c / (Pin * wavelength * 1e-9 * e);
    IQE_calc = EQE_calc * Pin /(Pin - Plost);
    EQE_array(n) = EQE_calc;
    IQE_array(n) = IQE_calc;
end
EQE_sol.EQE = EQE_array;
EQE_sol.Jsc = Jsc_array;
EQE_sol.wavelengths = Wavelength_scan;
EQE_sol.IQE = IQE_array;
figure
plot(EQE_sol.wavelengths,EQE_sol.IQE)
xlabel('Wavelength (nm)')
ylabel('IQE (%)')
if save_TF == true
    rootFilePath = 'Results/';
    full_file_path = fullfile(rootFilePath,save_file_name);
    save(full_file_path,"EQE_sol")
end
