% clear all
file_name = 'input_files/pn_junction_var10';
save_TF = false;                             % save the output in a file true/false
save_file_name = 'test';
tetracene = false;
% Calculate PL w.r.t the wavelength, calculates equilibrium solution
% as well as performing a CV to obtain Jsc. For wavelengths between 300 and 1100 nm.
% Input:    file_name       = file name of the input file one usually
%           save_TF         = logic if the output data should be saved
%           save_file_name  = file name under which the data will be
%                             saved. Change the rootFilePath in the code if
%                             you want to save it somewhere else.
% Output:   PL.wavelengths     = array of scanned wavelengths [nm]
%           PL.I= array of resulting PL intensity [some unit]

%% initialize system
initialise_df;
par = pc(file_name);


scan_start = 400; %nm (must be between 301 - 1200 nm)
scan_end = 800;  %nm (must be between 301 - 1200 nm)
scan_interval = 20;

Wavelength_scan = scan_start:scan_interval:scan_end;
number_of_wavelengths = length(Wavelength_scan);
I_array = zeros(1,number_of_wavelengths);

%% Calculate equilibirum solution:
soleq = equilibrate(par); % calculates equilibrium solution

for n = 1 : number_of_wavelengths
    wavelength = Wavelength_scan(n); % [nm]
    disp(wavelength)

    par_temp = soleq.el.par;
    par_temp.light_source1 = 'laser';
    par_temp.laser_lambda1 = wavelength;
    par_temp.Tetracene_TF = tetracene;
    par_temp.pulsepow = 1; %1e24*h*c/wavelength * 1e3;        % [mW cm-2]
    % par.g_inj = 1e24;       % [cm-3s-1]
    par_temp.gx1 = generation(par_temp, par_temp.light_source1, par_temp.laser_lambda1); % calculates generation profile based on the used wavelength
    soleq_temp = soleq.el;
    soleq_temp.par = par_temp;

    %% Turn on light:
    Vmin = -0.3;
    Vmax = 0.6;
    scan_rate = 50e-3;
    cycles = 1;
    Rs = 0;
    mobseti = 0;
    light_intensity = 1; % [suns]
    tpoints = 10;
    stable_time = -1;
    sol = lightonRs(soleq_temp, light_intensity, stable_time, mobseti, Rs, tpoints);
    %% Extract r_rad and calculate EQE for specific wavelength:
    mesh_option = 'whole'; % Choose 'whole' or 'sub'
    [r, ns, ps, alpha_xn, beta_xp] = dfana.calcr(sol, mesh_option);

    r_rad = r.btb;
    x = par.xx;
    I_PL = trapz(x, r_rad(end,:));
    I_array(n) = I_PL;
end
PL.I= I_array;
PL.wavelengths = Wavelength_scan;

%% Plot solution spectrum:
figure
plot(PL.wavelengths,PL.I)
xlabel('Wavelength (nm)')
ylabel('PL (a.u.)')

%% Save data:
if save_TF == true
    rootFilePath = 'Results/';
    full_file_path = fullfile(rootFilePath,save_file_name);
    save(full_file_path,"PL")
end
