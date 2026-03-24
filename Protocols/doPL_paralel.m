clear all
system('taskkill /F /IM EXCEL.EXE');
delete(gcp('nocreate'))
% Input:
file_name = 'input_files/pn_junction_var25';
scan_start = 490; %nm (must be between 301 - 1200 nm)
scan_end = 640;  %nm (must be between 301 - 1200 nm)
scan_interval = 3; % must be an integer value
Tetracene_TF = false;
% Save options:
save_TF = false;                             
save_file_name = '20240723 pn_junction_var10 cations';
% Paralel  calculations on no. cores:
number_of_workers = 24;

%% Initialisation:
initialise_df;
Wavelength_scan = scan_start:scan_interval:scan_end;
number_of_wavelengths = length(Wavelength_scan);
PL_array = zeros(1,number_of_wavelengths);

%% Calculate equilibrium solution:
par = pc(file_name);
soleq = equilibrate(par,1); % calculates equilibrium solution

%% Calculate the number of complete chunks and the size of the remaining chunk
numCompleteChunks = floor(number_of_wavelengths / number_of_workers);
remainingIterations = mod(number_of_wavelengths, number_of_workers);

% Loop through each complete chunk
for chunk = 1:numCompleteChunks
    startIdx = (chunk - 1) * number_of_workers + 1;
    endIdx = chunk * number_of_workers;
    disp(['Processing chunk ', num2str(chunk), ': n = ', num2str(startIdx), ' to ', num2str(endIdx)]);

    parfor n = startIdx:endIdx
        % Calculate PL data for the required wavelengths:
        wavelength = Wavelength_scan(n);
        disp(wavelength);
        PL_array(n) = PL_paralel_loop(soleq, wavelength, Tetracene_TF);
    end
    delete(gcp('nocreate'))
    % Close all excel stuff
    system('taskkill /F /IM EXCEL.EXE');
end

delete(gcp('nocreate'))
% Store data in struct:
PL.I_PL = PL_array;
PL.wavelengths = Wavelength_scan;

% Plot the PL spectrum:
figure
plot(PL.wavelengths,PL.I_PL)
xlabel('Wavelength (nm)')
ylabel('Intensity (a.u.)')

% Save the data:
if save_TF == true
    rootFilePath = 'Results/';
    full_file_path = fullfile(rootFilePath,save_file_name);
    save(full_file_path,"EQE_sol")
end


function I_PL = PL_paralel_loop(soleq, wavelength, tetracene_TF)
    par_temp = soleq.el.par;
    par_temp.light_source1 = 'laser';
    par_temp.laser_lambda1 = wavelength;
    par_temp.Tetracene_TF = tetracene_TF;
    par_temp.pulsepow = 1; %1e24*h*c/wavelength * 1e3;        % [mW cm-2]
    par_temp.gx1 = generation(par_temp, par_temp.light_source1, par_temp.laser_lambda1); % calculates generation profile based on the used wavelength
    soleq_temp = soleq.el;
    soleq_temp.par = par_temp;

    %% Turn on light:
    Rs = 0;
    mobseti = 0;
    int1 = 1; % Maximum fraction of used generation rate i.e. how bright you set your light source.
    tpoints = 10;
    stable_time = -1;
    sol = lightonRs(soleq_temp, int1, stable_time, mobseti, Rs, tpoints);
    %% Extract r_rad and calculate EQE for specific wavelength:
    mesh_option = 'whole'; % Choose 'whole' or 'sub'
    [r, ns, ps, alpha_xn, beta_xp] = dfana.calcr(sol, mesh_option);

    r_rad = r.btb;
    x = par_temp.xx;
    I_PL = trapz(x, r_rad(end,:));

end