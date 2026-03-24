%% Combined JV + profiles + recombination + EQE script
% Uses ONE parameter set for both JV and EQE and saves everything in 'runs'

delete(gcp('nocreate'));
system('taskkill /F /IM EXCEL.EXE');  % optional cleanup of Excel

%% ------------------------------------------------------------------------
%  Settings
% -------------------------------------------------------------------------
file_name = 'Input_files/pn_junction_chargedlayer.csv';
run_name = 'ch0_SF1_eta2_ref0p15_sp42500_noSRH';

% SF/charged layer settings
charges          = 0;        % [cm^-3]
thickness        = 1.5;      % [nm]
Singlet_Fission  = true;
SF_efficiency    = 2;        % par.eta
reflection       = 0.15;      % par.kappa (JV uses this constant value)

% JV settings
light_intensity_JV = 1;    % suns
Vmin_JV            = -0.3;
Vmax_JV            =  0.7;
scan_rate_JV       = 50e-3;
cycles_JV          = 2;
tpoints_JV         = 300;

% EQE settings
scan_start        = 400;    % nm
scan_end          = 800;   % nm
scan_interval     = 5;      % nm
save_TF           = true;
number_of_workers = 12;     % desired workers for parpool

% Reflection spectrum for EQE (κ(λ)):
kappa_lambda = lambda_tmm;
% kappa_vals   = Reflection_tmm;
kappa_vals = reflection * ones(size(lambda_tmm));

%% ------------------------------------------------------------------------
%  Initialise Driftfusion and single parameter set
% -------------------------------------------------------------------------
initialise_df;

par = pc(file_name);
par.AbsTol        = 5e-6;
par.RelTol        = 5e-3;
par.MaxStepFactor = 1;

% Use ONE consistent parameter set for both JV and EQE:
par.eta          = SF_efficiency;     % singlet-fission-related parameter
par.Tetracene_TF = Singlet_Fission;   % same SF flag for JV and EQE
par.kappa        = reflection;        % JV uses constant kappa

% Apply charged layer once to this single par
par = add_charged_layer(par, charges, thickness);

% Generation profile for JV run (and default light source)
par.gx1 = generation(par, par.light_source1, par.laser_lambda1);

txt  = ['Injected hole flux = ', num2str(par.extra_holes)];
txt2 = ['Injected electron flux = ', num2str(par.extra_electrons)];

%% ------------------------------------------------------------------------
%  1) Equilibrium and JV / profiles / recombination
% -------------------------------------------------------------------------

% Equilibrium solution (used by both JV and EQE)
solEq = equilibrate(par);

% JV sweep using same par
solCV_JV = doCV(solEq.el, light_intensity_JV, ...
                0, Vmax_JV, Vmin_JV, scan_rate_JV, cycles_JV, tpoints_JV);

% Position index (middle of active layer)
xmesh_JV = solCV_JV.x;
ppos_JV  = getpointpos(par.d_midactive, xmesh_JV);

% Current density / Vapp
J_CV_JV = dfana.calcJ(solCV_JV, "sub");
Vapp_JV = dfana.calcVapp(solCV_JV);

% Time and x-range for profiles
tarr_JV   = solEq.el.t(end);
xrange_JV = [1.795*1e5, 1.803*1e5];  % adjust to taste

% --- Profiles ---

dfplot.npx(solEq.el, tarr_JV, xrange_JV);
dfplot.rhoxFxVx(solEq.el, tarr_JV, xrange_JV);
dfplot.ELx(solEq.el, tarr_JV, xrange_JV);
dfplot.rx(solCV_JV, tarr_JV, xrange_JV);

% --- JV and PCE ---

power_JV = Vapp_JV .* transpose(J_CV_JV.tot(:, ppos_JV));  % W cm^-2
MPP_JV   = min(power_JV);
Pin_JV   = dfana.calcPin(solCV_JV);
eff_JV   = abs(100 * (MPP_JV / Pin_JV));                   % [%]

figure;
hold on;
txt3 = ['PCE = ', num2str(eff_JV), ' %'];
plot(Vapp_JV, 1000 * J_CV_JV.tot(:, ppos_JV), 'DisplayName', txt3);
xlabel('Applied Voltage, V_{app} [V]');
ylabel('Current Density, J [mA cm^{-2}]');
set(legend, 'FontSize', 16);
set(legend, 'EdgeColor', [1 1 1]);
legend show;
hold off;

% Voc
Voc_JV = interp1(J_CV_JV.tot(:, ppos_JV), Vapp_JV, 0);
disp(['JV -> Voc: ', num2str(Voc_JV), ' V']);
disp(['JV -> PCE: ', num2str(eff_JV), ' %']);

%% ------------------------------------------------------------------------
%  2) EQE / IQE vs wavelength using the SAME parameter set
% -------------------------------------------------------------------------

Wavelength_scan       = scan_start:scan_interval:scan_end;
number_of_wavelengths = numel(Wavelength_scan);

EQE_array = zeros(1, number_of_wavelengths);
IQE_array = zeros(1, number_of_wavelengths);
Jsc_array = zeros(1, number_of_wavelengths);

% Start parallel pool
if isempty(gcp('nocreate'))
    parpool(number_of_workers);
end

disp('Starting EQE wavelength loop...');
parfor n = 1:number_of_wavelengths
    wavelength = Wavelength_scan(n);

    [Jsc_array(n), EQE_array(n), IQE_array(n)] = ...
        EQE_parallel_loop(solEq, wavelength, Singlet_Fission, ...
                          kappa_lambda, kappa_vals);
end

% Store EQE results in a struct
EQE_sol.EQE         = EQE_array;
EQE_sol.Jsc         = Jsc_array;
EQE_sol.wavelengths = Wavelength_scan;
EQE_sol.IQE         = IQE_array;

% Plot EQE
figure;
hold on;
plot(EQE_sol.wavelengths, EQE_sol.EQE, 'DisplayName', 'EQE');
xlabel('Wavelength (nm)');
ylabel('External Quantum Efficiency (%)');
legend show;
title('EQE spectrum');
hold off;

%% ------------------------------------------------------------------------
%  3) Append thisRun into 'runs' (workspace + Driftfusion/Results/runs_log.mat)
% -------------------------------------------------------------------------

% Build thisRun as you already do above
thisRun = struct();

% JV info
thisRun.Vapp      = Vapp_JV;
thisRun.Jtot      = J_CV_JV.tot(:, ppos_JV);
thisRun.solEq     = solEq;          % equilibrium used for both JV & EQE
thisRun.solCV     = solCV_JV;
thisRun.par       = par;
thisRun.Voc       = Voc_JV;
thisRun.PCE       = eff_JV;
thisRun.MPP_power = MPP_JV;
thisRun.Pin       = Pin_JV;

% Geometry / parameter (example: d(1,3))
thisRun.d_13  = par.d(1,3);  % [cm]
thisRun.label = sprintf('Run at d(1,3) = %.0f nm', thisRun.d_13 * 1e7);
thisRun.run_name = run_name;

% Metadata
thisRun.time = datetime('now');

% EQE info
thisRun.EQE_wavelengths = EQE_sol.wavelengths;
thisRun.EQE             = EQE_sol.EQE;
thisRun.IQE             = EQE_sol.IQE;
thisRun.Jsc_vs_lambda   = EQE_sol.Jsc;

%% ---- Determine where to start 'runs' from (workspace or file) ----

rootFilePath = fullfile('Driftfusion','Results');
if ~exist(rootFilePath,'dir')
    mkdir(rootFilePath);
end
runsLogFile = fullfile(rootFilePath,'runs_log.mat');

if evalin('base','exist(''runs'',''var'')')
    % 1) Use runs from workspace if it exists
    runs = evalin('base','runs');
elseif exist(runsLogFile,'file')
    % 2) If no runs in workspace, but log file exists -> load from file
    S    = load(runsLogFile,'runs');
    runs = S.runs;
else
    % 3) No runs in workspace, no log file -> start fresh
    runs = [];
end

%% ---- Append thisRun to runs, harmonising fields if needed ----

if isempty(runs)
    % First ever run
    runs = thisRun;

elseif ~isstruct(runs)
    warning('Existing ''runs'' is not a struct. Resetting it with thisRun.');
    runs = thisRun;

else
    % Harmonise fieldnames between existing runs and thisRun
    fExisting = fieldnames(runs);
    fNew      = fieldnames(thisRun);

    % Ensure thisRun has all fields existing runs have
    for k = 1:numel(fExisting)
        fn = fExisting{k};
        if ~isfield(thisRun, fn)
            thisRun.(fn) = [];
        end
    end

    % Ensure existing runs have all fields thisRun has
    for k = 1:numel(fNew)
        fn = fNew{k};
        if ~isfield(runs, fn)
            [runs(1:end).(fn)] = deal([]);
        end
    end

    % Now they share the same field set -> safe to append
    runs(end+1) = thisRun;  %#ok<SAGROW>
end

%% ---- Push to workspace and save to MAT file ----

assignin('base','runs',runs);
disp(['Stored run #', num2str(numel(runs)), ' in variable ''runs'' in workspace.']);

save(runsLogFile,'runs','-v7.3');
disp(['Saved updated runs (', num2str(numel(runs)), ' entries) to ', runsLogFile]);


%% ------------------------------------------------------------------------
%  Local helper functions
% -------------------------------------------------------------------------

function par_out = add_charged_layer(par_in, charges, thickness_nm)
    par = par_in;
    layer_points_p_layer = par.parr(1);
    x = par.xx;
    
    layer_index = find(x - (x(end) - thickness_nm*1e-7) > 0, 1, 'first');
    
    actual_thickness = (x(end) - x(layer_index)) * 1e7;
    compensation_charges = -1 * actual_thickness / (x(layer_points_p_layer)*1e7) * charges;
    
    par.dev_sub.Extra_charge(layer_index:end)        = charges;
    par.dev_sub.Extra_charge(1:layer_points_p_layer) = compensation_charges;
    
    par_out = par;
end

function [Jsc, EQE_calc, IQE_calc] = EQE_parallel_loop( ...
            soleq, wavelength, Tetracene_TF, kappa_lambda, kappa_vals)
    % EQE_parallel_loop: compute Jsc, EQE, IQE at a given wavelength
    %
    % Uses the SAME base parameter set stored in soleq.el.par,
    % but with wavelength-dependent reflection kappa(λ).

    % Constants
    h = 6.62607015e-34;      % Planck [J s]
    e = 1.602176634e-19;     % elementary charge [C]
    c = 2.99792458e8;        % speed of light [m s^-1]

    % Local copy of parameters
    par_temp = soleq.el.par;
    par_temp.Tetracene_TF = Tetracene_TF;

    % --- wavelength-dependent kappa(λ) for THIS wavelength ---
    kappa_eff = interp1(kappa_lambda, kappa_vals, ...
                        wavelength, 'linear', 'extrap');
    par_temp.kappa = kappa_eff;

    disp(['EQE wavelength = ', num2str(wavelength), ...
          ' nm, kappa_eff = ', num2str(kappa_eff)]);

    par_temp.light_source1 = 'laser';
    par_temp.laser_lambda1 = wavelength;
    par_temp.pulsepow      = 1;  % [mW cm^-2]

    % Generation at this wavelength (uses par_temp.kappa inside generation.m)
    par_temp.gx1 = generation(par_temp, 'laser', wavelength);

    % Pack into temporary equilibrium solution
    soleq_temp     = soleq.el;
    soleq_temp.par = par_temp;

    % CV for specific wavelength
    Vmin     = -0.3;
    Vmax     =  0.6;
    scan_rate = 50e-3;
    cycles    = 1;
    tpoints   = 400;
    CVsol = doCV(soleq_temp, 1, 0, Vmax, Vmin, scan_rate, cycles, tpoints);

    % Extract Jsc and calculate EQE
    V = dfana.calcVapp(CVsol);
    J = dfana.calcJ(CVsol, "sub").tot;
    
    change_sweep_direction_index = find(V == max(V), 1);
    J_f = J(1:change_sweep_direction_index);
    V_f = V(1:change_sweep_direction_index);
    
    Jsc = abs(interp1(V_f, J_f, 0, 'linear'));  % [A cm^-2]

    Pin = par_temp.pulsepow * 1e-3;            % [W cm^-2]

    % Residual power at left side (same as before)
    I = beerlambertI(par_temp, par_temp.xx, par_temp.light_source1, wavelength, 0);
    Plost = I(1, wavelength-300);  % beerlambertI uses lambda=301:1200

    % External QE
    EQE_calc = 100 * Jsc * h * c / (Pin * wavelength * 1e-9 * e);
    % Internal QE
    IQE_calc = EQE_calc * Pin / (Pin - Plost);
end