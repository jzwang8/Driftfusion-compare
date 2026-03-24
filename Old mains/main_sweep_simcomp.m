%% Combined JV + profiles + recombination + EQE script
% Uses ONE parameter set for both JV and EQE and saves everything in 'runs'
% Now with a 1D parameter sweep capability.

delete(gcp('nocreate'));
system('taskkill /F /IM EXCEL.EXE');  % optional cleanup of Excel

%% ------------------------------------------------------------------------
%  Base Settings (defaults that will be used unless swept)
% -------------------------------------------------------------------------
file_name = 'Input_files/pn_junction_chargedlayer.csv';
run_name = 'ch0_SF1_eta2_ref0p15_sp1000';

% --- Base SF/charged layer settings (defaults) ---
charges_base          = 0;        % [cm^-3]
junction_depth_base   = 300; % [nm]
Singlet_Fission_base  = true;
SF_efficiency_base    = 2;        % par.eta
sp_r_base             = 1000; % SRV [cm/s]

% other parameters to change but not sweep
thickness        = 1.5;      % [nm]
reflection       = 0.15;      % par.kappa (JV uses this constant value)
% TODO: figure out kappa(JV lambda)

% Reflection spectrum for EQE (κ(λ)):
kappa_lambda = lambda_tmm;                 % from your TMM
kappa_vals = reflection * ones(size(lambda_tmm));
% kappa_vals = Reflection_tmm;

% JV settings (kept constant in this script)
light_intensity_JV = 1;    % suns
Vmin_JV            = -0.3;
Vmax_JV            =  0.7;
scan_rate_JV       = 50e-3;
cycles_JV          = 2;
tpoints_JV         = 300;

% EQE settings (kept constant)
scan_start        = 400;    % nm
scan_end          = 800;    % nm
scan_interval     = 5;      % nm
save_TF           = true;
number_of_workers = 8;     % desired workers for parpool

%% ------------------------------------------------------------------------
%  Choose which variable to sweep
% -------------------------------------------------------------------------
% Options (as implemented below): 'charges', 'junction_depth', 'SF_flag',
% 'eta_sf', 'reflection', 'sp_r'
sweepParam  = 'junction_depth';              % <--- change this
sweepValues = [200,300,400];              % <--- and this (example: sweep charges 0 and 1e19)

%% ------------------------------------------------------------------------
%  Initialise Driftfusion once
% -------------------------------------------------------------------------
initialise_df;

%% ------------------------------------------------------------------------
%  MAIN SWEEP LOOP
% -------------------------------------------------------------------------
for iSweep = 1:numel(sweepValues)
    fprintf('\n===== Sweep %d / %d: %s = %g =====\n', ...
            iSweep, numel(sweepValues), sweepParam, sweepValues(iSweep));

    % --- Reset to base values for this sweep iteration ---
    charges         = charges_base;
    junction_depth  = junction_depth_base;
    Singlet_Fission = Singlet_Fission_base;
    SF_efficiency   = SF_efficiency_base;
    sp_r            = sp_r_base;

    % --- Override just the swept parameter ---
    switch sweepParam
        case 'charges'
            charges = sweepValues(iSweep);
        
        case 'junction_depth'
            junction_depth = sweepValues(iSweep);

        case 'SF_flag'
            % Expect sweepValues to be [false true] or [0 1]
            Singlet_Fission = logical(sweepValues(iSweep));

        case 'eta_sf'
            SF_efficiency = sweepValues(iSweep);

        case 'sp_r'
            sp_r = sweepValues(iSweep);

        otherwise
            error('Unknown sweepParam "%s".', sweepParam);
    end

    %% --------------------------------------------------------------------
    %  Build parameter set for THIS sweep value
    % ---------------------------------------------------------------------
    par = pc(file_name);
    par.AbsTol        = 5e-6;
    par.RelTol        = 5e-3;
    par.MaxStepFactor = 1;
    
    par.d(1,3) = junction_depth*1e-9*1e2;% [nm]
    par.sp_r = sp_r;

    % Use ONE consistent parameter set for both JV and EQE:
    par.eta          = SF_efficiency;     % singlet-fission-related parameter
    par.Tetracene_TF = Singlet_Fission;   % same SF flag for JV and EQE
    par.kappa        = reflection;        % JV uses constant kappa

    % Apply charged layer once to this single par
    par = add_charged_layer(par, charges, thickness);

    % Generation profile for JV run (and default light source)
    par.gx1 = generation(par, par.light_source1, par.laser_lambda1);

    par = refresh_device(par);

    txt  = ['Injected hole flux = ', num2str(par.extra_holes)];
    txt2 = ['Injected electron flux = ', num2str(par.extra_electrons)];

    %% --------------------------------------------------------------------
    %  1) Equilibrium and JV / profiles / recombination
    % ---------------------------------------------------------------------
    solEq = equilibrate(par);

    solCV_JV = doCV(solEq.el, light_intensity_JV, ...
                    0, Vmax_JV, Vmin_JV, scan_rate_JV, cycles_JV, tpoints_JV);

    xmesh_JV = solCV_JV.x;
    ppos_JV  = getpointpos(par.d_midactive, xmesh_JV);

    J_CV_JV = dfana.calcJ(solCV_JV, "sub");
    Vapp_JV = dfana.calcVapp(solCV_JV);

    tarr_JV   = solEq.el.t(end);
    xrange_JV = [1.795*1e5, 1.803*1e5];  % adjust to taste

    % Profiles (comment out if you don’t want lots of figures)
    dfplot.npx    (solEq.el, tarr_JV, xrange_JV);
    dfplot.rhoxFxVx(solEq.el, tarr_JV, xrange_JV);
    dfplot.ELx    (solEq.el, tarr_JV, xrange_JV);
    dfplot.rx     (solCV_JV, tarr_JV, xrange_JV);

    power_JV = Vapp_JV .* transpose(J_CV_JV.tot(:, ppos_JV));  % W cm^-2
    MPP_JV   = min(power_JV);
    Pin_JV   = dfana.calcPin(solCV_JV);
    eff_JV   = abs(100 * (MPP_JV / Pin_JV));                   % [%]

    figure;
    hold on;
    txt3 = sprintf('PCE = %.3f %%', eff_JV);
    plot(Vapp_JV, 1000 * J_CV_JV.tot(:, ppos_JV), 'DisplayName', txt3);
    xlabel('Applied Voltage, V_{app} [V]');
    ylabel('Current Density, J [mA cm^{-2}]');
    set(legend, 'FontSize', 16, 'EdgeColor', [1 1 1]);
    legend show;
    hold off;

    Voc_JV = interp1(J_CV_JV.tot(:, ppos_JV), Vapp_JV, 0);
    disp(['JV -> Voc: ', num2str(Voc_JV), ' V']);
    disp(['JV -> PCE: ', num2str(eff_JV), ' %']);

    %% --------------------------------------------------------------------
    %  2) EQE / IQE vs wavelength using the SAME parameter set
    % ---------------------------------------------------------------------
    Wavelength_scan       = scan_start:scan_interval:scan_end;
    number_of_wavelengths = numel(Wavelength_scan);

    EQE_array = zeros(1, number_of_wavelengths);
    IQE_array = zeros(1, number_of_wavelengths);
    Jsc_array = zeros(1, number_of_wavelengths);

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

    EQE_sol.EQE         = EQE_array;
    EQE_sol.Jsc         = Jsc_array;
    EQE_sol.wavelengths = Wavelength_scan;
    EQE_sol.IQE         = IQE_array;

    figure;
    hold on;
    plot(EQE_sol.wavelengths, EQE_sol.EQE, 'DisplayName', 'EQE');
    xlabel('Wavelength (nm)');
    ylabel('External Quantum Efficiency (%)');
    legend show;
    title(sprintf('EQE spectrum (%s = %g)', sweepParam, sweepValues(iSweep)));
    hold off;

    %% --------------------------------------------------------------------
    %  Build parameter set for THIS sweep value
    % ---------------------------------------------------------------------
    par = pc(file_name);
    par.AbsTol        = 5e-6;
    par.RelTol        = 5e-3;
    par.MaxStepFactor = 1;
    
    par.d(1,3) = junction_depth*1e-9*1e2;% [nm]
    par.sp_r = sp_r;

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

    %% --------------------------------------------------------------------
    %  1) Equilibrium and JV / profiles / recombination
    % ---------------------------------------------------------------------
    solEq = equilibrate(par);

    solCV_JV = doCV(solEq.el, light_intensity_JV, ...
                    0, Vmax_JV, Vmin_JV, scan_rate_JV, cycles_JV, tpoints_JV);

    xmesh_JV = solCV_JV.x;
    ppos_JV  = getpointpos(par.d_midactive, xmesh_JV);

    J_CV_JV = dfana.calcJ(solCV_JV, "sub");
    Vapp_JV = dfana.calcVapp(solCV_JV);

    tarr_JV   = solEq.el.t(end);
    xrange_JV = [1.795*1e5, 1.803*1e5];  % adjust to taste

    % Profiles (comment out if you don’t want lots of figures)
    dfplot.npx    (solEq.el, tarr_JV, xrange_JV);
    dfplot.rhoxFxVx(solEq.el, tarr_JV, xrange_JV);
    dfplot.ELx    (solEq.el, tarr_JV, xrange_JV);
    dfplot.rx     (solCV_JV, tarr_JV, xrange_JV);

    power_JV = Vapp_JV .* transpose(J_CV_JV.tot(:, ppos_JV));  % W cm^-2
    MPP_JV   = min(power_JV);
    Pin_JV   = dfana.calcPin(solCV_JV);
    eff_JV   = abs(100 * (MPP_JV / Pin_JV));                   % [%]

    figure;
    hold on;
    txt3 = sprintf('PCE = %.3f %%', eff_JV);
    plot(Vapp_JV, 1000 * J_CV_JV.tot(:, ppos_JV), 'DisplayName', txt3);
    xlabel('Applied Voltage, V_{app} [V]');
    ylabel('Current Density, J [mA cm^{-2}]');
    set(legend, 'FontSize', 16, 'EdgeColor', [1 1 1]);
    legend show;
    hold off;

    Voc_JV = interp1(J_CV_JV.tot(:, ppos_JV), Vapp_JV, 0);
    disp(['JV -> Voc: ', num2str(Voc_JV), ' V']);
    disp(['JV -> PCE: ', num2str(eff_JV), ' %']);

    %% --------------------------------------------------------------------
    %  2) EQE / IQE vs wavelength using the SAME parameter set
    % ---------------------------------------------------------------------
    Wavelength_scan       = scan_start:scan_interval:scan_end;
    number_of_wavelengths = numel(Wavelength_scan);

    EQE_array = zeros(1, number_of_wavelengths);
    IQE_array = zeros(1, number_of_wavelengths);
    Jsc_array = zeros(1, number_of_wavelengths);

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

    EQE_sol.EQE         = EQE_array;
    EQE_sol.Jsc         = Jsc_array;
    EQE_sol.wavelengths = Wavelength_scan;
    EQE_sol.IQE         = IQE_array;

    figure;
    hold on;
    plot(EQE_sol.wavelengths, EQE_sol.EQE, 'DisplayName', 'EQE');
    xlabel('Wavelength (nm)');
    ylabel('External Quantum Efficiency (%)');
    legend show;
    title(sprintf('EQE spectrum (%s = %g)', sweepParam, sweepValues(iSweep)));
    hold off;

    %% --------------------------------------------------------------------
    %  3) Build thisRun for this sweep value and append to runs
    % ---------------------------------------------------------------------
    thisRun = struct();

    % JV info
    thisRun.Vapp      = Vapp_JV;
    thisRun.Jtot      = J_CV_JV.tot(:, ppos_JV);
    thisRun.solEq     = solEq;
    thisRun.solCV     = solCV_JV;
    thisRun.par       = par;

    thisRun.Voc       = Voc_JV;
    thisRun.PCE       = eff_JV;
    thisRun.MPP_power = MPP_JV;
    thisRun.Pin       = Pin_JV;

    % Human-readable identifiers
    thisRun.run_name = run_name;
    thisRun.label = sprintf('%s = %g', sweepParam, sweepValues(iSweep));

    % Geometry / parameter
    thisRun.d_13  = par.d(1,3);  % [cm]

    % Metadata
    thisRun.time        = datetime('now');
    thisRun.sweepParam  = sweepParam;
    thisRun.sweepValue  = sweepValues(iSweep);

    % EQE info
    thisRun.EQE_wavelengths = EQE_sol.wavelengths;
    thisRun.EQE             = EQE_sol.EQE;
    thisRun.IQE             = EQE_sol.IQE;
    thisRun.Jsc_vs_lambda   = EQE_sol.Jsc;

       % ---- Append thisRun to / create 'runs' in base workspace (field-agnostic) ----
    if evalin('base','exist(''runs_simcomp'',''var'')')
        runs = evalin('base','runs_simcomp');
    else
        runs = struct([]);   % 0x0 empty struct array
    end

    if isempty(runs)
        % First ever run, just start the array
        runs = thisRun;

    elseif ~isstruct(runs)
        % Something weird in workspace, reset
        warning('Variable "runs_simcomp" exists but is not a struct; resetting it with thisRun.');
        runs = thisRun;

    else
        % --- 1) Build the union of all field names ---
        fRuns = fieldnames(runs);
        fNew  = fieldnames(thisRun);
        allFields = union(fRuns, fNew);   % cell array of all unique field names

        % --- 2) Make sure *every* field exists in both runs and thisRun ---
        for k = 1:numel(allFields)
            fn = allFields{k};

            % If runs is missing this field, add it as [] to all existing elements
            if ~isfield(runs, fn)
                [runs(1:end).(fn)] = deal([]);
            end

            % If thisRun is missing this field, add it as []
            if ~isfield(thisRun, fn)
                thisRun.(fn) = [];
            end
        end

        % --- 3) Force the SAME field order in both struct arrays ---
        % Build a template struct with fields in the desired order
        template = cell2struct(cell(size(allFields)), allFields, 1);
        runs    = orderfields(runs, template);
        thisRun = orderfields(thisRun, template);

        % --- 4) Safe append ---
        runs(end+1) = thisRun; %#ok<SAGROW>
    end

    % Push back to base workspace
    assignin('base','runs_simcomp',runs);
    fprintf('Stored run #%d in variable ''runs_simcomp'' in workspace.\n', numel(runs));
end

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

    h = 6.62607015e-34;
    e = 1.602176634e-19;
    c = 2.99792458e8;

    par_temp = soleq.el.par;
    par_temp.Tetracene_TF = Tetracene_TF;

    kappa_eff    = interp1(kappa_lambda, kappa_vals, wavelength, 'linear', 'extrap');
    par_temp.kappa = kappa_eff;

    disp(['EQE wavelength = ', num2str(wavelength), ...
          ' nm, kappa_eff = ', num2str(kappa_eff)]);

    par_temp.light_source1 = 'laser';
    par_temp.laser_lambda1 = wavelength;
    par_temp.pulsepow      = 1;  % [mW cm^-2]

    par_temp.gx1 = generation(par_temp, 'laser', wavelength);

    soleq_temp     = soleq.el;
    soleq_temp.par = par_temp;

    Vmin     = -0.3;
    Vmax     =  0.6;
    scan_rate = 50e-3;
    cycles    = 1;
    tpoints   = 400;
    CVsol = doCV(soleq_temp, 1, 0, Vmax, Vmin, scan_rate, cycles, tpoints);

    V = dfana.calcVapp(CVsol);
    J = dfana.calcJ(CVsol, "sub").tot;
    
    idx_maxV = find(V == max(V), 1);
    J_f = J(1:idx_maxV);
    V_f = V(1:idx_maxV);
    
    Jsc = abs(interp1(V_f, J_f, 0, 'linear'));  % [A cm^-2]

    Pin = par_temp.pulsepow * 1e-3;            % [W cm^-2]

    I = beerlambertI(par_temp, par_temp.xx, par_temp.light_source1, wavelength, 0);
    Plost = I(1, wavelength-300);

    EQE_calc = 100 * Jsc * h * c / (Pin * wavelength * 1e-9 * e);
    IQE_calc = EQE_calc * Pin / (Pin - Plost);
end