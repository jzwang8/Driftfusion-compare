%% Combined JV + profiles + recombination + EQE script
delete(gcp('nocreate'));
system('taskkill /F /IM EXCEL.EXE');  % optional cleanup of Excel
set(gcf, 'WindowStyle', 'normal')
% system('caffeinate')

% Add ../tmm to the path relative to this file
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Transfer matrix modeling'));
TcZnPCSi_thicknesses;
lambda_tmm=lambda;
Reflection_tmm=Reflection;

% Input parameters file
file_name = 'Input_files/pn_junction_nochargedlayer.csv';

% --- Base settings (defaults) ---
SF_flag          = true;
SF_efficiency    = 2; % par.eta
junction_depth   = 300; % [nm]
sp_r             = 100000; % SRV [cm/s]
nemitter_EF0     = -4.103947060219521; % [eV]
rightelectrode   = -4.103947060219521; % [eV]

%% ------------------------------------------------------------------------
%  Choose which variable to sweep
% -------------------------------------------------------------------------
% Options (as implemented below): 'SF_flag', 'eta_sf', 'junction_depth', 
% 'sp_r', 'nemitter_EF0', 'rightelectrode
sweepParam  = 'sp_r'; % <--- change this
sweepValues = [1e5]; % <--- and this (example: sweep charges 0 and 1e19)

% Reflection spectrum for EQE (κ(λ)):
reflection = 0.15; % par.kappa (JV uses this constant value)
kappa_lambda = lambda_tmm; % from TMM
%kappa_vals = reflection * ones(size(lambda_tmm)); % constant reflection
kappa_vals = Reflection_tmm; % TMM reflection

% JV settings
light_intensity_JV = 1;    % suns
Vmin_JV            = -0.3;
Vmax_JV            =  0.7;
scan_rate_JV       = 50e-3;
cycles_JV          = 2;
tpoints_JV         = 300;

% EQE settings
scan_start        = 405;    % nm
scan_end          = 800;    % nm
scan_interval     = 5;      % nm
number_of_workers = 12;     % desired workers for parpool

%% ------------------------------------------------------------------------
%  Initialise Driftfusion once
% -------------------------------------------------------------------------
initialise_df;

%% ------------------------------------------------------------------------
%  MAIN SWEEP LOOP
% -------------------------------------------------------------------------
%% ------------------------------------------------------------------------
%  Load existing runs from MAT file (if it exists)
% -------------------------------------------------------------------------
runsFile = 'runs.mat';

if isfile(runsFile)
    S = load(runsFile,'runs');
    if isfield(S,'runs')
        runs = S.runs;
    else
        runs = struct([]);
    end
else
    runs = struct([]);
end

% also push into base workspace so the existing code below still works
assignin('base','runs',runs);

for iSweep = 1:numel(sweepValues)
    fprintf('\n===== Sweep %d / %d: %s = %g =====\n', ...
            iSweep, numel(sweepValues), sweepParam, sweepValues(iSweep));
    %% --------------------------------------------------------------------
    %  Build parameter set for THIS sweep value
    % ---------------------------------------------------------------------
    par = pc(file_name);
    par.AbsTol        = 5e-6;
    par.RelTol        = 5e-3;
    par.MaxStepFactor = 1;
    par.kappa         = reflection; % JV uses constant kappa

    % --- Override just the swept parameter ---
    switch sweepParam

        case 'SF_flag'
            SF_flag = logical(sweepValues(iSweep));
            par.Tetracene_TF = SF_flag;

        case 'eta_sf'
            SF_efficiency = sweepValues(iSweep);
            par.eta = SF_efficiency;
        
        case 'junction_depth'
            junction_depth = sweepValues(iSweep);
            par.d(1,3) = junction_depth*1e-9*1e2;% [nm]

        case 'sp_r'
            sp_r = sweepValues(iSweep);
            par.sp_r = sp_r;

        case 'nemitter_EF0'
            nemitter_EF0 = sweepValues(iSweep);
            par.EF0(3) = nemitter_EF0;
            N_d = getDopingConcFromEf(nemitter_EF0);

            % recalculate mobilities
            [mu_n, mu_p] = getMobilitiesFromDopingConc(N_d);
            par.mu_n(3) = mu_n;
            par.mu_p(3) = mu_p;

            % recalculate minority carrier lifetime
            if N_d > 1e18
                tau_p = getTaupFromDopingConc(N_d);
                par.taup(3) = tau_p;
            end

        case 'rightelectrode'
            rightelectrode = sweepValues(iSweep);
            par.Phi_right = rightelectrode;

        otherwise
            error('Unknown sweepParam "%s".', sweepParam);
    end

    % Refresh params of this device
    par = refresh_device(par);

    % Generation profile for JV run
    par.gx1 = generation(par, par.light_source1, par.laser_lambda1);
    
    % txt  = ['Injected hole flux = ', num2str(par.extra_holes)];
    % txt2 = ['Injected electron flux = ', num2str(par.extra_electrons)];

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
    xrange_JV = [1.795*1e5, 1.803*1e5];

    % Profiles
    dfplot.npx    (solEq.el, tarr_JV, xrange_JV);
    dfplot.rhoxFxVx(solEq.el, tarr_JV, xrange_JV);
    dfplot.ELx    (solEq.el, tarr_JV, xrange_JV);
    dfplot.rx     (solCV_JV, tarr_JV, xrange_JV);

    power_JV = Vapp_JV .* transpose(J_CV_JV.tot(:, ppos_JV)); % W cm^-2
    MPP_JV   = min(power_JV);
    Pin_JV   = dfana.calcPin(solCV_JV);
    eff_JV   = abs(100 * (MPP_JV / Pin_JV)); % [%]

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
    %  2) EQE / IQE vs wavelength
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
            EQE_parallel_loop(solEq, wavelength, SF_flag, ...
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

    thisRun.par       = par;

    % Parameters sweeping
    thisRun.SF_flag          = SF_flag;
    thisRun.SF_efficiency    = SF_efficiency;
    thisRun.junction_depth   = junction_depth;
    thisRun.sp_r             = sp_r;
    thisRun.nemitter_EF0     = nemitter_EF0;
    thisRun.rightelectrode   = rightelectrode;

    % JV info
    thisRun.Vapp      = Vapp_JV;
    thisRun.Jtot      = J_CV_JV.tot(:, ppos_JV);
    thisRun.solEq     = solEq;
    thisRun.solCV     = solCV_JV;
    thisRun.Voc       = Voc_JV;
    thisRun.PCE       = eff_JV;
    thisRun.MPP_power = MPP_JV;
    thisRun.Pin       = Pin_JV;

    % Human-readable identifiers
    thisRun.label = sprintf('%s = %g', sweepParam, sweepValues(iSweep));

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
    if evalin('base','exist(''runs'',''var'')')
        runs = evalin('base','runs');
    else
        runs = struct([]);   % 0x0 empty struct array
    end

    if isempty(runs)
        % First ever run, just start the array
        runs = thisRun;

    elseif ~isstruct(runs)
        % Something weird in workspace, reset
        warning('Variable "runs" exists but is not a struct; resetting it with thisRun.');
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
    assignin('base','runs',runs);
    fprintf('Stored run #%d in variable ''runs'' in workspace.\n', numel(runs));

    %% ------------------------------------------------------------------------
    %  Save updated runs to MAT file
    % -------------------------------------------------------------------------
    runs = evalin('base','runs');   % get runs from base
    save(runsFile,'runs','-v7.3');
    disp(['Saved ', num2str(numel(runs)), ' runs to ', runsFile]);
end


%% ------------------------------------------------------------------------
%  Local helper function
% -------------------------------------------------------------------------
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

% system('killall caffeinate');