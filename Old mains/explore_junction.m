%% Example script: 1-parameter exploration and plotting with explore class
% Demonstrates how to run explore.explore1par and then use all
% plotting functions that are usable from the resulting exsol.
%
% LICENSE
% (same as original Driftfusion license…)

clear; clc;

%% Initialise driftfusion
initialise_df;   % make sure this is on your path

%% Base parameters
par_explore = pc('input_files/pn_junction_jzw.csv');
idx_active  = par_explore.active_layer; %#ok<NASGU>

%% Run 1-parameter exploration over active layer thickness d(1,3)

thickness_vec = [2e-5 3e-5 4e-5];   % cm (200, 300, 400 nm)

parex_dactive = explore.explore1par( ...
    par_explore, ...          % base par
    {'d(1,3)'}, ...           % parameter name
    {thickness_vec}, ...      % values
    200 );                    % JV points

%% Convenience logicals for plots
par1logical = true(1, numel(thickness_vec));  % plot all d(1,3)
par2logical = true;                           % only one dummy parval2 for explore1par

%% 1) JV curves (forward) for each thickness
figure;
explore.plotJV(parex_dactive, par1logical, par2logical);
legend('d_{active} = 200 nm', ...
       'd_{active} = 300 nm', ...
       'd_{active} = 400 nm', ...
       'Location','best');
title('JV (forward scan) vs active layer thickness');

%% 2) PL surface (if PL data available)
% For explore1par this will just be a "strip" in the second param direction.
try
    explore.plotPL(parex_dactive);
    title('PL intensity vs parameter(s)');
catch ME
    warning('explore.plotPL could not be plotted: %s', ME.message);
end

%% 3) Surfaces of key JV statistics vs parameter(s)
% These come from exsol.stats and work with the dummy second parameter.

% Voc_reverse surface
figure;
explore.plotsurf(parex_dactive, 'Voc_r', 0, 0, 0);
title('V_{oc} (reverse scan) vs parameters');

% Voc_forward surface
figure;
explore.plotsurf(parex_dactive, 'Voc_f', 0, 0, 0);
title('V_{oc} (forward scan) vs parameters');

% Jsc_forward surface
figure;
explore.plotsurf(parex_dactive, 'Jsc_f', 0, 0, 0);
title('J_{sc} (forward scan) vs parameters');

% Fill factor (forward)
figure;
explore.plotsurf(parex_dactive, 'FF_f', 0, 0, 0);
title('Fill factor (forward) vs parameters');

% mpp (forward)
figure;
explore.plotsurf(parex_dactive, 'mpp_f', 0, 0, 0);
title('Maximum power point (forward) vs parameters');

%% 4) 2D plots vs parval2 (x-axis) and parval1 (families)
% For explore1par, parval2 is just a single dummy point, so these are
% not super informative but they *do* run and are consistent with the
% explore API.

% y vs parval2, families of parval1
explore.plotstat_2D_parval1(parex_dactive, 'Voc_r', 0, 0);
title('V_{oc,r} vs parval2 (families over parval1)');

% y vs parval1, families of parval2
explore.plotstat_2D_parval2(parex_dactive, 'Voc_r', 0, 0);
title('V_{oc,r} vs parval1 (families over parval2)');

%% 5) Profile and recombination plots (only if profile fields exist)
% These require exsol.n_f, exsol.p_f, exsol.x, etc.
% In many Driftfusion builds they are commented out in explore2par /
% explore1par. We try them, but catch errors gracefully.

% 5a) Final energy level diagrams
try
    explore.plotfinalELx(parex_dactive);
    title('Final energy level diagrams');
catch ME
    warning('explore.plotfinalELx could not be plotted: %s', ME.message);
end

% 5b) 1D profiles of a chosen quantity (e.g. ''n_f'' or ''p_f'')
%     Here we just try ''n_f'' as an example, if it exists.
if isfield(parex_dactive, 'n_f')
    try
        figure;
        explore.plotprof_2D(parex_dactive, 'n_f', par1logical, par2logical, 0, 1);
        title('Electron density profiles');
    catch ME
        warning('explore.plotprof_2D(n_f) failed: %s', ME.message);
    end
end

% 5c) Recombination profiles and integrated recombination
try
    explore.plotU(parex_dactive, par1logical, par2logical, 0, 1);
    % plotU creates two figures: recombination vs x, and integrated vs thickness
catch ME
    warning('explore.plotU could not be plotted: %s', ME.message);
end

%% 6) CE (charge extraction) plots – require a second exsol (equilibrium)
% plotCE(exsol_Voc, exsol_eq, ...) needs:
%   - exsol_Voc : exploration at Voc (light on)
%   - exsol_eq  : exploration at equilibrium (dark)
% If you have those, you would do something like:
%
%   exsol_Voc = parex_dactive_Voc;   % from a separate explore run at Voc
%   exsol_eq  = parex_dactive_eq;    % from explore run for equilibrium
%   explore.plotCE(exsol_Voc, exsol_eq, 1, 0, 0, false);
%
% Since in this simple example we only have parex_dactive, we leave plotCE
% commented as a template.

%{
% Example template (requires separate exsol_Voc and exsol_eq):
explore.plotCE(exsol_Voc, exsol_eq, 1, 0, 0, false);
%}

%% Done
disp('Finished parameter exploration and plotting for d(1,3).');