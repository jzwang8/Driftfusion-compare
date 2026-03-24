function config = sweep_config_profiles()
%SWEEP_CONFIG_PROFILES Single source of truth for profile sweep parameters
%
%   config = sweep_config_profiles()
%
%   Returns a struct with all sweep parameter definitions and fixed settings.
%   Used by both run_sweep_chunk_profiles.m and postprocess_sweep_profiles.m
%   to guarantee identical parameter spaces.
%
%   MODIFY THESE TO ADD NEW VALUES - existing runs won't be repeated.

    %% Sweep parameter values
    config.junction_depth_vals = [200, 300, 400];       % nm
    config.sp_r_vals           = [10, 100, 1000, 10000, 100000];  % cm/s
    config.eta_sf_with_SF      = [0, 1, 2];            % SF efficiency when SF_flag=1

    % EF_vals represents PEAK doping concentration
    % (surface concentration for exponential; center peak for Gaussian)
    config.EF_vals = [-4.05, -4.104, -4.223, -4.342];  % eV (~8e19, 1e19, 1e17, 1e15 cm^-3)

    % Phi_R offset relative to the SURFACE EF (EF of outermost sublayer)
    config.Phi_R_offset_vals = [0.054, 0, -0.054];   % eV

    % Physical constraint: Phi_R cannot exceed conduction band edge
    config.E_C = -4.05;  % Si conduction band edge (eV)

    % Doping profile shapes
    % 'uniform'     - constant doping (matches previous sweep behavior)
    % 'gaussian'    - ion implantation profile, peak at x0=d/2, sigma=d/5
    % 'exponential' - diffusion doping, dose-matched to Gaussian
    config.profile_shapes = {'uniform', 'gaussian', 'exponential'};

    % Number of emitter sublayers for profile discretization
    % (uniform uses 1 layer; gaussian/exponential use N_sublayers)
    config.N_sublayers = 50;

    %% Reflection model
    %  'spectral' - wavelength-dependent R(lambda) from transfer matrix model
    %               (get_TMM_reflection.m); uses Tc/ZnPc/AlOx/Si stack when
    %               SF_flag=true, bare AlOx/Si when SF_flag=false.
    %  'constant' - flat reflectance = config.reflection at all wavelengths.
    config.reflection_model = 'spectral';   % 'spectral' | 'constant'
    config.reflection       = 0.15;         % constant kappa used for:
                                            %   - JV (par.kappa, always)
                                            %   - EQE (only when reflection_model='constant')

    %% Spectral response under bias settings
    %  Set to [] to disable reverse-bias sweeps (EQE at V=0 only).
    %  Add values like [-0.5, -1.0] to re-enable.
    config.SR_bias_vals = [];   % V (applied external bias)

    %% Fixed parameters
    % JV settings
    config.light_intensity_JV = 1;
    config.Vmin_JV   = -0.3;
    config.Vmax_JV   = 0.7;
    config.scan_rate_JV = 50e-3;
    config.cycles_JV = 2;
    config.tpoints_JV = 300;

    % EQE settings
    config.scan_range_1 = 405:5:600;
    config.scan_range_2 = 610:10:850;
    config.number_of_workers = 8;

    % Output directory
    config.outputDir = 'sweep_results_profiles';

    % Chunking
    config.runs_per_chunk = 20;

    % Input file
    config.file_name = 'Input_files/pn_junction_nochargedlayer.csv';
end