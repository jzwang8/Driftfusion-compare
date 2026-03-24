function [combinations, combo_profile_names, EF_to_ND, n_skipped_invalid] = ...
        generate_sweep_combinations(config)
%GENERATE_SWEEP_COMBINATIONS Generate all valid parameter combinations
%
%   [combinations, combo_profile_names, EF_to_ND, n_skipped_invalid] = ...
%       generate_sweep_combinations(config)
%
%   Single source of truth for combination generation. Used by both
%   run_sweep_chunk_profiles.m and postprocess_sweep_profiles.m.
%
%   Requires the Doping folder on the MATLAB path (for getDopingConcFromEf
%   and getEfFromDopingConc). If those functions are unavailable, falls back
%   to Boltzmann approximations with a warning.
%
%   Output columns (10):
%     1:SF_flag  2:eta_sf  3:jd  4:sp  5:N0_peak  6:EF_peak
%     7:Phi_R    8:Phi_R_offset  9:profile_idx  10:EF_surface

    %% Unpack config
    junction_depth_vals = config.junction_depth_vals;
    sp_r_vals           = config.sp_r_vals;
    eta_sf_with_SF      = config.eta_sf_with_SF;
    EF_vals             = config.EF_vals;
    Phi_R_offset_vals   = config.Phi_R_offset_vals;
    E_C                 = config.E_C;
    profile_shapes      = config.profile_shapes;
    N_sublayers         = config.N_sublayers;

    %% Pre-compute N_D for each EF value
    use_exact = exist('getDopingConcFromEf', 'file') == 2;
    if ~use_exact
        warning('sweep:approx', ...
            ['getDopingConcFromEf not on path — using Boltzmann approximation.\n' ...
             'Add the Doping folder to your MATLAB path for exact values.']);
    end

    EF_to_ND = containers.Map('KeyType', 'double', 'ValueType', 'double');

    fprintf('\nFermi level (peak) -> Doping mapping:\n');
    for i = 1:length(EF_vals)
        EF = EF_vals(i);
        if use_exact
            [N_D, ~] = getDopingConcFromEf(EF);
        else
            N_D = boltzmann_ND_from_EF(EF);
        end
        EF_to_ND(EF) = N_D;
        fprintf('  E_F_peak = %.4f eV -> N0_peak = %.2e cm^-3\n', EF, N_D);
    end

    %% Choose EF-from-ND function
    if exist('getEfFromDopingConc', 'file') == 2
        ef_from_nd = @getEfFromDopingConc;
    else
        ef_from_nd = @boltzmann_EF_from_ND;
    end

    %% Generate all parameter combinations
    combinations = [];
    combo_profile_names = {};
    n_skipped_invalid = 0;

    % --- SF_flag = 0, eta_sf = 0 only ---
    for i_prof = 1:length(profile_shapes)
        profile_type = profile_shapes{i_prof};
        for i_jd = 1:length(junction_depth_vals)
            for i_sp = 1:length(sp_r_vals)
                for i_EF = 1:length(EF_vals)
                    for i_offset = 1:length(Phi_R_offset_vals)

                        junction_depth = junction_depth_vals(i_jd);
                        sp_r = sp_r_vals(i_sp);
                        EF_peak = EF_vals(i_EF);
                        N0_peak = EF_to_ND(EF_peak);
                        Phi_R_offset = Phi_R_offset_vals(i_offset);

                        EF_surface = get_surface_EF(N0_peak, junction_depth, ...
                            profile_type, N_sublayers, ef_from_nd);
                        Phi_R = EF_surface + Phi_R_offset;

                        if Phi_R > E_C
                            n_skipped_invalid = n_skipped_invalid + 1;
                            continue;
                        end

                        combinations = [combinations; ...
                            0, 0, junction_depth, sp_r, N0_peak, EF_peak, ...
                            Phi_R, Phi_R_offset, i_prof, EF_surface]; %#ok<AGROW>
                        combo_profile_names{end+1} = profile_type; %#ok<AGROW>
                    end
                end
            end
        end
    end

    % --- SF_flag = 1, all eta_sf values ---
    for i_eta = 1:length(eta_sf_with_SF)
        for i_prof = 1:length(profile_shapes)
            profile_type = profile_shapes{i_prof};
            for i_jd = 1:length(junction_depth_vals)
                for i_sp = 1:length(sp_r_vals)
                    for i_EF = 1:length(EF_vals)
                        for i_offset = 1:length(Phi_R_offset_vals)

                            eta_sf = eta_sf_with_SF(i_eta);
                            junction_depth = junction_depth_vals(i_jd);
                            sp_r = sp_r_vals(i_sp);
                            EF_peak = EF_vals(i_EF);
                            N0_peak = EF_to_ND(EF_peak);
                            Phi_R_offset = Phi_R_offset_vals(i_offset);

                            EF_surface = get_surface_EF(N0_peak, junction_depth, ...
                                profile_type, N_sublayers, ef_from_nd);
                            Phi_R = EF_surface + Phi_R_offset;

                            if Phi_R > E_C
                                n_skipped_invalid = n_skipped_invalid + 1;
                                continue;
                            end

                            combinations = [combinations; ...
                                1, eta_sf, junction_depth, sp_r, N0_peak, EF_peak, ...
                                Phi_R, Phi_R_offset, i_prof, EF_surface]; %#ok<AGROW>
                            combo_profile_names{end+1} = profile_type; %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end
end


%% ========================================================================
%  PRIVATE HELPER FUNCTIONS
%% ========================================================================

function EF_surface = get_surface_EF(N0_peak, junction_depth_nm, profile_type, N_sublayers, ef_from_nd_fn)
%GET_SURFACE_EF Return EF at the outermost (surface) sublayer for Phi_R reference

    d = junction_depth_nm;
    sublayer_thickness = d / N_sublayers;
    x_surface = sublayer_thickness / 2;  % center of first sublayer (nm)

    switch lower(profile_type)
        case 'uniform'
            N_surface = N0_peak;
        case 'gaussian'
            x0    = d / 2;
            sigma = d / 5;
            N_surface = N0_peak * exp(-((x_surface - x0)^2) / (2 * sigma^2));
            N_surface = max(N_surface, 1e15);
        case 'exponential'
            N_surface = N0_peak;
        otherwise
            N_surface = N0_peak;
    end

    EF_surface = ef_from_nd_fn(N_surface);
end


function EF = boltzmann_EF_from_ND(N_D)
%BOLTZMANN_EF_FROM_ND Boltzmann approximation for n-type Si Fermi level
%   Fallback when getEfFromDopingConc is not available.
%   Matches to within ~1 meV for the standard doping levels.

    E_C = -4.05;     % eV
    N_C = 2.86e19;   % cm^-3
    kT  = 0.02585;   % eV at 300 K

    EF = E_C + kT * log(N_D / N_C);
end


function N_D = boltzmann_ND_from_EF(EF)
%BOLTZMANN_ND_FROM_EF Inverse Boltzmann approximation
%   Fallback when getDopingConcFromEf is not available.

    E_C = -4.05;     % eV
    N_C = 2.86e19;   % cm^-3
    kT  = 0.02585;   % eV at 300 K

    N_D = N_C * exp((EF - E_C) / kT);
end
