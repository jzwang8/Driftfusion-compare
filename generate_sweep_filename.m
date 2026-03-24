function filename = generate_sweep_filename(SF_flag, SF_efficiency, junction_depth, sp_r, N0_peak, Phi_R_offset, profile_type)
%GENERATE_SWEEP_FILENAME Create deterministic filename from parameters
%
%   Shared by run_sweep_chunk_profiles.m and postprocess_sweep_profiles.m
%   to guarantee filenames always match.
%
%   Format: SF{0/1}_eta{X}_jd{X}_sp{X}_N0{X}_PhiOff{X}_prof{X}.mat

    sf_str   = sprintf('SF%d', SF_flag);
    eta_str  = sprintf('eta%.1f', SF_efficiency);
    jd_str   = sprintf('jd%d', junction_depth);
    sp_str   = sprintf('sp%.0e', sp_r);
    n0_str   = sprintf('N0%.0e', N0_peak);
    phi_str  = sprintf('PhiOff%.4f', Phi_R_offset);
    prof_str = sprintf('prof%s', profile_type);

    sp_str  = strrep(sp_str,  '+', '');
    n0_str  = strrep(n0_str,  '+', '');
    phi_str = strrep(phi_str, '-', 'n');
    phi_str = strrep(phi_str, '.', 'p');

    filename = sprintf('%s_%s_%s_%s_%s_%s_%s.mat', ...
                       sf_str, eta_str, jd_str, sp_str, n0_str, phi_str, prof_str);
end
