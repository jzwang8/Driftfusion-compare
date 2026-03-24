function temp_csv_path = write_emitter_profile_csv(base_csv, EF_profile, mu_n_profile, mu_p_profile, taup_profile, junction_depth_nm)
%WRITE_EMITTER_PROFILE_CSV Generate a temporary CSV with expanded emitter sublayers
%
%   Reads the base CSV, replaces the n-emitter row with N_sublayers rows
%   (one per profile sublayer), writes to a temp file, and returns the path.
%
%   The temp file is named after the parameters to avoid race conditions
%   when running in parallel (each parfor worker gets its own file).
%
%   Inputs:
%     base_csv          - path to original CSV (e.g. 'Input_files/pn_junction_nochargedlayer.csv')
%     EF_profile        - [1 x N_sub] Fermi levels per sublayer [eV]
%     mu_n_profile      - [1 x N_sub] electron mobilities [cm^2/Vs]
%     mu_p_profile      - [1 x N_sub] hole mobilities [cm^2/Vs]
%     taup_profile      - [1 x N_sub] hole lifetimes [s]
%     junction_depth_nm - total emitter thickness [nm]
%
%   Output:
%     temp_csv_path - path to the written temp CSV file

    N_sub = length(EF_profile);
    sublayer_thickness_cm = (junction_depth_nm * 1e-7) / N_sub;  % nm -> cm

    %% Read base CSV as raw text (preserves formatting)
    fid = fopen(base_csv, 'r');
    lines = {};
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line)  % fgetl returns -1 (numeric) at EOF - skip it
            line = strrep(line, sprintf('\r'), '');  % strip \r from Windows line endings
            lines{end+1} = line; %#ok<AGROW>
        end
    end
    fclose(fid);

    %% Find the n-emitter row (last 'active' row before the right electrode)
    % Strategy: find the row with the highest (least negative) EF0 among active rows
    % In this CSV: row 5 (1-indexed after header) is the n-emitter
    % More robustly: find the last 'active' row before the right 'electrode' row

    header = lines{1};
    header = regexprep(header, '^\xEF\xBB\xBF', '');  % strip UTF-8 BOM if present
    header = regexprep(header, '^\xFF\xFE', '');        % strip UTF-16 BOM if present
    lines{1} = header;
    col_names = regexp(header, ',', 'split');

    % Find column indices we need
    col_layer_type = find(strcmp(col_names, 'layer_type'));
    col_thickness  = find(strcmp(col_names, 'thickness'));
    col_points     = find(strcmp(col_names, 'layer_points'));
    col_EF0        = find(strcmp(col_names, 'EF0'));
    col_mu_n       = find(strcmp(col_names, 'mu_n'));
    col_mu_p       = find(strcmp(col_names, 'mu_p'));
    col_taup       = find(strcmp(col_names, 'taup'));

    % Sanity check - fail clearly if columns not found
    required = {'col_layer_type', 'col_thickness', 'col_points', 'col_EF0', 'col_mu_n', 'col_mu_p', 'col_taup'};
    required_vals = {col_layer_type, col_thickness, col_points, col_EF0, col_mu_n, col_mu_p, col_taup};
    for ri = 1:length(required)
        if isempty(required_vals{ri})
            error('write_emitter_profile_csv: column "%s" not found in CSV header.\nHeader: %s', ...
                  required{ri}, header);
        end
    end
    emitter_line_idx = [];
    for i = length(lines):-1:2
        fields = strsplit(lines{i}, ',');
        if length(fields) >= col_layer_type
            lt = strtrim(fields{col_layer_type});
            if strcmpi(lt, 'active')
                emitter_line_idx = i;
                break;
            end
        end
    end

    if isempty(emitter_line_idx)
        error('write_emitter_profile_csv: could not find n-emitter (last active) row in CSV');
    end

    %% Parse emitter row into fields
    emitter_fields = regexp(lines{emitter_line_idx}, ',', 'split');

    % Points per sublayer: distribute original total evenly
    orig_points = str2double(strtrim(emitter_fields{col_points}));
    points_per_sub = max(6, round(orig_points / N_sub));  % at least 6 per sublayer

    %% Flip profiles: x_sublayers convention is surface->junction,
    %  but CSV sublayer order is junction->surface (left to right)
    EF_profile   = fliplr(EF_profile);
    mu_n_profile = fliplr(mu_n_profile);
    mu_p_profile = fliplr(mu_p_profile);
    taup_profile = fliplr(taup_profile);

    %% Build replacement rows (one per sublayer)
    sublayer_lines = cell(N_sub, 1);
    for i = 1:N_sub
        f = emitter_fields;  % copy all fields from original emitter row

        % Override doping-dependent fields
        f{col_thickness} = sprintf('%.4E', sublayer_thickness_cm);
        f{col_points}    = sprintf('%d', points_per_sub);
        f{col_EF0}       = sprintf('%.6f', EF_profile(i));
        f{col_mu_n}      = sprintf('%.4E', mu_n_profile(i));
        f{col_mu_p}      = sprintf('%.4E', mu_p_profile(i));
        f{col_taup}      = sprintf('%.4E', taup_profile(i));

        sublayer_lines{i} = strjoin(f, ',');
    end

    %% Assemble new CSV lines
    new_lines = [lines(1:emitter_line_idx-1), ...  % header + layers before emitter
                 sublayer_lines',               ...  % expanded emitter sublayers
                 lines(emitter_line_idx+1:end)];     % right electrode

    %% Write to temp file (unique name to avoid parallel conflicts)
    hash = num2str(sum(EF_profile * (1:N_sub)'), 'hash%+.0f');
    hash = strrep(hash, '-', 'n');
    hash = strrep(hash, '.', 'p');
    hash = strrep(hash, '+', '');
    temp_csv_path = fullfile(tempdir, sprintf('df_emitter_%dsub_%s.csv', N_sub, hash));

    fid = fopen(temp_csv_path, 'w');
    for i = 1:length(new_lines)
        fprintf(fid, '%s\n', new_lines{i});
    end
    fclose(fid);

    fprintf('  [csv] Wrote %d-sublayer emitter CSV: %s\n', N_sub, temp_csv_path);
end