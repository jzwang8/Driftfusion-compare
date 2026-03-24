% figure;
% plot(thisRun.SR_EQE);
% hold on;
% plot(thisRun.EQE,'LineStyle','--');

% % Generate profile
% N0 = 8.06e19;
% jd = 200;
% N_sub = 50;
% 
% [EF, mu_n, mu_p, taup, ND, x] = generate_emitter_EF_profile(N0, jd, 'gaussian', N_sub);
% 
% % Check: x(1) should be near surface, ND(1) should be highest
% fprintf('x(1)=%.1f nm, ND(1)=%.2e  (should be surface, highest)\n', x(1), ND(1));
% fprintf('x(end)=%.1f nm, ND(end)=%.2e  (should be junction, lowest)\n', x(end), ND(end));
% 
% % Now check what gets written to CSV
% temp_csv = write_emitter_profile_csv('Input_files/pn_junction_nochargedlayer.csv', EF, mu_n, mu_p, taup, jd);
% 
% % Read back and check EF0 column order
% T = readtable(temp_csv);
% % Find active rows
% active_mask = strcmp(T.layer_type, 'active');
% EF0_in_csv = T.EF0(active_mask);
% 
% fprintf('\nIn CSV (left=junction, right=surface):\n');
% fprintf('  First active row  EF0=%.4f (should be LOW doping = more negative EF)\n', EF0_in_csv(1)); %todo this is wrong 
% fprintf('  Last active row   EF0=%.4f (should be HIGH doping = less negative EF)\n', EF0_in_csv(end));

%% Plot doping profiles for EF_peak = -4.342 runs
runs = summary.runs;
T = summary.table;

%% Filter: successful runs with EF_peak = -4.342 and profile data
tol = 0.01;
EF_target = -4.05;

success_mask = strcmp(T.status, 'success');
ef_mask = abs(T.EF_peak_eV - EF_target) < tol;
mask = success_mask & ef_mask;

% Further filter to one set of conditions (e.g., SF=0, SRV=1000, PhiOff=0)
mask = mask & (T.SF_flag == 0) & ...
       (T.SRV_cm_s == 1000) & ...
       (abs(T.Phi_R_offset_eV) < tol);

idx_found = find(mask);
fprintf('Found %d runs matching EF_peak=%.3f, SF=0, SRV=1000, PhiOff=0\n', ...
        numel(idx_found), EF_target);

%% Plot
profile_types = {'uniform', 'gaussian', 'exponential'};
profile_labels = {'Uniform', 'Gaussian', 'Exponential'};
profile_colors = [0.27 0.33 0.80;    % blue
                  0.90 0.60 0.00;    % orange
                  0.13 0.55 0.13];   % green

figure('Name', 'Doping_Profile_EF4342', 'Position', [100 100 700 500], 'Color', 'w');
hold on;

for j = 1:3
    prof_mask = mask & strcmp(T.profile_type_col, profile_types{j});
    idx = find(prof_mask, 1);
    if ~isempty(idx)
        r = runs(T.run_number(idx));
        if ~isempty(r.profile.x_nm) && ~isempty(r.profile.N_D)
            plot(r.profile.x_nm, r.profile.N_D, '-', ...
                 'LineWidth', 2.5, 'Color', profile_colors(j,:), ...
                 'DisplayName', profile_labels{j});
        else
            fprintf('  %s: no profile arrays stored\n', profile_types{j});
        end
    else
        fprintf('  %s: no matching run found\n', profile_types{j});
    end
end

xlabel('Depth from surface (nm)', 'FontSize', 20);
ylabel('N_D (cm^{-3})', 'FontSize', 20);
set(gca, 'YScale', 'log', 'FontSize', 16, 'LineWidth', 1.2, 'Box', 'on');
grid on;
legend('Location', 'best', 'FontSize', 14);
title(sprintf('Emitter Doping Profile (N_0 \\approx 8e19, EF_{peak} = %.3f eV)', EF_target), ...
      'FontSize', 16);
hold off;