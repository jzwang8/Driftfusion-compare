if ~exist('runs','var')
    error('Variable ''runs'' not found in workspace.');
end

runs_idx = [20,22,19,21];
nRuns    = numel(runs_idx);

% -------------------------------------------------------------------------
% Toggles to include / exclude measured cell EQE
% -------------------------------------------------------------------------
plot_meas_with = false;   % measured cell with Tc/ZnPc/Tc
plot_meas_no   = false;   % measured cell without Tc/ZnPc

%% ------------------------------------------------------------------------
% 1) Plot EQE of all runs on top of each other
% -------------------------------------------------------------------------
figure;
ax = gca;

% Give this axes nRuns distinct colors
% colororder(ax, turbo(nRuns));   % or lines(nRuns), parula(nRuns), hsv(nRuns), ...

hold(ax, 'on');

for k = runs_idx
    if ~isfield(runs(k),'EQE_wavelengths') || ~isfield(runs(k),'EQE')
        warning('Run %d does not have EQE data – skipping EQE.', k);
        continue;
    end

    wl  = runs(k).EQE_wavelengths;
    eqe = runs(k).EQE;

    if isfield(runs(k),'label')
        lab = runs(k).label;
    else
        lab = sprintf('Run %d', k);
    end

    % No 'Color' argument – MATLAB will take the next color from ColorOrder
    plot(ax, wl, eqe, 'DisplayName', lab);
end

% --- measured cell with Tc/ZnPc/Tc ---
if plot_meas_with
    meas_with = readmatrix('250718_planar_4_DIW5min_ZnPc_Tc.txt');   % [λ(nm), EQE(%), ...]
    wl_with   = meas_with(:,1);
    eqe_with  = meas_with(:,2);

    plot(ax, wl_with, eqe_with, 'k', 'LineWidth', 2, ...
         'DisplayName', 'Tc/ZnPc/n^+-p MW');
end

% --- measured cell without Tc/ZnPc ---
if plot_meas_no
    meas_no = readmatrix('250717_diw_5min_4_modified.txt');   % <-- put your filename here
    wl_no   = meas_no(:,1);
    eqe_no  = meas_no(:,2);

    plot(ax, wl_no, eqe_no, 'Color', [0.4 0.4 0.4], 'LineWidth', 2, ...
         'LineStyle', '--', ...
         'DisplayName', 'n^+-p MW');
end

xlabel('Wavelength (nm)');
xlim([400,800])
ylabel('EQE (%)');
% title('EQE spectra for all runs');
legend('show','Location','best');
grid on;
hold(ax, 'off');
% %% ------------------------------------------------------------------------
% % 2) Plot electric field profiles for all runs
% % -------------------------------------------------------------------------
% figure;
% hold on;
% 
% for ii = 1:numel(runs_idx)
%     k = runs_idx(ii);
% 
%     if ~isfield(runs(k),'solEq') || ~isfield(runs(k).solEq,'el')
%         warning('Run %d has no solEq.el – skipping E-field.', k);
%         continue;
%     end
% 
%     solEq_el = runs(k).solEq.el;
% 
%     [u, t, x, par, dev, n, p, a, c, V] = dfana.splitsol(solEq_el); %#ok<ASGLU>
% 
%     x_vec = x(:).';      
% 
%     if ndims(V) == 2
%         V_end = V(end,:);
%     else
%         warning('Run %d: unexpected V size in splitsol – skipping E-field.', k);
%         continue;
%     end
% 
%     L      = min(numel(x_vec), numel(V_end));
%     x_plot = x_vec(1:L);
%     V_plot = V_end(1:L);
% 
%     dx = diff(x_plot);
%     dV = diff(V_plot);
%     E  = -dV ./ dx;
% 
%     x_mid = x_plot(1:end-1) + dx/2;
% 
%     if isfield(runs(k),'label')
%         lab = runs(k).label;
%     else
%         lab = sprintf('Run %d', k);
%     end
% 
%     plot(x_mid*1e7, E, 'DisplayName', lab, 'Color', colors(ii,:));   % x in nm
% end
% 
% xlabel('Position x [nm]');
% ylabel('Electric field E(x) [V/cm]');
% title('Electric field profiles (final time) for all runs');
% legend('show','Location','best');
% grid on;
% hold off;
% 
% %% ------------------------------------------------------------------------
% % 3) Plot recombination rate profiles for all runs
% % -------------------------------------------------------------------------
% figure;
% hold on;
% 
% for ii = 1:numel(runs_idx)
%     k = runs_idx(ii);
% 
%     if ~isfield(runs(k),'solCV')
%         warning('Run %d has no solCV – skipping recombination.', k);
%         continue;
%     end
% 
%     solCV = runs(k).solCV;
% 
%     [u, t, x, par, dev, n, p, a, c, V] = dfana.splitsol(solCV); %#ok<ASGLU>
% 
%     if isfield(par, 'x_sub')
%         x_sub = par.x_sub;   % [cm]
%     else
%         x_sub = x;
%     end
%     x_sub = x_sub(:).';
% 
%     r = dfana.calcr(solCV, "sub");
% 
%     fields = intersect({'btb','srh','vsr','tot'}, fieldnames(r));
%     if isempty(fields)
%         warning('Run %d: calcr returned no recognised fields – skipping.', k);
%         continue;
%     end
% 
%     if isfield(r,'tot')
%         r_final = r.tot(end,:);
%     else
%         r_sum = zeros(1, size(r.(fields{1}),2));
%         for f = 1:numel(fields)
%             r_sum = r_sum + r.(fields{f})(end,:);
%         end
%         r_final = r_sum;
%     end
% 
%     r_final = r_final(:).';
%     L      = min(numel(x_sub), numel(r_final));
%     x_plot = x_sub(1:L);
%     r_plot = r_final(1:L);
% 
%     if isfield(runs(k),'label')
%         lab = runs(k).label;
%     else
%         lab = sprintf('Run %d', k);
%     end
% 
%     semilogy(x_plot*1e7, r_plot, 'DisplayName', lab, 'Color', colors(ii,:));
% end
% 
% xlabel('Position x [nm]');
% ylabel('Recombination rate r(x) [cm^{-3} s^{-1}]');
% title('Recombination profiles (final time) for all runs');
% legend('show','Location','best');
% grid on;
% hold off;