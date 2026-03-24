%% Overlay EQE, E-field and recombination profiles for all runs
% Assumes 'runs' struct array is already in the base workspace,
% each element having: runs(k).solEq, runs(k).solCV, runs(k).label, etc.

if ~exist('runs_simcomp','var')
    error('Variable ''runs'' not found in workspace.');
end


runs_idx = [1,2,3,4,5,6];
%% ------------------------------------------------------------------------
% 1) Plot EQE of all runs on top of each other
% -------------------------------------------------------------------------
figure;
hold on;

for k = runs_idx
    if ~isfield(runs_simcomp(k),'EQE_wavelengths') || ~isfield(runs_simcomp(k),'EQE')
        warning('Run %d does not have EQE data – skipping EQE.', k);
        continue;
    end

    wl  = runs_simcomp(k).EQE_wavelengths;
    eqe = runs_simcomp(k).EQE;

    if isfield(runs_simcomp(k),'label')
        lab = runs_simcomp(k).label;
    else
        lab = sprintf('Run %d', k);
    end

    plot(wl, eqe, 'DisplayName', lab);
end

xlabel('Wavelength (nm)');
ylabel('EQE (%)');
title('EQE spectra for all runs');
legend('show','Location','best');
grid on;
hold off;

%% ------------------------------------------------------------------------
% 2) Plot electric field profiles for all runs
%    E(x) is computed from potential V(x) using dfana.splitsol
% -------------------------------------------------------------------------
figure;
hold on;

for k = runs_idx

    if ~isfield(runs_simcomp(k),'solEq') || ~isfield(runs_simcomp(k).solEq,'el')
        warning('Run %d has no solEq.el – skipping E-field.', k);
        continue;
    end

    solEq_el = runs_simcomp(k).solEq.el;

    % Split solution (same helper used inside dfplot functions)
    [u, t, x, par, dev, n, p, a, c, V] = dfana.splitsol(solEq_el); %#ok<ASGLU>

    % Use the full x mesh from splitsol
    x_vec = x(:).';      % [cm], row vector

    % Potential at final time step
    if ndims(V) == 2
        V_end = V(end,:);        % 1 x Nx
    else
        warning('Run %d: unexpected V size in splitsol – skipping E-field.', k);
        continue;
    end

    % Make sure same length
    L = min(numel(x_vec), numel(V_end));
    x_plot = x_vec(1:L);
    V_plot = V_end(1:L);

    % Electric field E = -dV/dx (finite differences)
    % V in [V], x in [cm] -> E in [V/cm]
    dx = diff(x_plot);
    dV = diff(V_plot);
    E  = -dV ./ dx;              % length L-1

    % For plotting, align E with midpoints of x
    x_mid = x_plot(1:end-1) + dx/2;

    if isfield(runs_simcomp(k),'label')
        lab = runs_simcomp(k).label;
    else
        lab = sprintf('Run %d', k);
    end

    plot(x_mid*1e7, E, 'DisplayName', lab);   % x in nm
end

xlabel('Position x [nm]');
ylabel('Electric field E(x) [V/cm]');
title('Electric field profiles (final time) for all runs');
legend('show','Location','best');
grid on;
hold off;

%% ------------------------------------------------------------------------
% 3) Plot recombination rate profiles for all runs
%    Uses the same approach as dfplot.rx:
%    r = dfana.calcr(solCV, "sub");
% -------------------------------------------------------------------------
figure;
hold on;

for k = runs_idx

    if ~isfield(runs_simcomp(k),'solCV')
        warning('Run %d has no solCV – skipping recombination.', k);
        continue;
    end

    solCV = runs_simcomp(k).solCV;

    % Split solution to get subdomain positions
    [u, t, x, par, dev, n, p, a, c, V] = dfana.splitsol(solCV); %#ok<ASGLU>
    % par.x_sub is what dfplot.rx uses
    if isfield(par, 'x_sub')
        x_sub = par.x_sub;   % [cm]
    else
        % fallback: use full x if x_sub not present
        x_sub = x;
    end
    x_sub = x_sub(:).';      % row vector

    % Recombination rates on substrate ("sub") domain
    r = dfana.calcr(solCV, "sub");   % struct: r.btb, r.srh, r.vsr, (maybe) r.tot

    % Each r.field is (nt × Nx_sub). We want the final time slice.
    fields = intersect({'btb','srh','vsr','tot'}, fieldnames(r));
    if isempty(fields)
        warning('Run %d: calcr returned no recognised fields – skipping.', k);
        continue;
    end

    % If r.tot exists, use it; otherwise sum the mechanisms.
    if isfield(r,'tot')
        r_final = r.tot(end,:);   % 1 × Nx_sub
    else
        r_sum = zeros(1, size(r.(fields{1}),2));
        for f = 1:numel(fields)
            r_sum = r_sum + r.(fields{f})(end,:);
        end
        r_final = r_sum;          % 1 × Nx_sub
    end

    % Ensure matching lengths for x and r
    r_final = r_final(:).';
    L = min(numel(x_sub), numel(r_final));
    x_plot = x_sub(1:L);
    r_plot = r_final(1:L);

    % Label
    if isfield(runs_simcomp(k),'label')
        lab = runs_simcomp(k).label;
    else
        lab = sprintf('Run %d', k);
    end

    semilogy(x_plot*1e7, r_plot, 'DisplayName', lab);  % x in nm, log-y
end

xlabel('Position x [nm]');
ylabel('Recombination rate r(x) [cm^{-3} s^{-1}]');
title('Recombination profiles (final time) for all runs');
legend('show','Location','best');
grid on;
hold off;