% Plot EQE of all runs on top of each other

if ~exist('runs','var')
    error('Variable ''runs'' not found in workspace.');
end

figure;
hold on;

for k = 1:numel(runs)
    % Basic safety checks
    if ~isfield(runs(k),'EQE_wavelengths') || ~isfield(runs(k),'EQE')
        warning('Run %d does not have EQE data – skipping.', k);
        continue;
    end

    wl  = runs(k).EQE_wavelengths;
    eqe = runs(k).EQE;

    % Fallback label if none stored
    if isfield(runs(k),'label')
        lab = runs(k).label;
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

% % Plot electric fields E(x) of all runs on top of each other
% 
% if ~exist('runs','var')
%     error('Variable ''runs'' not found in workspace.');
% end
% 
% figure;
% hold on;
% 
% for k = 1:numel(runs)
%     % Check that this run has an equilibrium solution
%     if ~isfield(runs(k),'solEq') || ~isfield(runs(k).solEq,'el')
%         warning('Run %d has no solEq.el – skipping.', k);
%         continue;
%     end
% 
%     solEq = runs(k).solEq.el;
% 
%     % Get x mesh
%     x = solEq.x;      % [cm]
% 
%     % Compute electric field vs x, t
%     % (If your Driftfusion version uses a different name, change this line)
%     F = dfana.calcFx(solEq);   % e.g. returns F(t,x)
% 
%     % Take electric field at final time
%     % (if F is [nt x nx], this picks last time slice)
%     if ndims(F) == 2
%         E = F(end,:);   % [V/cm]
%     else
%         % If your version returns something else, adjust here accordingly
%         error('Unexpected shape returned by dfana.calcFx.');
%     end
% 
%     % Label for legend
%     if isfield(runs(k),'label')
%         lab = runs(k).label;
%     else
%         lab = sprintf('Run %d', k);
%     end
% 
%     plot(x*1e7, E, 'DisplayName', lab);   % x*1e7 to convert cm -> nm (if you like)
% end
% 
% xlabel('Position x [nm]');
% ylabel('Electric field F(x) [V/cm]');
% title('Electric field profiles for all runs');
% legend('show','Location','best');
% grid on;
% hold off;