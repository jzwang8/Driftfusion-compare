function plot_caughey_thomas()
%PLOT_CAUGHEY_THOMAS  Plot electron and hole mobility vs doping concentration
%   using the Caughey-Thomas empirical model for Si at 300K.
%
%   Parameters from Sze, Physics of Semiconductor Devices (also Wikipedia):
%     mu_n = 65 + 1265 / (1 + (N/8.5e16)^0.72)
%     mu_p = 48 + 447  / (1 + (N/6.3e16)^0.76)

N = logspace(14, 20, 500);

% Caughey-Thomas parameters for Si at 300K (Sze / Wikipedia)
[mu_n, mu_p] = getMobilitiesFromDopingConc(N);

figure('Name', 'Caughey-Thomas', 'Position', [100 100 700 500], 'Color', 'w');
loglog(N, mu_n, 'b-', 'LineWidth', 2, 'DisplayName', 'Electrons');
hold on;
loglog(N, mu_p, 'r-', 'LineWidth', 2, 'DisplayName', 'Holes');
hold off;

xlabel('Doping concentration N [cm^{-3}]', 'FontSize', 14);
ylabel('Mobility \mu [cm^2/Vs]', 'FontSize', 14);
% title('Caughey--Thomas mobility vs. doping (Sze)', 'FontSize', 15);
legend('Location', 'northeast', 'FontSize', 12);
set(gca, 'FontSize', 13, 'LineWidth', 1.2, 'Box', 'on');
xlim([1e14, 1e20]);
ylim([30, 3000]);
grid on;

end