clear all; close all;

%% ---------- DIFFUSION DOPING : exponential fit -----------------------
data = readmatrix("n_diffdope.txt");
x = data(:,1);
y = data(:,2);

fit_idx    = (x>=0)&(x<=300);
p          = polyfit(x(fit_idx), log(y(fit_idx)), 1);
B_fit      = p(1);
A_fit      = exp(p(2));

x_fit      = linspace(0,300,500);
y_exp_fit  = A_fit .* exp(B_fit .* x_fit);

figure('Color','w');
semilogy(x, y, 'ko', 'MarkerSize',6, 'DisplayName','Data');  hold on
semilogy(x_fit, y_exp_fit, 'r-', 'LineWidth',3, 'DisplayName','Exponential Fit');

% purple guide line WITHOUT legend entry
xline(300, 'Color',[0.5 0 0.7], 'LineWidth',3, 'HandleVisibility','off','LineStyle','--');

xlabel('x (nm)','FontSize',20)
ylabel('N_D (cm^{-3})','FontSize',20)
legend('Location','southwest','FontSize',18)
set(gca,'FontSize',18)
xlim([0 400])                        % plot only 0–400 nm
grid on

%% ---------- ION IMPLANTATION : Gaussian fit --------------------------
data = readmatrix("n_ionimplant.txt");
x = data(:,1);  y = data(:,2);

gaussEqn = 'N0*exp(-((x-x0)^2)/(2*sigma^2))';
gFit     = fit(x, y, gaussEqn, 'StartPoint',[max(y) mean(x) std(x)]);

x_fit    = linspace(0,400,500);
y_gauss  = gFit.N0 .* exp(-((x_fit - gFit.x0).^2) / (2*gFit.sigma^2));

figure('Color','w');
semilogy(x, y, 'ko', 'MarkerSize',6, 'DisplayName','Data'); hold on
semilogy(x_fit, y_gauss, 'r-', 'LineWidth',3, 'DisplayName','Gaussian Fit');

xline(300, 'Color',[0.5 0 0.7], 'LineWidth',3, 'HandleVisibility','off','LineStyle','--');

xlabel('x (nm)','FontSize',20)
ylabel('N_D (cm^{-3})','FontSize',20)
legend('Location','southwest','FontSize',18)
set(gca,'FontSize',18)
xlim([0 400])
grid on