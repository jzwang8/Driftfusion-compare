%Sweep additional different hole flux densities and calculate, plot JV
%curves

clear all;
% j_array = zeros(1,10);
% n_array = 1:5;
hold on
holes_array = [0,1e3,1e6,1e9];
for n = 1:length(holes_array)
    initialise_df;
    par = pc('input_files/pn_junction_var6');
    par.extra_holes = holes_array(n);
    par.extra_electrons = 0;
    par.Vcharging = 0;
    soleq = equilibrate(par,1);
    
    txt = ['Injected hole flux = ',num2str(holes_array(n))];
    txt2 = ['Injected electron flux = ',num2str(par.extra_electrons)];
    dfplot.npx(soleq.el)
    % JVsol = doJV(soleq.el, 1e-2, 100, 1, 1, -0.5, .5, 3);
    %  j_array(n) = j.p(end,200);

%% Plot hole flux density as a function of distance:
    figure
    hold on
    [J, j, x] = dfana.calcJ(soleq.el, "sub");
    yyaxis left
    semilogy(x,j.p(end,:),'DisplayName',txt)
    ylabel('Hole flux density (cm^{-2}s^{-1})')
    yyaxis right
    semilogy(x,j.n(end,:),'DisplayName',txt2)
    xlabel('Position (cm)')
    ylabel('Electron flux density (cm^{-2}s^{-1})')
    legend show
    hold off
%% OLD: Plot JV curve (adapted from dfplot)
%     J.ill.f = dfana.calcJ(JVsol.ill.f, "sub");
%     Vapp.ill.f = dfana.calcVapp(JVsol.ill.f);
%     J.ill.r = dfana.calcJ(JVsol.ill.r, "sub");
%     Vapp.ill.r = dfana.calcVapp(JVsol.ill.r);
% 
%     txt = ['injected holes = 10^{',num2str(n), '}'];
%     % plot(Vapp.ill.f, J.ill.f.tot(:,end),'--')%, 'Color', [0, 0.4470, 0.7410]);
% 
%     plot(Vapp.ill.r, J.ill.r.tot(:,end),'DisplayName',txt);%,'Color', [0, 0.4470, 0.7410]);
% xlabel('Applied voltage [V]')
% ylabel('Current density [Acm-2]');
%     legend show

    %% Plot JV curve in the middle of the active area:
    % Parameters for the sweep:
    figure
    Vmin = -0.5;
    Vmax = 1;
    scan_rate = 50e-3;
    cycles = 2;
    tpoints = 300;
    % Perform the sweep:
    CVsol = doCV(soleq.el, 1, 0, Vmax, Vmin, scan_rate, cycles, tpoints);
    % Get the x position to plot:
    xmesh = CVsol.x;
    ppos = getpointpos(0, xmesh); %getpointpos(par.d_midactive, xmesh);

    J = dfana.calcJ(CVsol, "sub");
    Vapp = dfana.calcVapp(CVsol);
    % Plot:

    plot(Vapp, J.tot(:, ppos), 'DisplayName',txt);
    xlabel('Applied Voltage, Vapp [V]');
    ylabel('Current Density, J [A cm^{-2}]');
    set(legend,'FontSize',16);
    set(legend,'EdgeColor',[1 1 1]);
    legend show

    % Calculate PCE:
    power = Vapp.*transpose(J.tot(:, ppos));
    MPP = min(power);
    Pin = dfana.calcPin(CVsol);
    efficiency(n) = abs(100*(MPP/Pin));
end

hold off
