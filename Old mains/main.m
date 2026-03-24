%% Input parameters:
% file_name = 'Input_files/No SRH';
file_name = 'Input_files/pn_junction_chargedlayer_nosrh.csv';

%% Calculate solutions:
% initialize system
initialise_df;
par = pc(file_name);
par.AbsTol = 5e-6; % solver parameters to tweak when script won't run (default = 1e-6)
par.RelTol = 5e-3; % solver parameters to tweak when script won't run (default = 1e-3)
par.MaxStepFactor = 1; % solver parameters to tweak when script won't run (default = 1)

par.eta = 2;
par.Tetracene_TF = true;
par.kappa = 0.15; % auto 0

%% Charged layer stuff (make into function)
layer_points_p_layer = par.parr(1);
% Get thickness of charged layer
charges = 1e19; % [cm-3]
% charges = 0; % [cm-3]
thickness = 1.5;  % [nm]
x = par.xx;
layer_index = find(x - (x(end) - thickness*1e-7) > 0, 1, 'first'); %finds index for the thickness of charged layer
actual_thickness = (x(end)-x(layer_index))*1e7;%compensates for the numerical grid that might not take exactly thickness as thickness
compensation_charges = -1 * actual_thickness/(x(layer_points_p_layer)*1e7) * charges;

par.dev_sub.Extra_charge(layer_index:end) = charges;
par.dev_sub.Extra_charge(1:layer_points_p_layer) = compensation_charges;

%% Initialize system:
% par.light_source1 = 'laser';
% par.laser_lambda1 = 525;
par.gx1 = generation(par, par.light_source1, par.laser_lambda1);
txt = ['Injected hole flux = ',num2str(par.extra_holes)];
txt2 = ['Injected electron flux = ',num2str(par.extra_electrons)];

% Calculate equilibrium solution
solEq = equilibrate(par);

% Calculate JV-sweep
light_intensity = 1; % sun
Vmin = -0.3;
Vmax = 0.7;
scan_rate = 50e-3;
cycles = 2;
tpoints = 300;
solCV = doCV(solEq.el, light_intensity, 0, Vmax, Vmin, scan_rate, cycles, tpoints);

% Get the x position to plot:
xmesh = solCV.x;
ppos = getpointpos(par.d_midactive, xmesh);

% Calculate current density from CV data
J_CV = dfana.calcJ(solCV, "sub");
Vapp = dfana.calcVapp(solCV);

% Plot solutions:
tarr = solEq.el.t(end);
xrange = [1.795*1e5,1.803*1e5]; % Standard input params
% xrange = [1.803006*1e5,1.803019*1e5]; % No SRH
% Plot charge carrier densities as a function of distance
% dfplot.npx(solEq.el,tarr,xrange) % FIGURE

% Plot charge density, electric field, electrostatic potential as a function of distance
% dfplot.rhoxFxVx(solEq.el,tarr,xrange) % FIGURE

% Plot energy level diagram as a function of distance
% dfplot.ELx(solEq.el,tarr,xrange) % FIGURE

% % Plot flux densities as a function of distance:
% figure;
% clf
% hold on
% [J, j, x] = dfana.calcJ(solEq.el, "sub");
% yyaxis left
% semilogy(x,j.p(end,:),'DisplayName',txt)
% ylabel('Hole flux density (cm^{-2}s^{-1})')
% yyaxis right
% semilogy(x,j.n(end,:),'DisplayName',txt2)
% xlabel('Position (cm)')
% ylabel('Electron flux density (cm^{-2}s^{-1})')
% legend show
% hold off

% Calculate PCE:
power = Vapp.*transpose(J_CV.tot(:, ppos)); %W cm-2
MPP = min(power);
Pin = dfana.calcPin(solCV); %W cm-2
efficiency = abs(100*(MPP/Pin));

% Plot JV curve in the middle of the active area
% FIGURE
figure
hold on
txt3 = ['PCE = ',num2str(efficiency), ' %'];
plot(Vapp, 1000*J_CV.tot(:, ppos), 'DisplayName',txt3);
xlabel('Applied Voltage, Vapp [V]');
ylabel('Current Density, J [mA cm^{-2}]');
set(legend,'FontSize',16);
set(legend,'EdgeColor',[1 1 1]);
legend show

% Interpolate to find Voc (V when J = 0)
Voc = interp1(J_CV.tot(:, ppos), Vapp, 0); % interpolates V when J = 0
disp(['Voc: ', num2str(Voc), ' V']);
disp(['PCE: ', num2str(efficiency), ' %']);

% Plot recombination rates as a function of position:
dfplot.rx(solCV,tarr,xrange) % FIGURE

% SPV
% spvsol = doSPV(solEq.el, 1, 0, 200, 1e-3, 1e6, 1);
% spvsol.dat = spvana(spvsol);
% delta_sigma = max(spvsol.dat.sigmat) - min(spvsol.dat.sigmat);
