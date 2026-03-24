%% Input parameters:
file_name = 'Input_files/pn_junction_jzw.csv';

%% Calculate solutions:
% initialize system
initialise_df;
par = pc(file_name);
par.AbsTol = 1e-6; % solver parameters to tweak when script won't run (default = 1e-6)
par.RelTol = 1e-3; % solver parameters to tweak when script won't run (default = 1e-3)
par.MaxStepFactor = 1; % solver parameters to tweak when script won't run (default = 1)

par.eta = 1;
par.Tetracene_TF = false;
par.kappa = 0; % auto 0

%% Charged layer stuff (make into function)
% layer_points_p_layer = par.parr(1);
% % Get thickness of charged layer
% % charges = 1e19; % [cm-3]
% charges = 0; % [cm-3]
% thickness = 1.5;  % [nm]
% x = par.xx;
% layer_index = find(x - (x(end) - thickness*1e-7) > 0, 1, 'first'); %finds index for the thickness of charged layer
% actual_thickness = (x(end)-x(layer_index))*1e7;%compensates for the numerical grid that might not take exactly thickness as thickness
% compensation_charges = -1 * actual_thickness/(x(layer_points_p_layer)*1e7) * charges;
% 
% par.dev_sub.Extra_charge(layer_index:end) = charges;
% par.dev_sub.Extra_charge(1:layer_points_p_layer) = compensation_charges;

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
dfplot.npx(solEq.el,tarr,xrange) % FIGURE

% Plot charge density, electric field, electrostatic potential as a function of distance
dfplot.rhoxFxVx(solEq.el,tarr,xrange) % FIGURE

% Plot energy level diagram as a function of distance
dfplot.ELx(solEq.el,tarr,xrange) % FIGURE

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

%% Store this run in an array 'runs' in the workspace
% Build a struct for this simulation

thisRun = struct();

% Basic JV info
thisRun.Vapp       = Vapp;                      % voltage vector
thisRun.Jtot       = J_CV.tot(:, ppos);         % total J at active-layer point

% Full solutions & parameters (handy if you want profiles later)
thisRun.solEq      = solEq;                     % equilibrium solution
thisRun.solCV      = solCV;                     % CV / JV sweep solution
thisRun.par        = par;                       % full parameter object

% Key device metrics
thisRun.Voc        = Voc;                       % open-circuit voltage
thisRun.PCE        = efficiency;                % power conversion efficiency [%]
thisRun.MPP_power  = MPP;                       % absolute minimum of P(V)
thisRun.Pin        = Pin;                       % incident power density [W cm^-2]

% Example: store whatever parameter you are varying between runs
% (edit these to match what you are actually sweeping)
thisRun.d_13 = par.d(1,3);                  % junction / layer thickness [cm]


% You can also add your own label if you like
thisRun.label = sprintf('Run at d(1,3) = %.3g cm', thisRun.d_13);

% Append to (or create) the 'runs' array in the base workspace
if evalin('base','exist(''runs'',''var'')')
    runs = evalin('base','runs');               % pull current runs
    runs(end+1) = thisRun;                      % append
else
    runs = thisRun;                             % first run
end

assignin('base','runs',runs);                   % push back to base workspace

% (Optional) display which run number this is
disp(['Stored run #', num2str(numel(runs)), ' in variable ''runs''']);