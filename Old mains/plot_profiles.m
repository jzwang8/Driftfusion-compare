

load("Driftfusion/sweep_results/SF1_eta2.0_jd400_sp1e05_ND1e17_PhiOffn0p0539.mat")

this_solEq_el = thisRun.solEq.el;
this_solCV = thisRun.solCV;

tarr   = this_solEq_el.t(end);
xrange = [1.795*1e5, 1.803*1e5];

[u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(this_solEq_el);
% npx
figure;
dfplot.x2d(this_solEq_el, x, {n, p}, {'n', 'p'}, {'-','-'},'Carrier density [cm-3]', tarr, xrange, 0, 1)

% rhoxFxVx
rho = dfana.calcrho(this_solEq_el, "whole");
F = dfana.calcF(this_solEq_el, "whole");

figure;
subplot(3, 1, 1)
dfplot.x2d(this_solEq_el, x, {rho},{'\rho'},{'-'}, 'Charge density [cm-3]', tarr, xrange, 0, 0);

subplot(3, 1, 2)
dfplot.x2d(this_solEq_el, x, {F},{'F'},{'-'},'Electric field [Vcm-1]', tarr, xrange, 0, 0);

subplot(3, 1, 3)
dfplot.x2d(this_solEq_el, x, {V},{'V'},{'-'},'Electrostatic potential [V]', tarr, xrange, 0, 0);

% ELx
[Ecb, Evb, Efn, Efp] = dfana.calcEnergies(this_solEq_el);
figure;
dfplot.x2d(this_solEq_el, x, {Efn, Efp, Ecb, Evb}, {'E_{fn}', 'E_{fp}', 'E_{CB}', 'E_{VB}'},...
    {'--', '--', '-', '-'}, 'Energy [eV]', tarr, xrange, 0, 0)

% rx
x_sub = par.x_sub;
r = dfana.calcr(this_solEq_el, "sub");
figure;
dfplot.x2d(this_solEq_el, x_sub, {r.btb, r.srh, r.vsr},{'r_{btb}', 'r_{SRH}', 'r_{vsr}'},...
                {'-','-','-'}, 'Recombination rate (cm^{-3}s^{-1})', tarr, xrange, 0, 0);