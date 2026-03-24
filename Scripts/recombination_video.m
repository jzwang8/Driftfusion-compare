% MAKE SURE THAT DFPLOT.RX PLOTS IN FIGURE(17)!
% If you stop the run without finishing it, make sure to close(obj)
% otherwise it will not run again
% This script should make a video of the recombination as a function of
% position in the pn-junction
% Idea: plot df.plot for different times and stick them together to make a
% video

sol = solCV;
xmin = 2.80025e5;
xmax = 2.80033e5;
n_tot = length(sol.t);

obj = VideoWriter('20240830_recombinaion movie var32 neutral no Tc', 'MPEG-4');
obj.FrameRate = 15;
open(obj)
for n = 1 : 4 : n_tot
    dfplot.rx(sol,sol.t(n),[xmin, xmax])
    figure(17)
    txt = ['t = ',num2str(sol.t(n)), ' s'];
    Vapp = dfana.calcVapp(solCV);
    txt2 = ['V = ', num2str(Vapp(n)), ' V'];
    title(txt2)
    ylim([0,1e15])
    writeVideo(obj, getframe(figure(17)))
    % pause (0.1)
end
close(obj)