% This script tested the hypothesis of the equivalence of calculating the
% generation rate of the solar spectrum (choosing AM15 in the generation 
% function) to a sum of generation rates calculated by the laser light
% source option and modulated by the power of the AM15 spectrum.

clear all
initialise_df;
par = pc('input_files/pn_junction_var10');
par.pulsepow = 1;        % [mW cm-2 nm-1]
par.g_inj = 0;           % [cm-3s-1]
laserlambda = 301:1:600;

%% calculate the generation of the AM1.5G spectrum:
gxsun = generation(par, 'AM15', 0); % calculates generation profile of the solar spectrum

%% Calculate the generation of a laser tuned by the AM1.5G spectrum:
g_laser = zeros(length(laserlambda),length(par.x_sub));
Isun = lightsource('AM15', laserlambda);
Ilaser = par.pulsepow * 1e-3;
for i=1:length(laserlambda)
    g_laser(i,:) =  Isun(i)/Ilaser * generation(par, 'laser', laserlambda(i));
end
gxsun2 = trapz(laserlambda,g_laser,1);
%% Compare the two generation solutions to see if they are the same
figure
hold on
plot(par.x_sub, 100*(gxsun2-gxsun)./gxsun, 'DisplayName', 'AM1.5G generation rate')
% plot(par.x_sub, gxsun2, 'DisplayName','Generation rate from lasers')
xlabel('Position, x (cm)')
ylabel('g_{laser}/g_{1.5G}-1 (%)')
set(legend,'FontSize',16);
set(legend,'EdgeColor',[1 1 1]);
legend show
hold off