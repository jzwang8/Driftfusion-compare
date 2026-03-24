% Test of the #photons going in are conserved
file_name = 'Input_files/pn_junction_var10';
initialise_df
par = pc(file_name);
wavelengths = 400:10:600;

G = zeros(length(wavelengths),1);
GTc = zeros(length(wavelengths),1);
GSi = zeros(length(wavelengths),1);
T = zeros(length(wavelengths),1);
IN = zeros(length(wavelengths),1);
for index = 1 : length(wavelengths)

par.light_source1 = 'laser';
par.laser_lambda1 = wavelengths(index);
par.pulsepow = 1;
par.Tetracene_TF = true;

% Calculate generation profile (as a function of position)
gSi = beerlambert(par, par.xx, par.light_source1, par.laser_lambda1, 0);
gTc = tetracene(par);
gSi = interp1(par.xx, gSi, par.x_sub);
gTc = interp1(par.xx, gTc, par.x_sub);

% Integrate over the whole device length to get the total number of
% absorbed photons/unit area
GTc(index) = trapz(par.x_sub, gTc);
GSi(index) = trapz(par.x_sub, gSi);
G(index) = GTc(index) + GSi(index);

% Calculate how much the #absorbed photons/area differs from the #incomming
% photons
h = 6.626e-34;     % Js
c = 2.998e8;       % m/s

I = beerlambertI(par, par.xx, par.light_source1, par.laser_lambda1, 0);
Plost = I(1,par.laser_lambda1-300); % calculates the residue power
E_photon = h*c/(par.laser_lambda1*1e-9);
T(index) = Plost/E_photon;
IN(index) = par.pulsepow*1e-3/E_photon;
end
Difference = (G+T)./(IN); % [%]
figure
plot(wavelengths, Difference)
xlabel('Wavelength (nm)')
ylabel('(P_{abs}+P_{trans})/P_{in} (-)')