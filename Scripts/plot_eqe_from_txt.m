hold on
data = importdata('EQE_30um_woglass.xlsx') ;
wavelength = data.data(:,1);
eqe = data.data(:,2);
plot(wavelength,eqe)