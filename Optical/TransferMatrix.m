function [f_Si_interp, f_Tc_interp] = TransferMatrix(wavelengths)
% This funciton returns the fraction of absorbed light for silicon and
% tetracene (f_Si, f_Tc resp.) for the calculation of the generation rates.
% The data is calculated using the transfer matrix method and dependent on
% the wavelength of the illuminated light.
% The program uses linear interpolation/extrapolation to determine the 
% index of refraction for wavelengths not listed in the library.

%Data in IndRefr, Column names in IndRefr_names
load(['../Transfer matrix modeling/Results/20240907_TMatrixCalculation_30nmTc.mat'])
file_wavelengths = lambda;
file_wavelengths_begin = file_wavelengths(1);
file_wavelengths_end = file_wavelengths(end);

f_Si = Absorption(5,:)./sum(Absorption(:,:));
f_Tc = Absorption(2,:)./sum(Absorption(:,:));

% Interpolate/Extrapolate data linearly to desired wavelengths
f_Si_interp=interp1(file_wavelengths, f_Si, wavelengths, 'linear', 'extrap');
f_Tc_interp=interp1(file_wavelengths, f_Tc, wavelengths, 'linear', 'extrap');

% Set fraction of light absorbed in silicon outside of the data range to
% 1 and for Tc to 0
max_wavelength_index = find(wavelengths > file_wavelengths_end,1 , "first" );
min_wavelength_index = find(wavelengths < file_wavelengths_begin,1 , "last");
f_Si_interp(1:min_wavelength_index) = 1;
f_Si_interp(max_wavelength_index:end) = 1;
f_Tc_interp(1:min_wavelength_index) = 0;
f_Tc_interp(max_wavelength_index:end) = 0;
end