function alpha_interp = LoadAbsorbance(wavelengths)
% This function returns the absorption value at the wavelength laserlambda
% for the absorbance data specified in the xlsread(file name). Where the
% column of the absorbance data is called 'Absorbance'

%Data in IndRefr, Column names in IndRefr_names
[IndRefr,IndRefr_names]=xlsread('Tet_30.xlsx');

% Load index of refraction data in spread sheet, will crash if misspelled
file_wavelengths=IndRefr(:,strmatch('Wavelength',IndRefr_names));
file_wavelengths_begin = file_wavelengths(1);
file_wavelengths_end = file_wavelengths(end);
k = IndRefr(:,endsWith(IndRefr_names,'_k'));
alpha = 4*pi*k./(file_wavelengths *1e-7);

% Interpolate/Extrapolate data linearly to desired wavelengths
alpha_interp=interp1(file_wavelengths, alpha, wavelengths, 'linear', 'extrap');
% set absorbance data to zero if data is not available
max_wavelength_index = find(wavelengths > file_wavelengths_end,1 , "first" );
min_wavelength_index = find(wavelengths < file_wavelengths_begin,1 , "last");
alpha_interp(1:min_wavelength_index) = 0;
alpha_interp(max_wavelength_index:end) = 0;
end