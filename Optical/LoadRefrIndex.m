% % Function LoadRefrIndex
% % This function returns the complex index of refraction spectra, ntotal, for the
% % material called 'name' for each wavelength value in the wavelength vector
% % 'wavelengths'.  The material must be present in the index of refraction
% % library 'Index_of_Refraction_library.xls'.  The program uses linear
% % interpolation/extrapolation to determine the index of refraction for
% % wavelengths not listed in the library.
% function [n_interp, k_interp] = LoadRefrIndex(name,wavelengths)
% 
% %Data in IndRefr, Column names in IndRefr_names
% [IndRefr,IndRefr_names]=readtable('Index_of_Refraction_library.xls');
% disp(readtable('Index_of_Refraction_library.xls');)
% 
% % Load index of refraction data in spread sheet, will crash if misspelled
% file_wavelengths=IndRefr(:,strmatch('Wavelength',IndRefr_names));
% n=IndRefr(:,strmatch(strcat(name,'_n'),IndRefr_names));
% k=IndRefr(:,strmatch(strcat(name,'_k'),IndRefr_names));
% 
% % Interpolate/Extrapolate data linearly to desired wavelengths
% n_interp=interp1(file_wavelengths, n, wavelengths, 'linear', 'extrap');
% k_interp=interp1(file_wavelengths, k, wavelengths, 'linear', 'extrap');
% 
% %Return interpolated complex index of refraction data
% %ntotal = n_interp+1i*k_interp; 
% 
% end

% Function LoadRefrIndex
% This function returns the complex index of refraction spectra, ntotal, for the
% material called 'name' for each wavelength value in the wavelength vector
% 'wavelengths'.  The material must be present in the index of refraction
% library 'Index_of_Refraction_library.xls'.  The program uses linear
% interpolation/extrapolation to determine the index of refraction for
% wavelengths not listed in the library.
function [n_interp, k_interp] = LoadRefrIndex(name,wavelengths)

% Load data from the spreadsheet using readtable
IndRefr = readtable('Index_of_Refraction_library.xls', 'VariableNamingRule', 'preserve');
% disp(IndRefr);

% Extract the wavelength column
file_wavelengths = IndRefr.Wavelength;  % Assuming 'Wavelength' is the column name

% Extract the refractive index (n) and extinction coefficient (k) columns
n_column = strcat(name, '_n');  % Construct column name for n
k_column = strcat(name, '_k');  % Construct column name for k

% Check if the columns exist in the table
if ismember(n_column, IndRefr.Properties.VariableNames) && ismember(k_column, IndRefr.Properties.VariableNames)
    n = IndRefr.(n_column);  % Access the column data for n
    k = IndRefr.(k_column);  % Access the column data for k
else
    error('Material name not found in the index of refraction library.');
end

% Interpolate/Extrapolate data linearly to desired wavelengths
n_interp = interp1(file_wavelengths, n, wavelengths, 'linear', 'extrap');
k_interp = interp1(file_wavelengths, k, wavelengths, 'linear', 'extrap');

%Return interpolated complex index of refraction data
%ntotal = n_interp + 1i * k_interp; 

end
