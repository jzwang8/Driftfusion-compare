% This function calculates the generation by a tetracene (Tc)
% layer as singlet fission material. It utilizes the absorption coefficient
% data stored in the library folder 
function Gentot = tetracene(par)
    %% Constants
    h = 6.626e-34;     % Js
    c = 2.998e8;       % m/s
    thickness = 1;    % nm
    stepsize = 0.5;    % nm
    lambda = 280:stepsize:1200; % wavelengths range

    %% Import the absorption data (absorption will be set to 0 if the
    % wavelength is outside of the absorbance data file).
    laserlambda = par.laser_lambda1;
    alpha = LoadAbsorbance(lambda);
    source_type = par.light_source1;
    % alpha = LoadAbsorbance(laserlambda)/(thickness * log10(exp(1))); Only
    % use if the data is of absorbance
    
    switch source_type
        case 'AM15'
            I0 = lightsource('AM15', lambda);
        case 'laser'
            % Throw up error if LASERLAMBDA is not within the acceptable range
            if laserlambda < lambda(1) || laserlambda > lambda(end)
                msg = ['Error in tetracene - laser wavelength must be in the range ', num2str(lambda(1)), ' - ', num2str(lambda(end)), ' nm'];
                error(msg);
            end
            I0 = zeros(1,length(lambda));
            I0_index = lambda == laserlambda; 
            I0(I0_index) = 1e-3 * par.pulsepow / stepsize;   % convert to W cm-2 nm-1 (the pulse has a linewidth and is not an ideal delta function such that we have to devide by the stepsize)
        otherwise
            warning('Light source unknown. Usinng AM1.5 as default')
            I0 = lightsource('AM15', lambda);                % W cm-2 nm-1
    end
    % find x where rectangular pulse begins:
    x = par.xx;
    box_index = find(x - (x(end) - thickness*1e-7) > 0, 1, 'first');
    actual_thickness = (x(end)-x(box_index))*1e7;
    [fSi, fTc] = TransferMatrix(lambda);
    Eph = h * c ./ (lambda * 1e-9);
    % Gen = par.Tetracene_TF * fTc .*  alpha .* I0 ./ Eph ;
    Gen = par.Tetracene_TF .* fTc .* I0 ./ (Eph * actual_thickness * 1e-7);
    % Gen = par.Tetracene_TF *  alpha .* I0 ./ Eph;
    Gentot = trapz(lambda, Gen) * rectangularPulse(x(box_index), x(end), x); % is 0.5 at the edges!!!
end