% function [L, A] = exponentialFromGaussian(N0, x0, sigma)
%     % Define depth range
%     x = linspace(0, 300, 1000); % nm
% 
%     % Gaussian doping concentration function
%     Nd_gauss = N0 * exp(-((x - x0).^2) / (2 * sigma^2));
% 
%     % Find the minimum Gaussian doping at x = 300 nm
%     N_min_gauss = Nd_gauss(end); % Value at x = 300 nm
% 
%     % If max doping concentration A results in EF>0
%         % Set A to the required maximum value
%         A = getDopingConcFromEf(0); % cm^-3
% 
%         % Solve for L using total doping conservation
%         syms L x_sym
%         total_doping_gauss = trapz(x, Nd_gauss);  % nm*cm^-3
% 
%         % Solve for L using the constraint that total doping remains the same
%         eqn_L = int(A * exp(-x_sym/L), x_sym, 0, 300) == total_doping_gauss;
% 
%         % Solve for L numerically
%         L_sol = double(vpasolve(eqn_L, L, [1, 500])); 
% 
%     % Compute exponential doping profile
%     Nd_exp = A * exp(-x / L_sol);
% 
%     % Check total doping conservation
%     total_doping_exp = trapz(x, Nd_exp);  % nm*cm^-3
% 
%     % Ensure that the exponential profile matches Gaussian minimum at x = 300 nm
%     min_exp = A * exp(-300 / L_sol);
%     fprintf('Gaussian min: %e, Exponential min: %e\n', N_min_gauss, min_exp);
%     fprintf('Gaussian total doping conc: %e, Exponential total doping conc: %e\n', total_doping_gauss, total_doping_exp);
% 
%     % Plot Gaussian and Exponential profiles
%     figure;
%     plot(x, Nd_gauss, 'b', 'LineWidth', 2); hold on;
%     plot(x, Nd_exp, 'r--', 'LineWidth', 2);
%     xlabel('Depth (nm)');
%     ylabel('Doping Concentration (cm^{-3})');
%     legend('Gaussian Profile', 'Exponential Profile');
%     title('Gaussian vs Exponential Doping Profiles');
%     grid on;
% 
%     % Return L and A
%     L = L_sol;
% end

% function [L, A] = exponentialFromGaussian(N0, x0, sigma)
%     % Define depth range
%     x = linspace(0, 300, 1000); % nm
% 
%     % Gaussian doping concentration function
%     Nd_gauss = N0 * exp(-((x - x0).^2) / (2 * sigma^2));
% 
%     % Find the minimum Gaussian doping at x = 300 nm
%     N_min_gauss = Nd_gauss(end); % Value at x = 300 nm
% 
%     % Determine A based on doping concentration constraint
%     Nd
%     if N_min_gauss > getDopin0
%         A = getDopingConcFromEf(0); % cm^-3
%     else
%         A = N0; % Keep the original peak doping if no adjustment is needed
%     end
% 
%     % Solve for L using total doping conservation
%     syms L x_sym
%     total_doping_gauss = trapz(x, Nd_gauss);  % nm*cm^-3
% 
%     % Solve for L using the constraint that total doping remains the same
%     eqn_L = int(A * exp(-x_sym/L), x_sym, 0, 300) == total_doping_gauss;
% 
%     % Solve for L numerically within reasonable bounds
%     L_sol = double(vpasolve(eqn_L, L, [1, 500])); 
% 
%     % Ensure L is valid, otherwise, assign a fallback value
%     if isempty(L_sol) || L_sol <= 0
%         error('Could not find a valid solution for L');
%     end
% 
%     % Compute exponential doping profile
%     Nd_exp = A * exp(-x / L_sol);
% 
%     % Check total doping conservation
%     total_doping_exp = trapz(x, Nd_exp);  % nm*cm^-3
% 
%     % Ensure that the exponential profile matches Gaussian minimum at x = 300 nm
%     min_exp = A * exp(-300 / L_sol);
% 
%     fprintf('Gaussian min: %e, Exponential min: %e\n', N_min_gauss, min_exp);
%     fprintf('Gaussian total doping conc: %e, Exponential total doping conc: %e\n', total_doping_gauss, total_doping_exp);
% 
%     % Plot Gaussian and Exponential profiles
%     figure;
%     plot(x, Nd_gauss, 'b', 'LineWidth', 2); hold on;
%     plot(x, Nd_exp, 'r--', 'LineWidth', 2);
%     xlabel('Depth (nm)');
%     ylabel('Doping Concentration (cm^{-3})');
%     legend('Gaussian Profile', 'Exponential Profile');
%     title('Gaussian vs Exponential Doping Profiles');
%     grid on;
% 
%     % Return L and A
%     L = L_sol;
% end

function [L, A] = exponentialFromGaussian(N0, x0, sigma)
    % Define depth range
    x = linspace(0, 300, 1000); % nm
    
    % Gaussian doping concentration function
    Nd_gauss = N0 * exp(-((x - x0).^2) / (2 * sigma^2));
    
    % Find the minimum Gaussian doping at x = 300 nm
    N_min_gauss = Nd_gauss(end); % Value at x = 300 nm
    
    % Compute total doping under the Gaussian profile
    total_doping_gauss = trapz(x, Nd_gauss);  % nm*cm^-3
    
    % Define function for numerical root-finding to solve for L
    total_doping_exp_func = @(L) (N_min_gauss / exp(-300/L)) * L * (1 - exp(-300/L)) - total_doping_gauss;
    
    % Solve for L numerically using fzero
    L_sol = fzero(total_doping_exp_func, 100);  % Initial guess of 100 nm
    
    % Ensure L is valid
    if isempty(L_sol) || isnan(L_sol) || L_sol <= 0
        error('Could not find a valid solution for L');
    end
    
    % Solve for A using the minimum doping condition
    A_sol = N_min_gauss / exp(-300 / L_sol);
    
    % Ensure A is valid
    if isempty(A_sol) || isnan(A_sol) || A_sol <= 0
        error('Could not find a valid solution for A');
    end

    if getEfFromDopingConc(A_sol) > 0
        % Set A to the required maximum value
        A_sol = getDopingConcFromEf(0); % cm^-3

        % Solve for L using total doping conservation
        syms L x_sym
        total_doping_gauss = trapz(x, Nd_gauss);  % nm*cm^-3

        % Solve for L using the constraint that total doping remains the same
        eqn_L = int(A_sol * exp(-x_sym/L), x_sym, 0, 300) == total_doping_gauss;

        % Solve for L numerically
        L_sol = double(vpasolve(eqn_L, L, [1, 500])); 
    end
    
    % Compute exponential doping profile
    Nd_exp = A_sol * exp(-x / L_sol);
    
    % Compute total doping under the exponential profile for verification
    total_doping_exp = trapz(x, Nd_exp);  % nm*cm^-3
    
    % Ensure that the exponential profile matches Gaussian minimum at x = 300 nm
    min_exp = A_sol * exp(-300 / L_sol);
    
    fprintf('Gaussian min: %e, Exponential min: %e\n', N_min_gauss, min_exp);
    fprintf('Gaussian total doping conc: %e, Exponential total doping conc: %e\n', total_doping_gauss, total_doping_exp);
    
    % Plot Gaussian and Exponential profiles
    figure;
    plot(x, Nd_gauss, 'b', 'LineWidth', 2); hold on;
    plot(x, Nd_exp, 'r--', 'LineWidth', 2);
    xlabel('Depth (nm)');
    ylabel('Doping Concentration (cm^{-3})');
    legend('Gaussian Profile', 'Exponential Profile');
    title('Gaussian vs Exponential Doping Profiles');
    grid on;
    
    % Return L and A
    L = L_sol;
    A = A_sol;
end
