function [regres_vects, parameters] = r_least_square(y, u_t, t)
    % Initialize variables
    y_vector = -1 * y;
    u_vector = u_t';
    p0 = 1e10 * eye(9);
    theta0 = zeros(9, 1);
    pm = p0;
    theta_m = theta0;
    
    % Initialize result arrays
    parameters = [];
    regres_vects = [];
    
    % Main loop
    for i = 1:length(t)
        % Create regression vector based on the current time index
        if i < 6
            regres_vect_y = flip([zeros(1, 5 - i), y_vector(1:i-1)']);
            regres_vect_u = flip([zeros(1, 5 - i), u_vector(1:i)']);
        else
            regres_vect_y = flip(y_vector((i-4):(i-1))');
            regres_vect_u = flip(u_vector(i-4:i)');
        end
        
        % Combine regression vectors
        regres_vect = [regres_vect_y, regres_vect_u];
        
        % Update regression vectors and parameters
        regres_vects = [regres_vects; regres_vect];
        p_m1 = pm - (pm * regres_vect' * regres_vect * pm) / (1 + regres_vect * pm * regres_vect');
        k_m1 = p_m1 * regres_vect';
        e_m1 = -1 * y_vector(i) - regres_vect * theta_m;
        theta_m1 = theta_m + k_m1 * e_m1;
        
        % Update variables for the next iteration
        theta_m = theta_m1;
        pm = p_m1;
        
        % Store parameters
        parameters = [parameters, theta_m];
    end
end
