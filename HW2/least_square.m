function [parameters, regres_matrix] = least_square(y, u_t, t)
    y_vector = -1 * y;  % Negate y
    u_vector = u_t';    % Transpose u_t

    regres_matrix = [];

    for i = 1:length(t)
        if i < 6
            % Construct regression vectors for the first 5 points
            regres_vect_y = flip([zeros(1, 5 - i), y_vector(1:i - 1)']);
            regres_vect_u = flip([zeros(1, 5 - i), u_vector(1:i)']);
        else
            % Construct regression vectors for subsequent points
            regres_vect_y = flip(y_vector((i - 4):(i - 1))');
            regres_vect_u = flip(u_vector(i - 4:i)');
        end

        % Combine regression vectors
        regres_vect = [regres_vect_y, regres_vect_u];

        % Append to the regression matrix
        regres_matrix = [regres_matrix; regres_vect];
    end

    % Display the regression matrix
    disp('Regression Matrix:');
    disp(regres_matrix);

    % Calculate the parameters using the least squares method
    parameters = inv(regres_matrix' * regres_matrix) * regres_matrix' * y;
end
