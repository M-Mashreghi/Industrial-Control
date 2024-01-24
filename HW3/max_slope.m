function max_slope = find_max_slope(points, times)
    % Check if there are at least two points
    if length(points) < 2 || length(times) < 2 || length(points) ~= length(times)
        error('Invalid input: at least two points and times are required.');
    end

    max_slope = -Inf;  % Initialize with negative infinity

    for i = 2:length(points)
        % Calculate slope for each pair of consecutive points and times
        slope = (points(i) - points(i - 1)) / (times(i) - times(i - 1));

        % Update max_slope if the current slope is greater
        if slope > max_slope
            max_slope = slope;
        end
    end
end