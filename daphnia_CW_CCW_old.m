clear all;
% Set the directory containing the CSV files
% pre = 'D:\\David\\Output\\';
% base = 'single_7_23';
% csv = '_csv';

pre = 'C:\\Users\\swany\\Videos\\Daphnia\\6-23,24-25\\';
base = 'Multiple_6-23-24-25';
csv = '_csv';

inputDir = strcat(pre,base,csv);

% Number of daphnia
N = 50;

% Ask user for bin size
binSize = 0.5;

% Specific daphnia to focus on
d = 3;

timeWindow = 2; % seconds

%%  
% Initialize arrays to store the data for the specific daphnia
X = [];
Y = [];
Time = [];
FPS = [];

% Construct the file name for the specific daphnia
csvFile = fullfile(inputDir, sprintf('%s_daphnia%d.csv', base, d-1));

% Read the CSV file for the specific daphnia
data = readtable(csvFile);

% Store the data in the corresponding arrays
X = data.X;
Y = data.Y;
Time = data.Time;
FPS = data.fps;

% Determine the length of the Time field for the specific daphnia
minLength = numel(Time);
maxTau = minLength - 1;

% Calculate frames per second (FPS)
fps = mean(FPS);

% Initialize array for speeds
speed = zeros(maxTau, 1);
speedAll = []; % Array to store all speeds

% Calculate instantaneous speeds, skipping invalid values
for j = 1:maxTau
    dx = X(j+1) - X(j); % change in X
    dy = Y(j+1) - Y(j); % change in Y
    dt = Time(j+1) - Time(j); % change in Time
    
    % Skip any invalid calculations like inf or dt = 0
    if isfinite(dx) && isfinite(dy) && isfinite(dt) && dt ~= 0
        % Calculate speed (distance / time)
        speed(j) = sqrt(dx^2 + dy^2) / dt;
        
        % Store valid speeds
        speedAll = [speedAll; speed(j)];
    else
        % Skip this point if invalid
        continue;
    end
end

% Determine the number of bins based on user input for bin size
minSpeed = min(speedAll);
maxSpeed = max(speedAll);
sp_numBins = ceil((maxSpeed - minSpeed) / binSize);

% Calculate the probability distribution
[C, E] = histcounts(speedAll, sp_numBins);
PROB = C / sum(C); % Normalize counts to get probability

% Calculate the centers (median values) of the histogram bins
binCenters = (E(1:end-1) + E(2:end)) / 2;

% Apply a log transformation to the probability, but exclude any zero probabilities to avoid log(0) issues
validIdx = PROB > 0;
logProbability = log(PROB(validIdx));
binCentersValid = binCenters(validIdx);

% Perform a linear fit on the valid data points
fitCoeffs = polyfit(binCentersValid, logProbability, 1); % Linear fit (1st degree polynomial)

x0 = [-1, 0];
fit_function = fittype( @(a,b,x) a.* x + b);
[fitted_curve, gof] = fit(binCentersValid', logProbability', fit_function, 'StartPoint', x0);

% Extract the slope from the fit coefficients
slope = fitted_curve.a;

% Generate the linear fit values
fitLine = polyval(fitCoeffs, binCentersValid);

% % Plot the original data
% figure;
% bar(binCenters, log(PROB), 'histc');
% hold on;
% 
% % Plot the linear fit on top of the histogram
% plot(binCentersValid, fitLine, 'r-', 'LineWidth', 2);
% hold off;
% 
% xlabel('Speed (cm/s)');
% ylabel('Log Probability');
% title(sprintf('Linear Fit to Histogram for Daphnia %d | Slope: %.4f', d, slope));

%%

% Calculate the instantaneous velocity components for the specific daphnia
Vx = diff(X) ./ diff(Time); % Velocity in the X direction
Vy = diff(Y) ./ diff(Time); % Velocity in the Y direction

% Calculate the total velocity magnitude at each time step
V_total = sqrt(Vx.^2 + Vy.^2);

% Set v_star as the threshold (inverse of the slope)
v_star = 1 / abs(slope);  % v_star is the inverse of the slope

% Find the indices where the velocity is greater than or equal to v_star
valid_indices = (V_total >= v_star) & isfinite(V_total);

% Retain only the valid velocities, corresponding components, and times
Vx_filtered = Vx(valid_indices);
Vy_filtered = Vy(valid_indices);
V_total_filtered = V_total(valid_indices);

% The time for each velocity will be Time(1:end-1) because of the diff() calculation
Time_filtered = Time(1:end-1);  % Use time excluding the last point (since diff reduces the size by 1)
Time_filtered = Time_filtered(valid_indices);  % Filter the time based on valid velocities

%%%%%% Angle Calculation and CW/CCW Counting %%%%%%%
delta_phi = zeros(1, numel(Vx_filtered)-1);

for a = 1:(numel(delta_phi))
    delta_Vx = Vx_filtered(a+1) - Vx_filtered(a);
    delta_Vy = Vy_filtered(a+1) - Vy_filtered(a);
    numerator = (Vx_filtered(a) * delta_Vy) - (Vy_filtered(a) * delta_Vx);
    denominator = (Vx_filtered(a)^2 + Vy_filtered(a)^2);
    delta_phi(a) = numerator / denominator;
end

% Normalize the angle increments
norm_delta_phi = delta_phi ./ (2 * pi);

% Calculate the cumulative sum of the normalized angle increments
cumulative_angle = cumsum(norm_delta_phi);

%%
%%TRAJECTORY--COMMENT OUT

% Find the number of time windows
numWindows = floor((max(Time_filtered) - min(Time_filtered)) / timeWindow);

% Initialize array for slopes (+1 for CCW, -1 for CW, 0 for flat slope)
slopeSigns = zeros(numWindows, 1);

% Initialize waitbar
h = waitbar(0, 'Calculating CW/CCW turns...');

% Initialize counters for CW and CCW turns
cwTurns = 0;  % Clockwise (CW) turns
ccwTurns = 0; % Counterclockwise (CCW) turns

figure;
scatter(Time_filtered(1:end-1), cumulative_angle);
hold on;

% Loop over each time window
for w = 1:numWindows
    % Update the progress bar
    waitbar(w / numWindows, h, sprintf('Processing window %d of %d...', w, numWindows));

    % Get the start and end indices for this time window
    startTime = (w - 1) * timeWindow + min(Time_filtered);
    endTime = w * timeWindow + min(Time_filtered);

    % Extract the data in this time window
    windowIdx = (Time_filtered >= startTime) & (Time_filtered <= endTime);
    if sum(windowIdx) < 2  % Skip windows with less than two data points
        continue;
    end

    ID = find(windowIdx == 1);
    xline(Time_filtered(ID(end)))

    % Extract the cumulative angle data in this window
    timeWindowData = Time_filtered(windowIdx);
    angleWindowData = cumulative_angle(windowIdx);

    % Perform a linear fit to the cumulative angle in the current window
    fitCoeffs_window = polyfit(timeWindowData, angleWindowData, 1);

    % Assign +1 for CCW (positive slope), -1 for CW (negative slope), and 0 for flat slope
    if fitCoeffs_window(1) > 0
        slopeSigns(w) = 1; % CCW
    elseif fitCoeffs_window(1) < 0
        slopeSigns(w) = -1; % CW
    else
        slopeSigns(w) = 0; % Flat slope (zero slope)
    end
end

% Now check the difference between consecutive non-zero values in slopeSigns
prevNonZeroIndex = 0; % Initialize to track the last non-zero slope index

for w = 2:numWindows
    % Update the progress bar for the second pass
    waitbar((numWindows + w) / (2 * numWindows), h, sprintf('Analyzing slopes %d of %d...', w, numWindows));

    % Skip if the current or previous value is zero
    if slopeSigns(w) == 0
        continue; % Skip zero slopes
    end

    if prevNonZeroIndex == 0
        % Find the first non-zero value
        prevNonZeroIndex = w;
    else
        % Calculate the difference between the current and the previous non-zero slope
        diffSlope = slopeSigns(w) - slopeSigns(prevNonZeroIndex);

        % Check if the difference is +2 (indicating CW turn) or -2 (indicating CCW turn)
        if diffSlope == 2
            cwTurns = cwTurns + 1;  % Tally for CW turn
        elseif diffSlope == -2
            ccwTurns = ccwTurns + 1; % Tally for CCW turn
        end

        % Update the previous non-zero index
        prevNonZeroIndex = w;
    end
end

% Close the progress bar
close(h);

%% Identify the sign change index in original time array

% Initialize an array to store the indices in Time_filtered where sign changes occur
index_star = [];

% Loop over slopeSigns to detect changes in sign
for w = 2:numWindows
    % Check if there is a change in sign (product of consecutive elements is negative)
    if slopeSigns(w) * slopeSigns(w - 1) < 0
        % Find the start and end time for the current window
        startTime = (w - 1) * timeWindow + min(Time_filtered);
        endTime = w * timeWindow + min(Time_filtered);

        % Extract the time window data from Time_filtered
        windowIdx = (Time_filtered >= startTime) & (Time_filtered <= endTime);
        timeWindowData = Time_filtered(windowIdx);

        % Find the mid-point of this time window to represent the approximate time of sign change
        if ~isempty(timeWindowData)
            signChangeTime = mean(timeWindowData); % You can choose other methods if needed
            % Find the closest index in Time_filtered corresponding to the sign change time
            [~, closestIdx] = min(abs(Time_filtered - signChangeTime));
            index_star(end + 1) = closestIdx; % Store the index
        end
    end
end

index_star = [index_star, numel(Time_filtered)];

% Initialize an array to store the angle differences for each chunk
angleDifferences = zeros(length(index_star), 1);

% Set the first chunk starting at the first value of cumulative_angle
firstValue = cumulative_angle(1);

% Loop through each chunk based on the index_star values
for k = 1:length(index_star)
    if k == 1
        % For the first chunk, calculate the difference from the first value of cumulative_angle
        chunkStartIdx = 1; % First value of cumulative_angle
        chunkEndIdx = index_star(k);
    else
        % For subsequent chunks, use consecutive index_star values
        chunkStartIdx = index_star(k-1);
        chunkEndIdx = index_star(k);
    end

    % Ensure that end index does not exceed the array bounds
    if chunkEndIdx > numel(cumulative_angle)
        chunkEndIdx = numel(cumulative_angle); % Adjust to last index if it exceeds
    end

    % Compute the difference in cumulative angles between the start and end of the chunk
    angleDifferences(k) = cumulative_angle(chunkEndIdx) - cumulative_angle(chunkStartIdx);
end

angle_diff = angleDifferences ./ (2*pi);

CCW = sum(angle_diff(angle_diff > 0));
CW = sum(angle_diff(angle_diff < 0));

% Output the count of CW and CCW turns
fprintf('Total CW turns: %d\n', abs(CW));
fprintf('Total CCW turns: %d\n', CCW);

% Handedness
Handedness = (abs(CW)-CCW)/(abs(CW)+CCW);
fprintf('Handedness: %d\n', Handedness)

%% Color the CW and CCW regions

% Create a figure window for the plot
figureHandle = figure;

% Add "Pause" and "Resume" buttons to the figure
pauseButton = uicontrol('Style', 'pushbutton', 'String', 'Pause', ...
    'Position', [20, 20, 60, 30], 'Callback', @(src, event) setappdata(figureHandle, 'isPaused', true));

resumeButton = uicontrol('Style', 'pushbutton', 'String', 'Resume', ...
    'Position', [100, 20, 60, 30], 'Callback', @(src, event) setappdata(figureHandle, 'isPaused', false));

% Set initial state for 'isPaused'
setappdata(figureHandle, 'isPaused', false);

% Plot initial cumulative angle
h = plot(Time_filtered(1:end-1), cumulative_angle, 'k', 'LineWidth', 1.5); % Initial plot of cumulative_angle
hold on;

% Get original axis limits (xlim and ylim) from the first plot to keep consistent scaling
original_xlim = [min(Time_filtered), max(Time_filtered)];
original_ylim = [min(cumulative_angle) - 1, max(cumulative_angle) + 1]; % Adjust this based on your original plot's y limits

% Set axis limits
xlim(original_xlim);
ylim(original_ylim);

% Place the time label within the graph frame
time_label = text(original_xlim(1) + 0.05 * range(original_xlim), original_ylim(2) - 0.1 * range(original_ylim), '', ...
    'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 3);

% Loop through each window chunk to animate
for w = 1:numWindows
    % Define the start and end of the time window
    startTime = (w - 1) * timeWindow + min(Time_filtered);
    endTime = w * timeWindow + min(Time_filtered);

    % Determine the color based on the slope sign
    if slopeSigns(w) == 1  % CCW (Counterclockwise)
        rectangleColor = [1 0.8 0.8]; % Light red
    elseif slopeSigns(w) == -1  % CW (Clockwise)
        rectangleColor = [0.8 0.8 1]; % Light blue
    else
        rectangleColor = [1 1 1]; % White (no turn, flat slope)
    end

    % Draw the rectangle, using the original y-limits to match the original plot scale
    fill([startTime, endTime, endTime, startTime], [original_ylim(1), original_ylim(1), original_ylim(2), original_ylim(2)], ...
        rectangleColor, 'EdgeColor', 'none', 'FaceAlpha', 0.5);

    % Update the cumulative angle plot on top of the rectangles
    plot(Time_filtered(1:end-1), cumulative_angle, 'k', 'LineWidth', 1.5);

    % Update the real-time label
    set(time_label, 'String', sprintf('Time: %.2f s', endTime));

    % Pause to create an animation effect
    pause(0.9); % Adjust pause duration to control the animation speed

    % Check if the figure is still open (user did not close it)
    if ~isvalid(figureHandle)
        break;
    end

    % Check the pause state and wait while paused
    while getappdata(figureHandle, 'isPaused')
        pause(0.9); % Small pause to avoid busy-waiting
    end
end

% Set the axis limits back to the original
xlim(original_xlim);
ylim(original_ylim);

% Create invisible patches for the legend
p1 = patch(NaN, NaN, [1 0.8 0.8]); % CCW (light red)
p2 = patch(NaN, NaN, [0.8 0.8 1]); % CW (light blue)

% Add legend for the rectangles (CW and CCW)
legend([p1, p2], 'CCW (Red)', 'CW (Blue)', 'Location', 'NorthOutside', 'Orientation', 'horizontal');

xlabel('Time (s)');
ylabel('Cumulative Angle');
title('Cumulative Angle with CW and CCW Regions');
grid on;

hold off;
