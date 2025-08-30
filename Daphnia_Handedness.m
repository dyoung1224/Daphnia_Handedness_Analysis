clear all;

% ============================================================
%   Daphnia CW/CCW Rotational Bias Analysis
%
%   This script analyzes the trajectory of a *single* Daphnia 
%   to classify its turns into clockwise (CW) and counterclockwise (CCW). 
%   It also computes a handedness index describing directional bias.
%
%   HOW IT FITS IN THE PIPELINE:
%     Step 1: Convert TREX .npz tracking to CSV
%     Step 2: Calculate cumulative angular dynamics across all Daphnia
%     Step 3: (THIS SCRIPT) CW/CCW classification for an individual
% ============================================================


% ============ USER-DEFINED PARAMETERS =======================

% Dataset base name (e.g., experiment label)
base = 'Multiple_6-23-24-25';

% Directory containing CSV files for each tracked Daphnia
% NOTE: change 'path' and 'folder' to your actual location
inputDir = fullfile('path', 'to', 'your', 'csv', 'folder', base, '_csv');

% Number of tracked Daphnia in dataset
N = 50;

% Histogram bin size (cm/s) for estimating speed distribution
binSize = 0.5;

% Index of Daphnia to analyze (1-based index: d=3 means "Daphnia #3")
d = 3;

% Time window length (in seconds) for CW/CCW classification
timeWindow = 2;

% ============================================================



%% ========= Load Data for Selected Daphnia ==================

% Construct file name based on base name and chosen daphnia index
csvFile = fullfile(inputDir, sprintf('%s_daphnia%d.csv', base, d-1));

% Load trajectory data
data    = readtable(csvFile);

% Extract X, Y positions, time stamps, and frame rate
X    = data.X;
Y    = data.Y;
Time = data.Time;
FPS  = data.fps;

% Compute average frames per second (some datasets repeat fps values)
fps     = mean(FPS);
maxTau  = numel(Time) - 1;  % number of usable time intervals



%% ========= Calculate Speeds ================================

speedAll = [];  % container for all instantaneous speeds

% Loop through trajectory points and calculate speed = distance/time
for j = 1:maxTau
    dx = X(j+1) - X(j);
    dy = Y(j+1) - Y(j);
    dt = Time(j+1) - Time(j);

    % Only include valid values (exclude NaN, Inf, dt=0)
    if isfinite(dx) && isfinite(dy) && isfinite(dt) && dt ~= 0
        spd = sqrt(dx^2 + dy^2) / dt;  % speed in cm/s
        speedAll = [speedAll; spd];
    end
end

% Histogram the speed distribution
minSpeed = min(speedAll);
maxSpeed = max(speedAll);
sp_numBins = ceil((maxSpeed - minSpeed) / binSize);

[C, E] = histcounts(speedAll, sp_numBins);
PROB = C / sum(C);  % normalize histogram to probability

% Bin centers (for plotting/log transform)
binCenters = (E(1:end-1) + E(2:end)) / 2;

% Use only bins with nonzero probability for fitting
validIdx = PROB > 0;
logProbability = log(PROB(validIdx));
binCentersValid = binCenters(validIdx);

% Linear fit (log probability vs. speed) → slope = -1/v*
fitCoeffs = polyfit(binCentersValid, logProbability, 1);
slope = fitCoeffs(1);



%% ========= Apply Velocity Threshold v* =====================

% Compute velocity components
Vx = diff(X) ./ diff(Time);
Vy = diff(Y) ./ diff(Time);
V_total = sqrt(Vx.^2 + Vy.^2);

% Physics-based velocity threshold
v_star = 1 / abs(slope);

% Keep only velocities above threshold
valid_indices = (V_total >= v_star) & isfinite(V_total);

Vx_f = Vx(valid_indices);
Vy_f = Vy(valid_indices);

% Filtered time array (exclude last point b/c of diff())
Time_f = Time(1:end-1);
Time_f = Time_f(valid_indices);



%% ========= Angular Change per Step =========================

% Calculate incremental change in angle (delta_phi)
delta_phi = zeros(1, numel(Vx_f)-1);

for a = 1:numel(delta_phi)
    dVx = Vx_f(a+1) - Vx_f(a);
    dVy = Vy_f(a+1) - Vy_f(a);

    % Formula for angular velocity component
    numr = (Vx_f(a) * dVy) - (Vy_f(a) * dVx);
    denr = Vx_f(a)^2 + Vy_f(a)^2;

    delta_phi(a) = numr / denr;
end

% Normalize by 2π and accumulate angle over time
norm_delta_phi   = delta_phi ./ (2*pi);
cumulative_angle = cumsum(norm_delta_phi);



%% ========= Windowed CW/CCW Analysis ========================

% Divide trajectory into fixed-length time windows
numWindows = floor((max(Time_f) - min(Time_f)) / timeWindow);
slopeSigns = zeros(numWindows, 1);

% Progress bar for window analysis
h = waitbar(0, 'Calculating CW/CCW turns...');

% Track turn counts
cwTurns = 0;
ccwTurns = 0;

figure;
scatter(Time_f(1:end-1), cumulative_angle);
hold on;

% Loop over each window
for w = 1:numWindows
    waitbar(w / numWindows, h, sprintf('Processing window %d/%d...', w, numWindows));

    % Window boundaries
    startTime = (w - 1) * timeWindow + min(Time_f);
    endTime   = w * timeWindow + min(Time_f);

    windowIdx = (Time_f >= startTime) & (Time_f <= endTime);
    if sum(windowIdx) < 2
        continue; % skip windows with < 2 points
    end

    % Plot vertical line for each window boundary
    xline(Time_f(find(windowIdx,1,'last')));

    % Extract data for window
    timeWindowData = Time_f(windowIdx);
    angleWindowData = cumulative_angle(windowIdx);

    % Fit line: slope > 0 = CCW, slope < 0 = CW
    fitCoeffs_window = polyfit(timeWindowData, angleWindowData, 1);

    if fitCoeffs_window(1) > 0
        slopeSigns(w) = 1;   % CCW
    elseif fitCoeffs_window(1) < 0
        slopeSigns(w) = -1;  % CW
    end
end



%% ========= Detect CW ↔ CCW Switches ========================

% Compare consecutive slopes to detect changes
prevNonZeroIndex = 0;
for w = 2:numWindows
    waitbar((numWindows+w)/(2*numWindows), h, 'Analyzing slopes...');
    if slopeSigns(w) == 0, continue; end

    if prevNonZeroIndex == 0
        prevNonZeroIndex = w;
    else
        diffSlope = slopeSigns(w) - slopeSigns(prevNonZeroIndex);

        % +2 → CW turn, -2 → CCW turn
        if diffSlope == 2
            cwTurns = cwTurns + 1;
        elseif diffSlope == -2
            ccwTurns = ccwTurns + 1;
        end
        prevNonZeroIndex = w;
    end
end

close(h);



%% ========= Handedness Calculation ==========================

% Handedness = (CW - CCW) / (CW + CCW)

% Identify time indices of slope sign changes
index_star = [];
for w = 2:numWindows
    if slopeSigns(w) * slopeSigns(w-1) < 0
        % Midpoint of window = sign change time
        startTime = (w-1)*timeWindow + min(Time_f);
        endTime   = w*timeWindow + min(Time_f);
        windowIdx = (Time_f >= startTime) & (Time_f <= endTime);
        timeWindowData = Time_f(windowIdx);

        if ~isempty(timeWindowData)
            signChangeTime = mean(timeWindowData);
            [~, idx] = min(abs(Time_f - signChangeTime));
            index_star(end+1) = idx;
        end
    end
end

% Add last index
index_star = [index_star, numel(Time_f)];

% Measure angular change across each chunk
angleDifferences = zeros(length(index_star), 1);
for k = 1:length(index_star)
    if k == 1
        sIdx = 1; eIdx = index_star(k);
    else
        sIdx = index_star(k-1);
        eIdx = index_star(k);
    end
    eIdx = min(eIdx, numel(cumulative_angle));
    angleDifferences(k) = cumulative_angle(eIdx) - cumulative_angle(sIdx);
end

% Convert to cycles
angle_diff = angleDifferences ./ (2*pi);

% Sum positive vs negative cycles
CCW = sum(angle_diff(angle_diff > 0));
CW  = sum(angle_diff(angle_diff < 0));

fprintf('Total CW turns: %d\n', abs(CW));
fprintf('Total CCW turns: %d\n', CCW);

Handedness = (abs(CW)-CCW)/(abs(CW)+CCW);
fprintf('Handedness: %d\n', Handedness);



%% ========= Visualization with Color Coding =================

% New figure with interactive Pause/Resume buttons
figureHandle = figure;
setappdata(figureHandle, 'isPaused', false);

pauseButton  = uicontrol('Style','pushbutton','String','Pause', ...
    'Position',[20,20,60,30],'Callback',@(src,event)setappdata(figureHandle,'isPaused',true));
resumeButton = uicontrol('Style','pushbutton','String','Resume', ...
    'Position',[100,20,60,30],'Callback',@(src,event)setappdata(figureHandle,'isPaused',false));

% Plot base trajectory
plot(Time_f(1:end-1), cumulative_angle, 'k','LineWidth',1.5);
hold on;
xlim([min(Time_f), max(Time_f)]);
ylim([min(cumulative_angle)-1, max(cumulative_angle)+1]);

% Add live time label inside plot
time_label = text(min(Time_f)+0.05*range(Time_f), ...
    max(cumulative_angle), '', 'FontSize',12,'FontWeight','bold',...
    'BackgroundColor','w','EdgeColor','k','Margin',3);

% Animate window coloring
for w = 1:numWindows
    startTime = (w-1)*timeWindow + min(Time_f);
    endTime   = w*timeWindow + min(Time_f);

    % Color by slope sign
    if slopeSigns(w) == 1
        c = [1 0.8 0.8]; % Red = CCW
    elseif slopeSigns(w) == -1
        c = [0.8 0.8 1]; % Blue = CW
    else
        c = [1 1 1];     % White = no turn
    end

    % Shade window
    fill([startTime,endTime,endTime,startTime], ...
         [ylim(gca), fliplr(ylim(gca))], c, ...
         'EdgeColor','none','FaceAlpha',0.5);

    % Redraw trajectory
    plot(Time_f(1:end-1), cumulative_angle,'k','LineWidth',1.5);

    % Update live time label
    set(time_label,'String',sprintf('Time: %.2f s',endTime));
    pause(0.9);

    % Allow Pause/Resume
    if ~isvalid(figureHandle), break; end
    while getappdata(figureHandle,'isPaused')
        pause(0.9);
    end
end

% Legend (using invisible patches)
p1 = patch(NaN,NaN,[1 0.8 0.8]);
p2 = patch(NaN,NaN,[0.8 0.8 1]);
legend([p1,p2],'CCW (Red)','CW (Blue)','Location','NorthOutside','Orientation','horizontal');

xlabel('Time (s)');
ylabel('Cumulative Angle');
title('Cumulative Angle with CW and CCW Regions');
grid on;
hold off;
