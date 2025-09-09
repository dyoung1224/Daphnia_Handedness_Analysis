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


%% ================== USER SETTINGS ==================
% Dataset base name (e.g., experiment label)
base = 'Multiple_6-23-24-25';

% Directory containing CSV files for each tracked Daphnia
% NOTE: change 'path' and 'folder' to your actual location
inputDir = fullfile('path', 'to', 'your', 'csv', 'folder', base, '_csv');

d = 3;                        % 1-based selection (file uses d-1)
csvFile = fullfile(inputDir, sprintf('%s_daphnia%d.csv', base, d-1));

binSize    = 0.5;             % speed histogram bin size (cm/s)
timeWindow = 2.0;             % window size for CW/CCW (s)

stepK     = 5;                % decimation for animation (bigger = faster, fewer frames)
basePause = 0.05;             % base frame delay (seconds); slider scales this
%% ===================================================

%% -------- Load data --------
data = readtable(csvFile);
X = data.X; Y = data.Y; Time = data.Time;
if ismember('fps', data.Properties.VariableNames)
    fpsVec = data.fps; %#ok<NASGU>
end

%% -------- Speed histogram -> slope -> v_star --------
dX = diff(X); dY = diff(Y); dT = diff(Time);
valid = isfinite(dX) & isfinite(dY) & isfinite(dT) & (dT ~= 0);
speedInst = sqrt(dX(valid).^2 + dY(valid).^2) ./ dT(valid);
if isempty(speedInst), error('No valid speed samples.'); end

minSpeed = min(speedInst); maxSpeed = max(speedInst);
sp_numBins = max(5, ceil((maxSpeed - minSpeed) / binSize));
[counts, edges] = histcounts(speedInst, sp_numBins);
prob = counts / sum(counts);
centers = 0.5*(edges(1:end-1) + edges(2:end));
idxPos = prob > 0; logP = log(prob(idxPos)); cOk = centers(idxPos);
if numel(cOk) < 2, error('Not enough positive-prob bins to fit.'); end
p = polyfit(cOk, logP, 1); slope = p(1);
if slope == 0 || ~isfinite(slope)
    warning('Slope zero/non-finite; using 75th percentile fallback.');
    v_star = prctile(speedInst, 75);
else
    v_star = 1/abs(slope);
end

%% -------- Velocities + v_star filter --------
Vx = dX ./ dT; Vy = dY ./ dT; V_total = sqrt(Vx.^2 + Vy.^2);
Time_vel = Time(1:end-1);
valid_indices = (V_total >= v_star) & isfinite(V_total) & isfinite(Time_vel);

Vx_f = Vx(valid_indices); Vy_f = Vy(valid_indices);
T_f  = Time_vel(valid_indices);
X_pos_all = X(1:end-1); Y_pos_all = Y(1:end-1);
X_pos = X_pos_all(valid_indices); Y_pos = Y_pos_all(valid_indices);
if numel(T_f) < 3, error('Too few samples after v_* filter.'); end

%% -------- Cumulative angle (turns) --------
dVx = diff(Vx_f); dVy = diff(Vy_f);
Vx_k = Vx_f(1:end-1); Vy_k = Vy_f(1:end-1);
den = Vx_k.^2 + Vy_k.^2; den(den==0) = eps;
delta_phi = (Vx_k .* dVy - Vy_k .* dVx) ./ den; % radians
norm_delta_phi = delta_phi ./ (2*pi);           % turns
cumulative_angle = cumsum(norm_delta_phi);
T_angle = T_f(1:end-1);

%% -------- Windowing & slope signs (robust) --------
t0 = min(T_angle); t1 = max(T_angle);
numWindows = max(1, floor((t1 - t0) / timeWindow));
slopeSigns = zeros(1, numWindows);

binIdx = ceil((T_angle - t0) / timeWindow);
binIdx(binIdx < 1) = 1; binIdx(binIdx > numWindows) = numWindows;

for w = 1:numWindows
    mask = (binIdx == w) & isfinite(T_angle) & isfinite(cumulative_angle);
    if nnz(mask) < 2, slopeSigns(w) = 0; continue; end
    xw = T_angle(mask); yw = cumulative_angle(mask);
    n = numel(xw); Sx = sum(xw); Sy = sum(yw);
    Sxx = sum(xw.^2); Sxy = sum(xw.*yw);
    denom = n*Sxx - Sx^2;
    if ~isfinite(denom) || denom == 0
        m = 0;
    else
        m = (n*Sxy - Sx*Sy) / denom; if ~isfinite(m), m = 0; end
    end
    slopeSigns(w) = sign(m); % +1 CCW, -1 CW, 0 flat
end

%% -------- Optional: CW/CCW totals & handedness (quick report) --------
index_star = [];
for w = 2:numWindows
    if slopeSigns(w) * slopeSigns(w-1) < 0
        startTime = (w - 1)*timeWindow + t0;
        endTime   =  w     *timeWindow + t0;
        signChangeTime = 0.5*(startTime + endTime);
        [~, closestIdx] = min(abs(T_angle - signChangeTime));
        index_star(end+1) = closestIdx; %#ok<AGROW>
    end
end
index_star = [index_star, numel(T_angle)];
angleDifferences = zeros(numel(index_star),1);
firstValue = cumulative_angle(1);
for k = 1:numel(index_star)
    if k == 1
        sIdx = 1; eIdx = index_star(k);
        angleDifferences(k) = cumulative_angle(eIdx) - firstValue;
    else
        sIdx = index_star(k-1); eIdx = min(index_star(k), numel(cumulative_angle));
        angleDifferences(k) = cumulative_angle(eIdx) - cumulative_angle(sIdx);
    end
end
CCW_turns = sum(angleDifferences(angleDifferences > 0));
CW_turns  = sum(angleDifferences(angleDifferences < 0));
Handedness = (abs(CW_turns) - CCW_turns) / (abs(CW_turns) + CCW_turns + eps);
fprintf('Total CW turns:  %.3f\n', abs(CW_turns));
fprintf('Total CCW turns: %.3f\n', CCW_turns);
fprintf('Handedness:      %.3f\n', Handedness);

%% -------- FAST side-by-side animation with SPEED + SCRUB --------
pad = 0.05;
xlim_traj = [min(X_pos)-pad*range(X_pos), max(X_pos)+pad*range(X_pos)];
ylim_traj = [min(Y_pos)-pad*range(Y_pos), max(Y_pos)+pad*range(Y_pos)];
xlim_angle = [t0, t1];
ylim_angle = [min(cumulative_angle)-1, max(cumulative_angle)+1];

fig = figure('Name','Daphnia CW/CCW & Trajectory (Fast)','Color','w','Renderer','opengl');
tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');

% Left: trajectory
axL = nexttile; hold(axL,'on'); axis(axL,'equal'); box(axL,'on'); grid(axL,'on');
xlabel(axL,'X (cm)'); ylabel(axL,'Y (cm)'); title(axL,'2D trajectory');
xlim(axL, xlim_traj); ylim(axL, ylim_traj);
plot(axL, X_pos, Y_pos, '-', 'Color',[0.8 0.8 0.8], 'LineWidth',1);
hPath = plot(axL, nan, nan, '-', 'LineWidth', 2);
hDot  = plot(axL, nan, nan, 'o', 'MarkerSize', 6, 'MarkerFaceColor','k', 'MarkerEdgeColor','k');
trajLabel = text(axL, xlim_traj(1)+0.05*range(xlim_traj), ...
                      ylim_traj(2)-0.08*range(ylim_traj), '', ...
                      'FontSize',11,'FontWeight','bold', ...
                      'BackgroundColor','w','EdgeColor','k','Margin',3);

% Right: bands + cumulative angle
axR = nexttile; hold(axR,'on'); box(axR,'on'); grid(axR,'on');
xlabel(axR,'Time (s)'); ylabel(axR,'Cumulative angle (turns)');
title(axR,'CW/CCW windows + cumulative angle');
xlim(axR, xlim_angle); ylim(axR, ylim_angle);

stripe = slopeSigns;               % -1,0,+1
stripe(stripe==-1) = 1;            % 1 -> blue (CW)
stripe(stripe==0)  = 2;            % 2 -> white (flat)
stripe(stripe==+1) = 3;            % 3 -> red (CCW)
imagesc(axR, [t0 t1], ylim_angle, stripe);
colormap(axR, [0.8 0.8 1; 1 1 1; 1 0.8 0.8]);
set(axR,'YDir','normal');
alpha(axR.Children(1), 0.4);

plot(axR, T_angle, cumulative_angle, 'k', 'LineWidth', 1.2);
hNow = line(axR, [t0 t0], ylim_angle, 'Color',[0.2 0.2 0.2], 'LineStyle',':');
angLabel = text(axR, xlim_angle(1)+0.05*range(xlim_angle), ...
                     ylim_angle(2)-0.08*range(ylim_angle), '', ...
                     'FontSize',11,'FontWeight','bold', ...
                     'BackgroundColor','w','EdgeColor','k','Margin',3);

% --- Controls ---
uicontrol('Style','pushbutton','String','Pause','Position',[20,20,60,28], ...
    'Callback', @(~,~) setappdata(fig,'isPaused',true));
uicontrol('Style','pushbutton','String','Resume','Position',[90,20,60,28], ...
    'Callback', @(~,~) setappdata(fig,'isPaused',false));

% Speed slider
uicontrol('Style','text','String','Speed', 'Position',[165,20,40,18], ...
    'HorizontalAlignment','left','BackgroundColor',get(fig,'Color'));
speedSlider = uicontrol('Style','slider','Min',0.25,'Max',4,'Value',1, ...
    'Position',[205,20,140,18], ...
    'Callback', @(s,~) setappdata(fig,'speedFactor', get(s,'Value')));
speedReadout = uicontrol('Style','text','String','1.00x', ...
    'Position',[350,20,50,18], 'HorizontalAlignment','left', ...
    'BackgroundColor',get(fig,'Color'));

% Scrub slider (seek by *time*, not index)
uicontrol('Style','text','String','Start at', 'Position',[415,20,50,18], ...
    'HorizontalAlignment','left','BackgroundColor',get(fig,'Color'));
scrubSlider = uicontrol('Style','slider','Min',t0,'Max',t1,'Value',t0, ...
    'Position',[470,20,240,18]);
scrubReadout = uicontrol('Style','text','String',sprintf('%.2fs',t0), ...
    'Position',[715,20,60,18],'HorizontalAlignment','left', ...
    'BackgroundColor',get(fig,'Color'));

% Play-from-here button
uicontrol('Style','pushbutton','String','Play from here','Position',[785,16,110,28], ...
    'Callback', @(~,~) playFromHere());

% AppData defaults
setappdata(fig,'isPaused',false);
setappdata(fig,'speedFactor',1);
setappdata(fig,'currentIndex',1);
setappdata(fig,'isPlaying',false);

% Listeners to keep readouts synced and allow live preview while scrubbing
addlistener(speedSlider,'Value','PostSet', @(~,ev) set(speedReadout,'String',sprintf('%.2fx', ev.AffectedObject.Value)));
addlistener(scrubSlider,'Value','PostSet', @(~,ev) onScrub(ev.AffectedObject.Value));

% --- Precompute index sequence for animation ---
nSamp = numel(T_angle);
idxSeqAll = 1:nSamp;  % logical time index space

% --- Initial frame ---
updateFrame(1);

% --- Helper: handle scrubbing (show frame immediately, don't start play) ---
    function onScrub(tVal)
        % find nearest sample index for this time
        [~, idx] = min(abs(T_angle - tVal));
        set(scrubReadout,'String',sprintf('%.2fs', T_angle(idx)));
        setappdata(fig,'currentIndex', idx);
        updateFrame(idx);  % preview frame instantly
    end

% --- Helper: play from the current scrubbed time ---
    function playFromHere()
        if ~isvalid(fig), return; end
        setappdata(fig,'isPlaying', true);
        k0 = getappdata(fig,'currentIndex');
        if isempty(k0) || k0 < 1, k0 = 1; end
        % build decimated sequence starting from k0
        idxSeq = k0:stepK:nSamp;
        for ii = 1:numel(idxSeq)
            if ~isvalid(fig), break; end
            % allow Pause
            while getappdata(fig,'isPaused'); pause(0.05); end
            % stop if user starts scrubbing or toggles play state
            if ~getappdata(fig,'isPlaying'), break; end

            k = idxSeq(ii);
            setappdata(fig,'currentIndex', k);
            updateFrame(k);

            % pacing
            sf = getappdata(fig,'speedFactor'); if isempty(sf), sf = 1; end
            pause(max(0, basePause)/max(0.01, sf));
        end
        setappdata(fig,'isPlaying', false);
    end

% --- Helper: draw a given frame index k (sync both panels) ---
    function updateFrame(k)
        k = max(1, min(k, numel(T_angle)));
        % Left panel
        kk = min(k+1, numel(X_pos)); % guard for path indexing
        set(hPath, 'XData', X_pos(1:kk), 'YData', Y_pos(1:kk));
        set(hDot,  'XData', X_pos(kk),   'YData', Y_pos(kk));
        set(trajLabel, 'String', sprintf('Time: %.2f s', T_angle(k)));

        % Right panel
        tk = T_angle(k);
        set(hNow, 'XData', [tk tk]);
        set(angLabel, 'String', sprintf('Time: %.2f s', tk));

        drawnow limitrate
    end

end
