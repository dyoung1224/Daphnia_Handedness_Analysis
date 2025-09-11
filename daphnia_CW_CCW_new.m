%% NOTE: THIS IS THE NEW CODE THAT HAS NEW ANIMATION FEATURES (PAUSE, SKIP, TIME RANGE, etc.), BUT THE CW/CCW ANALYSIS IS NOT CORRECT, AND IT DOESN'T MARK NEUTRAL TIME PERIODS.

function daphnia_CW_CCW
% clc; close all; clearvars;

%% ================== USER SETTINGS ==================
pre  = 'C:\Users\swany\Videos\Daphnia\6-23,24-25\';
base = 'Multiple_6-23-24-25';
csvs = '_csv';
inputDir = fullfile(pre, [base csvs]);

d = 3;                        % 1-based selection (file uses d-1)
csvFile = fullfile(inputDir, sprintf('%s_daphnia%d.csv', base, d-1));

binSize    = 0.5;             % speed histogram bin size (cm/s)
timeWindow = 1.0;             % window size for CW/CCW (s)

stepK     = 2;                % decimation for animation (bigger = faster, fewer frames)
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

% Positions aligned to Time_vel
X_pos_full = X(1:end-1); Y_pos_full = Y(1:end-1);
X_pos_full = X_pos_full(valid_indices);
Y_pos_full = Y_pos_full(valid_indices);

if numel(T_f) < 3, error('Too few samples after v_* filter.'); end

%% -------- Cumulative angle (turns) over FULL series --------
dVx = diff(Vx_f); dVy = diff(Vy_f);
Vx_k = Vx_f(1:end-1); Vy_k = Vy_f(1:end-1);
den = Vx_k.^2 + Vy_k.^2; den(den==0) = eps;
delta_phi = (Vx_k .* dVy - Vy_k .* dVx) ./ den; % radians
norm_delta_phi = delta_phi ./ (2*pi);           % turns
CA_full = cumsum(norm_delta_phi);               % cumulative angle (turns)
T_ang_full = T_f(1:end-1);                      % time for CA_full

% Full time bounds
t0_full = min(T_ang_full);
t1_full = max(T_ang_full);

%% -------- UI FIGURE & STATIC OBJECTS --------
fig = figure('Name','Daphnia Handedness (Range + Speed + Scrub)','Color','w','Renderer','opengl','Position',[100 100 1200 650]);
tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');

% Left: Trajectory
axL = nexttile; hold(axL,'on'); axis(axL,'equal'); box(axL,'on'); grid(axL,'on');
xlabel(axL,'X (cm)'); ylabel(axL,'Y (cm)'); title(axL,'2D trajectory');

% Placeholders (updated on range apply)
hFullPath = plot(axL, nan, nan, '-', 'Color',[0.8 0.8 0.8], 'LineWidth',1);
hPath     = plot(axL, nan, nan, '-', 'LineWidth', 2);
hDot      = plot(axL, nan, nan, 'o', 'MarkerSize', 6, 'MarkerFaceColor','k', 'MarkerEdgeColor','k');
trajLabel = text(axL, 0,0,'', 'FontSize',11, 'FontWeight','bold', ...
    'BackgroundColor','w','EdgeColor','k','Margin',3);

% Right: bands + cumulative angle
axR = nexttile; hold(axR,'on'); box(axR,'on'); grid(axR,'on');
xlabel(axR,'Time (s)'); ylabel(axR,'Cumulative angle (turns)');
title(axR,'CW/CCW windows + cumulative angle');

% Stripe image (bands) and line; created once, data swapped on range apply
hStripe = imagesc(axR, [t0_full t1_full], [0 1], 2); % temp
colormap(axR, [0.8 0.8 1; 1 1 1; 1 0.8 0.8]); % [blue; white; red]
set(axR,'YDir','normal'); alpha(hStripe, 0.4);
hAngle = plot(axR, nan, nan, 'k', 'LineWidth', 1.2);
hNow   = line(axR, [t0_full t0_full], [0 1], 'Color',[0.2 0.2 0.2], 'LineStyle',':');
angLabel = text(axR, 0,0,'', 'FontSize',11,'FontWeight','bold', ...
    'BackgroundColor','w','EdgeColor','k','Margin',3);

%% -------- CONTROLS (bottom) --------
% Pause/Resume
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

% Scrub slider (seek by time)
uicontrol('Style','text','String','Scrub', 'Position',[415,20,40,18], ...
    'HorizontalAlignment','left','BackgroundColor',get(fig,'Color'));
scrubSlider = uicontrol('Style','slider','Min',t0_full,'Max',t1_full,'Value',t0_full, ...
    'Position',[455,20,220,18]);
scrubReadout = uicontrol('Style','text','String',sprintf('%.2fs',t0_full), ...
    'Position',[680,20,60,18],'HorizontalAlignment','left', ...
    'BackgroundColor',get(fig,'Color'));
uicontrol('Style','pushbutton','String','Play from here','Position',[745,16,110,28], ...
    'Callback', @(~,~) playFromHere());

% ---- Time Range controls ----
uicontrol('Style','text','String','Range Start', 'Position',[880,36,70,18], ...
    'HorizontalAlignment','left','BackgroundColor',get(fig,'Color'));
rangeStartSlider = uicontrol('Style','slider','Min',t0_full,'Max',t1_full,'Value',t0_full, ...
    'Position',[950,36,200,18]);

uicontrol('Style','text','String','Range End', 'Position',[880,16,70,18], ...
    'HorizontalAlignment','left','BackgroundColor',get(fig,'Color'));
rangeEndSlider = uicontrol('Style','slider','Min',t0_full,'Max',t1_full,'Value',t1_full, ...
    'Position',[950,16,200,18]);

rangeStartReadout = uicontrol('Style','text','String',sprintf('%.2fs',t0_full), ...
    'Position',[1155,36,60,18],'HorizontalAlignment','left','BackgroundColor',get(fig,'Color'));
rangeEndReadout   = uicontrol('Style','text','String',sprintf('%.2fs',t1_full), ...
    'Position',[1155,16,60,18],'HorizontalAlignment','left','BackgroundColor',get(fig,'Color'));

uicontrol('Style','pushbutton','String','Apply Range','Position',[1040,58,110,26], ...
    'Callback', @(~,~) applyRange());

% AppData defaults
setappdata(fig,'isPaused',false);
setappdata(fig,'speedFactor',1);
setappdata(fig,'currentIndex',1);
setappdata(fig,'isPlaying',false);

% Listeners for readouts / instant scrub-preview
addlistener(speedSlider,'Value','PostSet', @(~,ev) set(speedReadout,'String',sprintf('%.2fx', ev.AffectedObject.Value)));
addlistener(scrubSlider,'Value','PostSet', @(~,ev) onScrub(ev.AffectedObject.Value));
addlistener(rangeStartSlider,'Value','PostSet', @(~,ev) onRangeSlider('start', ev.AffectedObject.Value));
addlistener(rangeEndSlider,'Value','PostSet', @(~,ev) onRangeSlider('end',   ev.AffectedObject.Value));

%% -------- INITIAL RANGE = FULL --------
state = struct();
applyRange();             % computes state for full range and draws first frame
updateFrame(1);           % show first frame

%% =================== NESTED HELPERS ===================
    function onRangeSlider(which, val)
        % Clamp so start < end with small safety margin
        s = get(rangeStartSlider,'Value'); e = get(rangeEndSlider,'Value');
        epsT = 1e-9;
        if strcmp(which,'start')
            s = min(val, e - epsT);
            set(rangeStartSlider,'Value', s);
        else
            e = max(val, s + epsT);
            set(rangeEndSlider,'Value', e);
        end
        set(rangeStartReadout,'String',sprintf('%.2fs', get(rangeStartSlider,'Value')));
        set(rangeEndReadout,  'String',sprintf('%.2fs', get(rangeEndSlider,'Value')));
    end

    function applyRange()
        % Stop playing while recomputing
        setappdata(fig,'isPlaying',false);

        % Read range, clamp to full bounds
        tA = max(t0_full, get(rangeStartSlider,'Value'));
        tB = min(t1_full, get(rangeEndSlider,'Value'));
        if tB <= tA, tB = tA + eps; end

        % Subset cumulative angle to [tA,tB]
        maskCA = (T_ang_full >= tA) & (T_ang_full <= tB);
        if nnz(maskCA) < 2
            % if too small, expand slightly
            [~,iA] = min(abs(T_ang_full - tA));
            [~,iB] = min(abs(T_ang_full - tB));
            iA = max(1,iA); iB = max(iA+1, min(numel(T_ang_full), iB));
            maskCA = false(size(T_ang_full)); maskCA(iA:iB) = true;
        end
        T_angle = T_ang_full(maskCA);
        CA = CA_full(maskCA);

        % Subset positions aligned to T_angle (X_pos_full/Y_pos_full are 1 longer)
        idxA = find(maskCA,1,'first');
        idxB = find(maskCA,1,'last');
        X_pos = X_pos_full(idxA:idxB+1);
        Y_pos = Y_pos_full(idxA:idxB+1);

        % Compute slopeSigns over the selected range
        slopeSigns = computeSlopeSigns(T_angle, CA, timeWindow);

        % Update state struct shared by animator
        state.T_angle = T_angle;
        state.CA = CA;
        state.X_pos = X_pos;
        state.Y_pos = Y_pos;
        state.slopeSigns = slopeSigns;

        % Axes limits
        pad = 0.05;
        xlim(axL, [min(X_pos)-pad*range(X_pos), max(X_pos)+pad*range(X_pos)]);
        ylim(axL, [min(Y_pos)-pad*range(Y_pos), max(Y_pos)+pad*range(Y_pos)]);
        set(trajLabel, 'Position', [axL.XLim(1)+0.05*range(axL.XLim), axL.YLim(2)-0.08*range(axL.YLim), 0]);

        xlim(axR, [min(T_angle), max(T_angle)]);
        ylim(axR, [min(CA)-1, max(CA)+1]);
        set(angLabel, 'Position', [axR.XLim(1)+0.05*range(axR.XLim), axR.YLim(2)-0.08*range(axR.YLim), 0]);

        % Update full path (light), reset animated path
        set(hFullPath, 'XData', X_pos, 'YData', Y_pos);
        set(hPath,     'XData', nan,   'YData', nan);
        set(hDot,      'XData', nan,   'YData', nan);

        % Update cumulative angle line
        set(hAngle, 'XData', T_angle, 'YData', CA);

        % Update stripe bands (imagesc)
        stripe = slopeSigns;                 % -1,0,+1
        stripe(stripe==-1) = 1;              % 1 -> blue (CW)
        stripe(stripe== 0) = 2;              % 2 -> white (flat)
        stripe(stripe==+1) = 3;              % 3 -> red (CCW)
        % Use a 1xN image stretched to axes limits
        set(hStripe, 'XData', [min(T_angle), max(T_angle)], ...
                     'YData', [min(CA)-1, max(CA)+1], ...
                     'CData', stripe);

        % Reset scrub slider to this range
        set(scrubSlider,'Min', min(T_angle), 'Max', max(T_angle), 'Value', min(T_angle));
        set(scrubReadout,'String',sprintf('%.2fs', min(T_angle)));

        % Move "now" line to range start
        set(hNow, 'XData', [min(T_angle) min(T_angle)], 'YData', ylim(axR));

        % Reset play state/index
        setappdata(fig,'currentIndex', 1);
        drawnow;
    end

    function slopeSigns = computeSlopeSigns(T_angle, CA, dt)
        t0 = min(T_angle); t1 = max(T_angle);
        numW = max(1, floor((t1 - t0) / dt));
        slopeSigns = zeros(1, numW);

        binIdx = ceil((T_angle - t0) / dt);
        binIdx(binIdx < 1) = 1; binIdx(binIdx > numW) = numW;

        for w = 1:numW
            mask = (binIdx == w) & isfinite(T_angle) & isfinite(CA);
            if nnz(mask) < 2, slopeSigns(w) = 0; continue; end
            xw = T_angle(mask); yw = CA(mask);
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
    end

    function onScrub(tVal)
        T_angle = state.T_angle;
        [~, idx] = min(abs(T_angle - tVal));
        set(scrubReadout,'String',sprintf('%.2fs', T_angle(idx)));
        setappdata(fig,'currentIndex', idx);
        updateFrame(idx);  % preview
    end

    function playFromHere()
        if ~isvalid(fig), return; end
        setappdata(fig,'isPlaying', true);
        T_angle = state.T_angle;
        k0 = getappdata(fig,'currentIndex');
        if isempty(k0) || k0 < 1, k0 = 1; end

        nSamp = numel(T_angle);
        idxSeq = k0:stepK:nSamp;
        for ii = 1:numel(idxSeq)
            if ~isvalid(fig), break; end
            while getappdata(fig,'isPaused'); pause(0.05); end
            if ~getappdata(fig,'isPlaying'), break; end
            k = idxSeq(ii);
            setappdata(fig,'currentIndex', k);
            updateFrame(k);
            sf = getappdata(fig,'speedFactor'); if isempty(sf), sf = 1; end
            pause(max(0, basePause)/max(0.01, sf));
        end
        setappdata(fig,'isPlaying', false);
    end

    function updateFrame(k)
        T_angle = state.T_angle;
        CA      = state.CA;
        X_pos   = state.X_pos;
        Y_pos   = state.Y_pos;

        if isempty(T_angle), return; end
        k = max(1, min(k, numel(T_angle)));

        % Left panel
        kk = min(k+1, numel(X_pos)); % X_pos/Y_pos are one element longer
        set(hPath, 'XData', X_pos(1:kk), 'YData', Y_pos(1:kk));
        set(hDot,  'XData', X_pos(kk),   'YData', Y_pos(kk));
        set(trajLabel, 'String', sprintf('Time: %.2f s', T_angle(k)));

        % Right panel
        set(hAngle, 'XData', T_angle, 'YData', CA);
        set(hNow,   'XData', [T_angle(k) T_angle(k)], 'YData', ylim(axR));
        set(angLabel, 'String', sprintf('Time: %.2f s', T_angle(k)));

        drawnow limitrate
    end
end


