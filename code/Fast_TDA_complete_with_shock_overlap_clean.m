%% FAST TDA COMPLETE ANALYSIS
%
% Input columns:
%   Time, Amplitude, Shock
%
% Main outputs:
%   1) Signal + moving TDA window
%   2) Delay-embedded point cloud + fitted ellipse
%   3) Original conic coefficients
%   4) Normalised conic coefficients
%   5) Zoomed normalised coefficients (3 to 4 s)
%   6) Normalised absolute coefficients (stacked)
%   7) Normalised absolute coefficient response (overlap)
%   8) Normalised absolute coefficient response + pure shock (overlap)
%   9) Raw geometric descriptors
%  10) Normalised geometric descriptors
%  11) CSV + MAT files with numerical results

clear; clc; close all;

%% 1. READ DATA
input_file = 'AMP_10-VIB_0.4-SH_0.4.csv';
data = readtable(input_file);

required_columns = {'Time','Amplitude','Shock'};
for k = 1:numel(required_columns)
    if ~ismember(required_columns{k}, data.Properties.VariableNames)
        error('Required column "%s" was not found.', required_columns{k});
    end
end

time   = data.Time(:);
output = data.Amplitude(:);
shock  = data.Shock(:);

if numel(time) ~= numel(output) || numel(time) ~= numel(shock)
    error('Time, Amplitude, and Shock must have the same number of samples.');
end

valid = isfinite(time) & isfinite(output) & isfinite(shock);
time = time(valid);
output = output(valid);
shock = shock(valid);

fprintf('Data loaded: %d rows\n', numel(time));

%% 2. SAMPLING RATE
dt = median(diff(time));
Fs = 1/dt;
N = numel(output);

fprintf('Sampling frequency = %.4f Hz\n', Fs);
fprintf('Sampling interval  = %.8f s\n', dt);
fprintf('Signal start       = %.4f s\n', time(1));
fprintf('Signal end         = %.4f s\n', time(end));
fprintf('Signal duration    = %.4f s\n', time(end)-time(1));

%% 3. FFT
Y = fft(output);
Ym = abs(Y/N);

if mod(N,2) == 0
    f = (0:N/2)'*(Fs/N);
    Ym = Ym(1:N/2+1);
    Ym(2:end-1) = 2*Ym(2:end-1);
else
    f = (0:(N-1)/2)'*(Fs/N);
    Ym = Ym(1:(N+1)/2);
    Ym(2:end) = 2*Ym(2:end);
end

%% 4. FIND SPECTRAL PEAKS
max_mag = max(Ym);
[~,pk] = findpeaks(Ym, ...
    'MinPeakHeight',0.05*max_mag, ...
    'MinPeakProminence',0.10*max_mag, ...
    'MinPeakDistance',5);

if isempty(pk)
    error('No spectral peaks found. Adjust findpeaks thresholds.');
end

peak_freqs = f(pk);
Fmin = min(peak_freqs);
Fmax = max(peak_freqs);

fprintf('\nFmax = %.4f Hz\n', Fmax);
fprintf('Fmin = %.4f Hz\n', Fmin);

%% 5. FAST TDA PARAMETERS
time_delay = 0.25/Fmax;
time_delay_samp = max(round(time_delay/dt),1);

win_pts = round((1/Fmin)/dt);
win_pts = max(win_pts,time_delay_samp+5);
win_dur = (win_pts-1)*dt;
num_windows = N-win_pts+1;
step_size = 5;

time_assignment = 'end';

if num_windows < 1
    error('Window size is too large for the signal.');
end

fprintf('\nTDA parameters\n');
fprintf('---------------------------------\n');
fprintf('Window size     = %d samples\n', win_pts);
fprintf('Window duration = %.6f s\n', win_dur);
fprintf('Time delay      = %d samples\n', time_delay_samp);
fprintf('Time delay      = %.6f s\n', time_delay);
fprintf('Step size       = %d samples\n', step_size);

%% 6. PREALLOCATE
all_params = nan(6,num_windows);
all_desc = nan(12,num_windows);
valid_idx = [];

%% 7. PRE-COMPUTE ELLIPSES AND DESCRIPTORS
fprintf('\nPre-computing ellipses and descriptors...\n');

for i = 1:step_size:num_windows
    i_end = i + win_pts - 1 - time_delay_samp;
    i_del = i + time_delay_samp;
    i_del_end = i_del + win_pts - 1 - time_delay_samp;

    if i_end > N || i_del_end > N
        break;
    end

    P = [output(i:i_end), output(i_del:i_del_end)];
    params = fit_ellipse_hls_golden(P(:,1),P(:,2));

    if any(~isfinite(params)) || all(params == 0)
        continue;
    end

    all_params(:,i) = params;
    all_desc(:,i) = ellipse_descriptors(params);
    valid_idx(end+1) = i; %#ok<AGROW>
end

if isempty(valid_idx)
    error('No valid ellipse windows were found.');
end

fprintf('Pre-computation complete: %d valid windows.\n',numel(valid_idx));

%% 8. WINDOW TIME INFORMATION
n_valid = numel(valid_idx);
window_number = (1:n_valid)';
sample_index = valid_idx(:);

t_start = time(valid_idx);
end_indices = min(valid_idx(:)+win_pts-1,N);
t_end = time(end_indices);
t_center = 0.5*(t_start+t_end);

t_start = t_start(:);
t_end = t_end(:);
t_center = t_center(:);

switch lower(time_assignment)
    case 'start'
        t_plot = t_start;
    case 'center'
        t_plot = t_center;
    case 'end'
        t_plot = t_end;
    otherwise
        error('time_assignment must be start, center, or end.');
end

fprintf('\nTime assignment for plots = %s\n',time_assignment);
fprintf('First valid window: start %.6f | center %.6f | end %.6f s\n', ...
    t_start(1),t_center(1),t_end(1));
fprintf('Last valid window:  start %.6f | center %.6f | end %.6f s\n', ...
    t_start(end),t_center(end),t_end(end));

%% 9. COEFFICIENT ARRAYS
coeff = all_params(:,valid_idx);
coeff_abs = abs(coeff);

coeff_norm = nan(size(coeff));
coeff_abs_norm = nan(size(coeff_abs));
for j = 1:6
    coeff_norm(j,:) = safe_minmax(coeff(j,:));
    coeff_abs_norm(j,:) = safe_minmax(coeff_abs(j,:));
end

%% 10. BASELINE-REFERENCED ABSOLUTE RESPONSE
% Used for response-timing plots. This does not replace coeff_abs_norm.
response_baseline_duration = 1.0; % s
baseline_mask = t_plot <= t_plot(1)+response_baseline_duration;

if ~any(baseline_mask)
    error('No valid TDA windows are available in the baseline interval.');
end

coeff_abs_response = nan(size(coeff_abs));
coeff_abs_response_norm = nan(size(coeff_abs));

for pidx = 1:6
    v_abs = coeff_abs(pidx,:);
    baseline_value = median(v_abs(baseline_mask),'omitnan');
    v_response = abs(v_abs-baseline_value);
    coeff_abs_response(pidx,:) = v_response;

    max_response = max(v_response,[],'omitnan');
    if ~isfinite(max_response) || max_response < 1e-12
        coeff_abs_response_norm(pidx,:) = zeros(size(v_response));
    else
        coeff_abs_response_norm(pidx,:) = v_response/max_response;
    end
end

% Pure shock magnitude, normalised independently to [0,1].
shock_abs = abs(shock);
shock_norm = safe_minmax(shock_abs);

%% 11. GEOMETRIC DESCRIPTORS
desc = all_desc(:,valid_idx);
angle_unwrapped_rad = unwrap_orientation_pi(desc(9,:));
angle_unwrapped_deg = rad2deg(angle_unwrapped_rad);

desc_plot = [ ...
    desc(1,:); ...       % Center X
    desc(2,:); ...       % Center Y
    desc(3,:); ...       % R major
    desc(4,:); ...       % R minor
    desc(7,:); ...       % Axis ratio
    desc(8,:); ...       % Eccentricity
    angle_unwrapped_deg; ... % Orientation
    desc(11,:); ...      % Ellipse area
    desc(12,:)];         % Focal distance

desc_norm = nan(size(desc_plot));
for j = 1:size(desc_plot,1)
    desc_norm(j,:) = safe_minmax(desc_plot(j,:));
end

%% 12. COLOURS AND STACK SETTINGS
stdColors = [ ...
    0.0000 0.4470 0.7410; ...
    0.8500 0.3250 0.0980; ...
    0.9290 0.6940 0.1250; ...
    0.4940 0.1840 0.5560; ...
    0.4660 0.6740 0.1880];
param_colors = stdColors;
param_names = {'e','d','c','b','a'};
stack_offsets = (0:4)'*1.15;
stack_centers = stack_offsets+0.5;

%% 13. MAIN INTERACTIVE FIGURE
min_acc = min(output);
max_acc = max(output);

fig = figure( ...
    'Name','Fast TDA - Interactive Window Inspector', ...
    'Position',[40 20 1500 1000], ...
    'Color',[0.96 0.96 0.96]);

uicontrol(fig, ...
    'Style','pushbutton', ...
    'Units','normalized', ...
    'Position',[0.88 0.955 0.10 0.03], ...
    'String','Save SVG', ...
    'FontSize',9, ...
    'Callback',@(~,~) cb_save_svg(fig));

%% PLOT 1: FULL SIGNAL + MOVING WINDOW
ax_sig = axes('Parent',fig,'Position',[0.06 0.78 0.88 0.16]);
plot(ax_sig,time,output,'Color',[0.4 0.4 0.4],'LineWidth',0.8);
hold(ax_sig,'on');
win_patch = patch(ax_sig, ...
    [t_start(1) t_end(1) t_end(1) t_start(1)], ...
    [min_acc min_acc max_acc max_acc], ...
    [0.25 0.55 0.95], ...
    'FaceAlpha',0.25, ...
    'EdgeColor',[0.1 0.3 0.8], ...
    'LineWidth',1.5);
xlabel(ax_sig,'Time (s)');
ylabel(ax_sig,'Amplitude');
title(ax_sig,'Signal - drag slider below to move window');
xlim(ax_sig,[time(1) time(end)]);
ylim(ax_sig,[min_acc-0.05*abs(min_acc), max_acc+0.05*abs(max_acc)]);
grid(ax_sig,'on'); box(ax_sig,'on');
hold(ax_sig,'off');

%% PLOT 2: DELAY EMBEDDING + ELLIPSE
ax_ell = axes('Parent',fig,'Position',[0.06 0.47 0.38 0.23]);
xlabel(ax_ell,'Amplitude');
ylabel(ax_ell,'Delayed amplitude');
title(ax_ell,'Ellipse at selected window');
grid(ax_ell,'on'); box(ax_ell,'on'); axis(ax_ell,'equal');

%% PLOT 3: ORIGINAL CONIC COEFFICIENTS
ax_par = axes('Parent',fig,'Position',[0.54 0.47 0.40 0.23]);
hold(ax_par,'on');
for j = 1:5
    pidx = 6-j;
    plot(ax_par,t_plot,coeff(pidx,:), ...
        'Color',param_colors(pidx,:), ...
        'LineWidth',1.2, ...
        'DisplayName',param_names{j});
end
cur_line = xline(ax_par,t_plot(1),'--k','LineWidth',1.5,'HandleVisibility','off');
xlabel(ax_par,'Time (s)');
ylabel(ax_par,'Parameter value');
title(ax_par,'Conic parameters over time');
legend(ax_par,'Location','best');
grid(ax_par,'on'); box(ax_par,'on');
xlim(ax_par,[time(1) time(end)]);
hold(ax_par,'off');

%% PLOT 4: NORMALISED COEFFICIENTS
ax_stk = axes('Parent',fig,'Position',[0.54 0.13 0.40 0.24]);
hold(ax_stk,'on');
for j = 1:5
    pidx = 6-j;
    plot(ax_stk,t_plot,coeff_norm(pidx,:)+stack_offsets(j), ...
        'Color',param_colors(pidx,:),'LineWidth',1.2);
end
cur_line_stk = xline(ax_stk,t_plot(1),'--k','LineWidth',1.5,'HandleVisibility','off');
xlabel(ax_stk,'Time (s)');
ylabel(ax_stk,'Normalised value');
title(ax_stk,'Normalised individual trends (stacked)');
grid(ax_stk,'on'); box(ax_stk,'on');
set(ax_stk,'YTick',stack_centers,'YTickLabel',param_names,'FontWeight','bold');
xlim(ax_stk,[time(1) time(end)]); ylim(ax_stk,[-0.2 5.8]);
hold(ax_stk,'off');

%% PLOT 5: ZOOMED NORMALISED COEFFICIENTS, 3 TO 4 s
ax_stk_zoom = axes('Parent',fig,'Position',[0.06 0.13 0.38 0.24]);
hold(ax_stk_zoom,'on');
mask_zoom = t_plot >= 3 & t_plot <= 4;
if any(mask_zoom)
    for j = 1:5
        pidx = 6-j;
        plot(ax_stk_zoom,t_plot(mask_zoom), ...
            coeff_norm(pidx,mask_zoom)+stack_offsets(j), ...
            'Color',param_colors(pidx,:),'LineWidth',1.2);
    end
    xlabel(ax_stk_zoom,'Time (s)');
    ylabel(ax_stk_zoom,'Normalised value');
    title(ax_stk_zoom,'Normalised individual trends from 3 to 4 s');
    set(ax_stk_zoom,'YTick',stack_centers,'YTickLabel',param_names,'FontWeight','bold');
    xlim(ax_stk_zoom,[3 4]); ylim(ax_stk_zoom,[-0.2 5.8]);
    grid(ax_stk_zoom,'on'); box(ax_stk_zoom,'on');
else
    text(ax_stk_zoom,0.5,0.5,'No data between 3 and 4 s', ...
        'Units','normalized','HorizontalAlignment','center');
end
hold(ax_stk_zoom,'off');

%% 14. SLIDER AND INTERACTIVE STATE
n_wins = n_valid;
uicontrol(fig,'Style','text','Units','normalized', ...
    'Position',[0.01 0.005 0.09 0.025], ...
    'String','Window:','HorizontalAlignment','right', ...
    'BackgroundColor',fig.Color);

sld = uicontrol(fig,'Style','slider','Units','normalized', ...
    'Position',[0.11 0.005 0.65 0.025], ...
    'Min',1,'Max',n_wins,'Value',1, ...
    'SliderStep',[1/max(n_wins-1,1), min(10/max(n_wins-1,1),1)]);

t_lbl = uicontrol(fig,'Style','text','Units','normalized', ...
    'Position',[0.77 0.005 0.20 0.025], ...
    'String',sprintf('start %.3f | center %.3f | end %.3f s',t_start(1),t_center(1),t_end(1)), ...
    'HorizontalAlignment','left','BackgroundColor',fig.Color);

sld.Callback = @(~,~) cb_update(fig);

ud.valid_idx = valid_idx;
ud.all_params = all_params;
ud.output = output;
ud.win_pts = win_pts;
ud.time_delay_samp = time_delay_samp;
ud.t_start = t_start;
ud.t_center = t_center;
ud.t_end = t_end;
ud.t_plot = t_plot;
ud.min_acc = min_acc;
ud.max_acc = max_acc;
ud.ax_sig = ax_sig;
ud.ax_ell = ax_ell;
ud.ax_par = ax_par;
ud.ax_stk = ax_stk;
ud.win_patch = win_patch;
ud.cur_line = cur_line;
ud.cur_line_stk = cur_line_stk;
ud.sld = sld;
ud.t_lbl = t_lbl;
ud.stdColors = stdColors;
fig.UserData = ud;

cb_update(fig);

%% ================================================================
% 21. BASELINE-ALIGNED NORMALISED ABSOLUTE COEFFICIENT PLOT
%     (FIGURE 2 - STACKED)
%
% Purpose:
%
%   Show the absolute coefficient responses individually while
%   placing the normal pre-impact baseline directly on each
%   coefficient tick.
%
%   The original Figure 2 used:
%
%       coeff_abs_norm
%
%   which is a GLOBAL min-max normalisation. Therefore, zero
%   corresponded to the global minimum over the entire record,
%   not to the normal pre-impact state.
%
%   Here we use:
%
%       coeff_abs_response_norm
%
%   from Section 10, which is referenced to the pre-impact
%   baseline.
%
% ================================================================

fig_abs_norm = figure( ...
    'Name', ...
    'Normalised Absolute Ellipse Parameters', ...
    'Position', ...
    [180 120 1050 650], ...
    'Color', ...
    'white');

ax_abs_norm = ...
    axes('Parent',fig_abs_norm);

hold(ax_abs_norm,'on');


%% ================================================================
% 21.1 STACKING SETTINGS
% ================================================================

% Each plotted response is approximately in [0,1].
%
% A spacing of 1.6 gives a clear gap between neighbouring
% coefficient trends.

abs_stack_spacing = 1.6;

abs_tick_positions = ...
    (0:4)' * abs_stack_spacing;


%% ================================================================
% 21.2 PLOT BASELINE-ALIGNED RESPONSES
% ================================================================

for j = 1:5

    % ------------------------------------------------------------
    % Coefficient ordering:
    %
    % j = 1 -> e
    % j = 2 -> d
    % j = 3 -> c
    % j = 4 -> b
    % j = 5 -> a
    % ------------------------------------------------------------

    pidx = 6-j;


    % ------------------------------------------------------------
    % Take the baseline-referenced response calculated in
    % Section 10.
    % ------------------------------------------------------------

    v_response = ...
        coeff_abs_response_norm(pidx,:);


    % ------------------------------------------------------------
    % Estimate the remaining small baseline level.
    %
    % The median is used rather than the first point so that
    % normal baseline fluctuations do not determine the
    % vertical alignment.
    % ------------------------------------------------------------

    baseline_level = ...
        median( ...
            v_response(baseline_mask), ...
            'omitnan');


    % ------------------------------------------------------------
    % Shift the baseline to zero.
    % ------------------------------------------------------------

    v_response = ...
        v_response - baseline_level;


    % ------------------------------------------------------------
    % Prevent the normal baseline from falling below its tick.
    %
    % Small negative values here simply mean that the signal
    % is below the median baseline level.
    % For the stacked visualisation, these are clipped at zero.
    % ------------------------------------------------------------

    v_response = ...
        max(v_response,0);


    % ------------------------------------------------------------
    % Re-normalise the displayed response to [0,1].
    %
    % This is ONLY for the figure.
    % The saved numerical data are not changed.
    % ------------------------------------------------------------

    response_max = ...
        max(v_response,[],'omitnan');


    if isfinite(response_max) && ...
            response_max > 1e-12

        v_response = ...
            v_response / response_max;

    else

        v_response = ...
            zeros(size(v_response));

    end


    % ------------------------------------------------------------
    % Add vertical stacking offset.
    %
    % The baseline now corresponds directly to the tick.
    % ------------------------------------------------------------

    v_plot = ...
        v_response + abs_tick_positions(j);


    % ------------------------------------------------------------
    % Plot
    % ------------------------------------------------------------

    plot( ...
        ax_abs_norm, ...
        t_plot, ...
        v_plot, ...
        'Color', ...
        param_colors(pidx,:), ...
        'LineWidth', ...
        1.2);

end


%% ================================================================
% 21.3 AXIS LABELS
% ================================================================

xlabel( ...
    ax_abs_norm, ...
    'Time (s)');

ylabel( ...
    ax_abs_norm, ...
    'Normalised absolute response');


title( ...
    ax_abs_norm, ...
    'Normalised Absolute Individual Trends (stacked)');


%% ================================================================
% 21.4 Y-AXIS TICKS
%
% Each tick represents the normal pre-impact baseline level.
% ================================================================

set( ...
    ax_abs_norm, ...
    'YTick', ...
    abs_tick_positions, ...
    'YTickLabel', ...
    param_names, ...
    'FontWeight', ...
    'bold');


%% ================================================================
% 21.5 AXIS LIMITS
% ================================================================

xlim( ...
    ax_abs_norm, ...
    [time(1) time(end)]);

ylim( ...
    ax_abs_norm, ...
    [ ...
    -0.2, ...
    abs_tick_positions(end) + 1.4 ...
    ]);


%% ================================================================
% 21.6 GRID / BOX
% ================================================================

grid( ...
    ax_abs_norm, ...
    'on');

box( ...
    ax_abs_norm, ...
    'on');

hold( ...
    ax_abs_norm, ...
    'off');


%% ================================================================
% 21.A BASELINE-REFERENCED NORMALISED ABSOLUTE COEFFICIENTS
%     — OVERLAP PLOT
%
% Purpose:
%
%   All five coefficient responses are placed on one common
%   0-to-1 scale.
%
%   Baseline:
%       approximately 0
%
%   Impact:
%       departure from 0
% ================================================================

fig_abs_overlap = figure( ...
    'Name', ...
    'Normalised Absolute Coefficient Response - Overlap', ...
    'Position', ...
    [220 120 1050 650], ...
    'Color', ...
    'white');

ax_abs_overlap = ...
    axes('Parent',fig_abs_overlap);

hold(ax_abs_overlap,'on');


for j = 1:5

    % ------------------------------------------------------------
    % e, d, c, b, a
    % ------------------------------------------------------------

    pidx = 6-j;


    % ------------------------------------------------------------
    % Baseline-referenced response
    % ------------------------------------------------------------

    v_response = ...
        coeff_abs_response_norm(pidx,:);


    % ------------------------------------------------------------
    % Align normal baseline approximately to zero.
    % ------------------------------------------------------------

    baseline_level = ...
        median( ...
            v_response(baseline_mask), ...
            'omitnan');


    v_response = ...
        v_response - baseline_level;


    v_response = ...
        max(v_response,0);


    % ------------------------------------------------------------
    % Normalise displayed response
    % ------------------------------------------------------------

    response_max = ...
        max(v_response,[],'omitnan');


    if isfinite(response_max) && ...
            response_max > 1e-12

        v_response = ...
            v_response / response_max;

    else

        v_response = ...
            zeros(size(v_response));

    end


    % ------------------------------------------------------------
    % Plot
    % ------------------------------------------------------------

    plot( ...
        ax_abs_overlap, ...
        t_plot, ...
        v_response, ...
        'Color', ...
        param_colors(pidx,:), ...
        'LineWidth', ...
        1.3, ...
        'DisplayName', ...
        param_names{j});

end


%% ================================================================
% 21.A.1 ZERO BASELINE
% ================================================================

yline( ...
    ax_abs_overlap, ...
    0, ...
    '--k', ...
    'LineWidth', ...
    1.0, ...
    'HandleVisibility', ...
    'off');


%% ================================================================
% 21.A.2 AXIS SETTINGS
% ================================================================

xlabel( ...
    ax_abs_overlap, ...
    'Time (s)');

ylabel( ...
    ax_abs_overlap, ...
    'Normalised absolute response from baseline');

title( ...
    ax_abs_overlap, ...
    'Normalised Absolute Ellipse-Parameter Response (overlap)');

xlim( ...
    ax_abs_overlap, ...
    [time(1) time(end)]);

ylim( ...
    ax_abs_overlap, ...
    [0 1.05]);

legend( ...
    ax_abs_overlap, ...
    'Location', ...
    'best');

grid( ...
    ax_abs_overlap, ...
    'on');

box( ...
    ax_abs_overlap, ...
    'on');

hold( ...
    ax_abs_overlap, ...
    'off');


%% ================================================================
% 21.B BASELINE-REFERENCED NORMALISED ABSOLUTE COEFFICIENTS
%      + PURE SHOCK OVERLAP
%
% Purpose:
%
%   Compare the response timing of:
%
%       a, b, c, d, e
%
%   with the pure shock signal.
%
%   All quantities are displayed on the same 0-to-1 scale.
% ================================================================

fig_abs_shock = figure( ...
    'Name', ...
    'Normalised Absolute Coefficients and Pure Shock', ...
    'Position', ...
    [220 120 1100 650], ...
    'Color', ...
    'white');

ax_abs_shock = ...
    axes('Parent',fig_abs_shock);

hold(ax_abs_shock,'on');


%% ================================================================
% 21.B.1 PLOT COEFFICIENT RESPONSES
% ================================================================

for j = 1:5

    pidx = 6-j;


    % ------------------------------------------------------------
    % Baseline-referenced coefficient response
    % ------------------------------------------------------------

    v_response = ...
        coeff_abs_response_norm(pidx,:);


    % ------------------------------------------------------------
    % Remove the normal baseline level.
    % ------------------------------------------------------------

    baseline_level = ...
        median( ...
            v_response(baseline_mask), ...
            'omitnan');


    v_response = ...
        v_response - baseline_level;


    v_response = ...
        max(v_response,0);


    % ------------------------------------------------------------
    % Normalise displayed response
    % ------------------------------------------------------------

    response_max = ...
        max(v_response,[],'omitnan');


    if isfinite(response_max) && ...
            response_max > 1e-12

        v_response = ...
            v_response / response_max;

    else

        v_response = ...
            zeros(size(v_response));

    end


    % ------------------------------------------------------------
    % Plot coefficient
    % ------------------------------------------------------------

    plot( ...
        ax_abs_shock, ...
        t_plot, ...
        v_response, ...
        'Color', ...
        param_colors(pidx,:), ...
        'LineWidth', ...
        1.25, ...
        'DisplayName', ...
        param_names{j});

end


%% ================================================================
% 21.B.2 PURE SHOCK SIGNAL
% ================================================================
%
% shock_norm is the independently min-max normalised magnitude
% of the Shock column from the input CSV.
%
% This preserves the original shock timing while placing its
% magnitude on a 0-to-1 scale.
% ================================================================

plot( ...
    ax_abs_shock, ...
    time, ...
    shock_norm, ...
    '--k', ...
    'LineWidth', ...
    1.8, ...
    'DisplayName', ...
    'Pure shock');


%% ================================================================
% 21.B.3 ZERO BASELINE
% ================================================================

yline( ...
    ax_abs_shock, ...
    0, ...
    ':k', ...
    'LineWidth', ...
    0.9, ...
    'HandleVisibility', ...
    'off');


%% ================================================================
% 21.B.4 AXIS SETTINGS
% ================================================================

xlabel( ...
    ax_abs_shock, ...
    'Time (s)');

ylabel( ...
    ax_abs_shock, ...
    'Normalised response / shock magnitude');

title( ...
    ax_abs_shock, ...
    'Normalised Absolute Ellipse Parameters and Pure Shock');

xlim( ...
    ax_abs_shock, ...
    [time(1) time(end)]);

ylim( ...
    ax_abs_shock, ...
    [0 1.05]);


legend( ...
    ax_abs_shock, ...
    'Location', ...
    'best');


grid( ...
    ax_abs_shock, ...
    'on');

box( ...
    ax_abs_shock, ...
    'on');

hold( ...
    ax_abs_shock, ...
    'off');

%% 22. RAW GEOMETRIC DESCRIPTORS
fig_desc = figure( ...
    'Name','Geometric Ellipse Descriptors', ...
    'Position',[160 70 1250 850], ...
    'Color','white');
desc_colors = lines(9);
descriptor_names = {'Center X','Center Y','R major','R minor', ...
    'Axis ratio','Eccentricity','Orientation (deg)','Ellipse area','Focal distance'};

for j = 1:9
    ax = subplot(3,3,j);
    plot(ax,t_plot,desc_plot(j,:), ...
        'Color',desc_colors(j,:),'LineWidth',1.2);
    xlabel(ax,'Time (s)');
    ylabel(ax,descriptor_names{j});
    title(ax,descriptor_names{j});
    xlim(ax,[time(1) time(end)]);
    grid(ax,'on'); box(ax,'on');
end
sgtitle(fig_desc,'Geometric Ellipse Descriptors over Time');

%% 23. NORMALISED GEOMETRIC DESCRIPTORS - STACKED
fig_desc_norm = figure( ...
    'Name','Normalised Geometric Descriptors', ...
    'Position',[220 100 1150 750], ...
    'Color','white');
ax_desc_norm = axes('Parent',fig_desc_norm);
hold(ax_desc_norm,'on');

desc_offsets = (0:size(desc_norm,1)-1)'*1.15;
desc_centers = desc_offsets+0.5;

for j = 1:size(desc_norm,1)
    plot(ax_desc_norm,t_plot,desc_norm(j,:)+desc_offsets(j), ...
        'Color',desc_colors(j,:),'LineWidth',1.15);
end

xlabel(ax_desc_norm,'Time (s)');
ylabel(ax_desc_norm,'Normalised descriptor value');
title(ax_desc_norm,'Normalised Geometric Descriptors (stacked)');
set(ax_desc_norm,'YTick',desc_centers,'YTickLabel',descriptor_names,'FontWeight','bold');
xlim(ax_desc_norm,[time(1) time(end)]);
ylim(ax_desc_norm,[-0.2 desc_offsets(end)+1.4]);
grid(ax_desc_norm,'on'); box(ax_desc_norm,'on');
hold(ax_desc_norm,'off');

%% 24. SAVE COEFFICIENT NUMERICAL VALUES
T_coeff = table( ...
    window_number,sample_index,t_start,t_center,t_end, ...
    coeff(1,:)',coeff(2,:)',coeff(3,:)',coeff(4,:)',coeff(5,:)',coeff(6,:)', ...
    coeff_abs(1,:)',coeff_abs(2,:)',coeff_abs(3,:)',coeff_abs(4,:)',coeff_abs(5,:)',coeff_abs(6,:)', ...
    coeff_norm(1,:)',coeff_norm(2,:)',coeff_norm(3,:)',coeff_norm(4,:)',coeff_norm(5,:)',coeff_norm(6,:)', ...
    coeff_abs_norm(1,:)',coeff_abs_norm(2,:)',coeff_abs_norm(3,:)',coeff_abs_norm(4,:)',coeff_abs_norm(5,:)',coeff_abs_norm(6,:)', ...
    coeff_abs_response(1,:)',coeff_abs_response(2,:)',coeff_abs_response(3,:)',coeff_abs_response(4,:)',coeff_abs_response(5,:)',coeff_abs_response(6,:)', ...
    coeff_abs_response_norm(1,:)',coeff_abs_response_norm(2,:)',coeff_abs_response_norm(3,:)',coeff_abs_response_norm(4,:)',coeff_abs_response_norm(5,:)',coeff_abs_response_norm(6,:)', ...
    'VariableNames',{ ...
    'WindowNumber','SampleIndex','TimeStart_s','TimeCenter_s','TimeEnd_s', ...
    'a','b','c','d','e','f', ...
    'abs_a','abs_b','abs_c','abs_d','abs_e','abs_f', ...
    'norm_a','norm_b','norm_c','norm_d','norm_e','norm_f', ...
    'norm_abs_a','norm_abs_b','norm_abs_c','norm_abs_d','norm_abs_e','norm_abs_f', ...
    'response_abs_a','response_abs_b','response_abs_c','response_abs_d','response_abs_e','response_abs_f', ...
    'norm_response_abs_a','norm_response_abs_b','norm_response_abs_c', ...
    'norm_response_abs_d','norm_response_abs_e','norm_response_abs_f'});
writetable(T_coeff,'fast_tda_coefficients.csv');
fprintf('\nSaved coefficient values: fast_tda_coefficients.csv\n');

%% 25. SAVE GEOMETRIC DESCRIPTOR NUMERICAL VALUES
T_desc = table( ...
    window_number,sample_index,t_start,t_center,t_end, ...
    desc(1,:)',desc(2,:)',desc(3,:)',desc(4,:)',desc(5,:)',desc(6,:)', ...
    desc(7,:)',desc(8,:)',desc(9,:)',desc(10,:)',angle_unwrapped_deg', ...
    desc(11,:)',desc(12,:)', ...
    desc_norm(1,:)',desc_norm(2,:)',desc_norm(3,:)',desc_norm(4,:)',desc_norm(5,:)', ...
    desc_norm(6,:)',desc_norm(7,:)',desc_norm(8,:)',desc_norm(9,:)', ...
    'VariableNames',{ ...
    'WindowNumber','SampleIndex','TimeStart_s','TimeCenter_s','TimeEnd_s', ...
    'CenterX','CenterY','R_major','R_minor','MajorAxisLength','MinorAxisLength', ...
    'AxisRatio','Eccentricity','Angle_rad','Angle_deg_wrapped','Angle_deg_unwrapped', ...
    'EllipseArea','FocalDistance', ...
    'norm_CenterX','norm_CenterY','norm_R_major','norm_R_minor', ...
    'norm_AxisRatio','norm_Eccentricity','norm_Angle','norm_EllipseArea','norm_FocalDistance'});
writetable(T_desc,'fast_tda_geometric_descriptors.csv');
fprintf('Saved descriptor values: fast_tda_geometric_descriptors.csv\n');

%% 26. SAVE MATLAB DATA
save('fast_tda_complete_results.mat', ...
    'time','output','shock','shock_norm','Fmin','Fmax', ...
    'time_delay','time_delay_samp','win_pts','win_dur','step_size', ...
    'valid_idx','t_start','t_center','t_end','t_plot', ...
    'coeff','coeff_abs','coeff_norm','coeff_abs_norm', ...
    'coeff_abs_response','coeff_abs_response_norm', ...
    'desc','desc_norm','angle_unwrapped_deg');
fprintf('Saved MATLAB results: fast_tda_complete_results.mat\n');

%% 27. AUTOMATIC WINDOW SWEEP
automatic_sweep = true;
animation_pause = 0.02;

if automatic_sweep
    fprintf('\nStarting automatic TDA window sweep...\n');
    for k = 1:n_wins
        set(sld,'Value',k);
        cb_update(fig);
        drawnow;
        pause(animation_pause);
    end
    fprintf('Automatic sweep finished.\n');
end

fprintf('\nFast TDA analysis complete.\n');

%%
% LOCAL FUNCTION: UPDATE INTERACTIVE WINDOW
% =================================================================
function cb_update(fig)
ud = fig.UserData;
idx = round(get(ud.sld,'Value'));
idx = max(1,min(idx,numel(ud.valid_idx)));
i = ud.valid_idx(idx);

i1 = i;
i2 = i + ud.win_pts - 1 - ud.time_delay_samp;
i3 = i + ud.time_delay_samp;
i4 = i3 + ud.win_pts - 1 - ud.time_delay_samp;

if i2 > numel(ud.output) || i4 > numel(ud.output)
    return;
end

P = [ud.output(i1:i2),ud.output(i3:i4)];

x_l = ud.t_start(idx);
x_r = ud.t_end(idx);
set(ud.win_patch,'XData',[x_l x_r x_r x_l], ...
    'YData',[ud.min_acc ud.min_acc ud.max_acc ud.max_acc]);

set(ud.t_lbl,'String',sprintf( ...
    'start %.3f | center %.3f | end %.3f s', ...
    ud.t_start(idx),ud.t_center(idx),ud.t_end(idx)));

cla(ud.ax_ell);
plot_ellipse(P,ud.all_params(:,i),ud.ax_ell,ud.stdColors);
title(ud.ax_ell,sprintf( ...
    'Ellipse: start %.3f s, center %.3f s', ...
    ud.t_start(idx),ud.t_center(idx)));
xlabel(ud.ax_ell,'Amplitude');
ylabel(ud.ax_ell,'Delayed amplitude');
axis(ud.ax_ell,'equal');
grid(ud.ax_ell,'on');
box(ud.ax_ell,'on');

t_current = ud.t_plot(idx);
ud.cur_line.Value = t_current;
ud.cur_line_stk.Value = t_current;
drawnow limitrate;
end

%% 
% LOCAL FUNCTION: SAVE MAIN FIGURE AS SVG
% =================================================================
function cb_save_svg(fig)
ud = fig.UserData;
idx = round(get(ud.sld,'Value'));
idx = max(1,min(idx,numel(ud.valid_idx)));
fname = sprintf('fast_tda_window_%.3fs.svg',ud.t_center(idx));

try
    exportgraphics(fig,fname,'ContentType','vector','BackgroundColor','white');
catch
    saveas(fig,fname);
end

fprintf('Saved SVG: %s\n',fname);
end

%% 
% LOCAL FUNCTION: PLOT POINT CLOUD + IMPLICIT ELLIPSE
% =================================================================
function plot_ellipse(P,params,ax,cols)
plot(ax,P(:,1),P(:,2),'.','Color',cols(1,:),'MarkerSize',4);
hold(ax,'on');

if any(~isfinite(params))
    hold(ax,'off');
    return;
end

a = params(1); b = params(2); c = params(3);
d = params(4); e = params(5); f0 = params(6);

mins = min(P,[],1);
maxs = max(P,[],1);
xpad = 0.10*max(1e-9,maxs(1)-mins(1));
ypad = 0.10*max(1e-9,maxs(2)-mins(2));

xg = linspace(mins(1)-xpad,maxs(1)+xpad,350);
yg = linspace(mins(2)-ypad,maxs(2)+ypad,350);
[X,Y] = meshgrid(xg,yg);

Z = a*X.^2 + b*X.*Y + c*Y.^2 + d*X + e*Y + f0;

try
    contour(ax,X,Y,Z,[0 0],'Color',cols(2,:),'LineWidth',2);
catch
    warning('Ellipse contour failed for current window.');
end

axis(ax,'equal');
hold(ax,'off');
end

%% =================================================================
% LOCAL FUNCTION: HLS-GOLDEN ELLIPSE FIT
% =================================================================
function a = fit_ellipse_hls_golden(x,y)
if numel(x) < 5
    a = nan(6,1);
    return;
end

D_mat = [x.^2,x.*y,y.^2,x,y,ones(size(x))];
S_mat = D_mat'*D_mat;
Spp = S_mat(1:3,1:3);
Spq = S_mat(1:3,4:6);
Sqq = S_mat(4:6,4:6);

REG = 1e-6;
aq = Sqq(1,1)+REG;
bq = Sqq(1,2);
cq = Sqq(1,3);
eq = Sqq(2,2)+REG;
fq = Sqq(2,3);
gq = Sqq(3,3)+REG;

dv = aq*(eq*gq-fq^2) - bq*(bq*gq-fq*cq) + cq*(bq*fq-eq*cq);
if abs(dv) < 1e-12
    a = nan(6,1);
    return;
end

id = 1/dv;
Qinv = zeros(3);
Qinv(1,1) = (eq*gq-fq^2)*id;
Qinv(1,2) = -(bq*gq-cq*fq)*id;
Qinv(1,3) = (bq*fq-cq*eq)*id;
Qinv(2,1) = Qinv(1,2);
Qinv(2,2) = (aq*gq-cq^2)*id;
Qinv(2,3) = -(aq*fq-cq*bq)*id;
Qinv(3,1) = Qinv(1,3);
Qinv(3,2) = Qinv(2,3);
Qinv(3,3) = (aq*eq-bq^2)*id;

S_red = Spp-Spq*Qinv*Spq';

invD = [0 0.5 0; 0.5 0 0; 0 0 -1];
A_mat = invD*S_red;
[V,De] = eig(A_mat);
lambda = diag(De);

C = zeros(6,6);
C(1,3) = 2;
C(3,1) = 2;
C(2,2) = -1;

best_a = nan(6,1);
best_score = inf;

for k = 1:3
    ap = V(:,k);
    aq_v = -Qinv*(Spq'*ap);
    candidate = [ap;aq_v];

    if any(~isfinite(candidate))
        continue;
    end

    if candidate(2)^2 - 4*candidate(1)*candidate(3) >= 0
        continue;
    end

    if candidate'*C*candidate <= 0
        continue;
    end

    score = abs(lambda(k));
    if score < best_score
        best_score = score;
        best_a = candidate;
    end
end

if any(isnan(best_a))
    a = nan(6,1);
    return;
end

a = best_a;
scale = sqrt(abs(a'*C*a));
if scale > 1e-12
    a = a/scale;
end
end

%% =================================================================
% LOCAL FUNCTION: GEOMETRIC ELLIPSE DESCRIPTORS
% =================================================================
function desc = ellipse_descriptors(params)
desc = nan(12,1);
if any(~isfinite(params))
    return;
end

a = params(1); b = params(2); c = params(3);
d = params(4); e = params(5); f0 = params(6);

Q = [a b/2; b/2 c];
linear = [d;e];

if rcond(Q) < 1e-12
    return;
end

center = -0.5*(Q\linear);
xc = center(1);
yc = center(2);

Fc = a*xc^2 + b*xc*yc + c*yc^2 + d*xc + e*yc + f0;

Qsym = 0.5*(Q+Q');
[V,L] = eig(Qsym);
lambda = diag(L);

if any(abs(lambda) < 1e-12)
    return;
end

radius_squared = -Fc./lambda;
if any(~isfinite(radius_squared)) || any(radius_squared <= 0)
    return;
end

radii = sqrt(radius_squared);
[R_major,major_idx] = max(radii);
R_minor = min(radii);

if R_major <= 0 || R_minor <= 0
    return;
end

major_vector = V(:,major_idx);
angle_rad = atan2(major_vector(2),major_vector(1));
angle_rad = mod(angle_rad+pi/2,pi)-pi/2;
angle_deg = rad2deg(angle_rad);

axis_ratio = R_minor/R_major;
eccentricity = sqrt(max(0,1-axis_ratio^2));
ellipse_area = pi*R_major*R_minor;
focal_distance = sqrt(max(0,R_major^2-R_minor^2));

desc(1) = xc;
desc(2) = yc;
desc(3) = R_major;
desc(4) = R_minor;
desc(5) = 2*R_major;
desc(6) = 2*R_minor;
desc(7) = axis_ratio;
desc(8) = eccentricity;
desc(9) = angle_rad;
desc(10) = angle_deg;
desc(11) = ellipse_area;
desc(12) = focal_distance;
end

%% =================================================================
% LOCAL FUNCTION: SAFE MIN-MAX NORMALISATION
% =================================================================
function vn = safe_minmax(v)
vn = nan(size(v));
valid = isfinite(v);

if ~any(valid)
    return;
end

vmin = min(v(valid));
vmax = max(v(valid));
range_v = vmax-vmin;

if range_v < 1e-12
    vn(valid) = 0;
else
    vn(valid) = (v(valid)-vmin)/range_v;
end
end

%% =================================================================
% LOCAL FUNCTION: UNWRAP ELLIPSE ORIENTATION
% =================================================================
function angle_out = unwrap_orientation_pi(angle_in)
angle_out = nan(size(angle_in));
valid = isfinite(angle_in);

if ~any(valid)
    return;
end

idx = find(valid);
breaks = [0,find(diff(idx)>1),numel(idx)];

for k = 1:numel(breaks)-1
    group = idx(breaks(k)+1:breaks(k+1));
    angle_out(group) = 0.5*unwrap(2*angle_in(group));
end
end
