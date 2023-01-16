function times = detectThresholdCrossing_ov20220105(signal, direction, threshold)
% adapted from thresholdDetect by LLoyd Russell
% direction = 'above' or 'below' or 'flipflop'

% decide if looking for positive or negative threshold crossing
if strcmpi(direction, 'above')
    thresh_signal = signal > threshold;
elseif strcmpi(direction, 'below')
    thresh_signal = signal < threshold;
elseif strcmpi(direction, 'flipflop')
    thresh_signal = signal > threshold;
else
    warning('Unknown. Defaulting to signal > threshold')
    thresh_signal = signal > threshold;
end

% find times
% thresh_signal = signal > threshold;
% thresh_signal(thresh_signal(1:end-1) & thresh_signal(2:end)) = 0;
thresh_signal = diff(thresh_signal);
thresh_signal(2:end+1) = thresh_signal;
thresh_signal(1) = 0;
times = find(thresh_signal);

if strcmpi(direction, 'above') && (signal(1)<threshold)
    times = times(1:2:end);
elseif strcmpi(direction, 'above') && (signal(1)>threshold)
    times = times(2:2:end);
elseif strcmpi(direction, 'below') && (signal(1)<threshold)
    times = times(2:2:end);
elseif strcmpi(direction, 'below') && (signal(1)>threshold)
    times = times(1:2:end);
elseif strcmpi(direction, 'flipflop')
    times = times;
else
    warning('Unknown. Defaulting to signal > threshold')
    times = times(1:2:end);
end