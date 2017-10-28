function waveform = get_waveform(waveform, sps)
%GET_WAVEFORM Summary of this function goes here
%   Detailed explanation goes here

if ischar(waveform)
    switch waveform
        case 'gcamp6s'
            rise_peak = 0.4; % seconds
            decay_half = 1.0; % seconds
            
        case 'gcamp6m'
            rise_peak = 0.24; % seconds
            decay_half = 0.79; % seconds
            
        case 'gcamp6f'
            rise_peak = 0.14; % seconds
            decay_half = 0.32; % seconds
            
        otherwise
            error('Invalid waveform name: %s.', waveform);
    end
    
    % decay to
    decay_to = -5;
    
    % convert to samples
    rise_peak_smp = round(rise_peak * sps);
    decay_smp = round(decay_half * sps * decay_to / log(0.5)); % decay_to / log(0.5) converts from time to half until time to full decay
    
    % make waveform: rise
    if rise_peak_smp > 1
        rise = 1 - exp(linspace(0, -3, rise_peak_smp));
    else
        rise = [];
    end
    
    % make waveform: decay
    decay = exp(linspace(0, decay_to, decay_smp));
    
    % drop the first rise time
    waveform = [rise(2:end) decay];
elseif ~isvector(waveform)
    error('Waveform must be a vector or string.');
end

end
