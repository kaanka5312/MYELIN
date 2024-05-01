function [acw0s, acw50s] = acw_windowed(ts, fs, windowsize)
%% Calculate ACWs in sliding window fashion with no overlap
% ts is 1D vector
% Windowsize is in unit of seconds
% fs is sampling frequency in Hz

window_samples = windowsize * fs;
datalength = length(ts);

acw0s = [];
acw50s = [];
i = 0;
while true
    indx = ((i*window_samples)+1):((i+1)*window_samples);
    if indx(end) > datalength
        break
    end

    [acw0, acw50] = acw(ts(indx), fs);
    acw0s = [acw0s, acw0];
    acw50s = [acw50s, acw50];
    i = i + 1;
end
end