function [acw0mat, acw50mat] = acw_windowed_loop(ts, fs, windowsize)
%% Loop over channels
% ts is a matrix (channels x time)
% fs is sampling frequency in Hz
% windowsize is length of window in seconds

nchan = size(ts, 1);
acw0cell = cell(nchan, 1);
acw50cell = cell(nchan, 1);
for i = 1:nchan
    [acw0cell{i}, acw50cell{i}] = acw_windowed(ts(i, :), fs, windowsize);
end

acw0mat = cat(1, acw0cell{:});
acw50mat = cat(1, acw50cell{:});
end