function detrendedData = detrendData(data, fs, movingwin)
% Uses chronux's locdetrend: Remove running line fit (using local linear regression)-continuous
%  processes
%
% Inputs:
% -- required
%  data         (data as a matrix ch x times or a single vector)
%  fs           (sampling frequency)
%-- optional
%  movingwin    (length of moving window, and stepsize in seconds) [window winstep] - optional.
%                   Default. window=full length of data (global detrend).
%                   winstep=window -- global detrend
%
% Output:
% data:         (locally detrended data)


if nargin<2
    error('data and fs are required inputs');
end



fprintf('\ndetrend the zscored  signals\n');

if nargin<3 
    detrendedData   = (locdetrend(data',fs))';
else
    detrendedData   = (locdetrend(data',fs,movingwin))';
end
