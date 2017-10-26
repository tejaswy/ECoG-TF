function  filtSignal =  movAvgFilter(data, fs, windowWidth,method)
% smooth data using smooth of tsmovavg
% 
%  Inputs:
%--required
%  data         (data as a matrix ch x times or a single vector)
%  fs           (sampling frequency)
% windowWidth :    window width in s
% 
%--Optional:
% method : smooth or  tsmovavg % for differences look at doc
% smooth(centered window),tsmovavg (past points), Default: 'tsmovavg',
% remember smooth doesnt account for end points
% 
% Output:
%filtSignal:


if nargin<3
    error('\n Insufficient inputs')
end
if nargin<4|| isempty(method)
    method='tsmovavg';
end

numChannels= size(data,1);
windowWidth =  windowWidth  * fs ; % convert ms to samples

if strcmp(method,'smooth')
    filtSignal= smooth(data',windowWidth,'moving');
    filtSignal = reshape(filtSignal,size(data));
elseif strcmp(method,'tsmovavg')
    filtSignal= tsmovavg(data,'s',windowWidth);
else
    error('\n Invalid method')
end


end