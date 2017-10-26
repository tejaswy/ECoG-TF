function  ecogInd=TimeAlign(eventTime,ecogTimeStamp)
% gives ecog sample index closest to the event Times
% inputs: 
% required:
% eventTime : a vector of event times 
% ecogTimeStamp:  time stamps of ecog samples 
%NOTE: bothe eventTime and ecogTimeStamp should be in the same units
% Output:  ecog indices that match with eventTime

if nargin<2
    error('Insufficient number of inputs')
end

numEvents            = length(eventTime);
ecogInd              = zeros(size(eventTime));

% this script finds the matching index of eventTime in ecog timebase
for i=1:numEvents
    [~,ecogInd(i)]    = min(abs(eventTime(i) - ecogTimeStamp));   
end

end