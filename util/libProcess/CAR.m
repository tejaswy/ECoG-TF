function carData = CAR(data, channel_groups, bad_channels)
% Common Average Referencing
% Input:
%       -- required
%       data (in form  channels x samples)
%       -- optional
%           channelGroups: structure with channels to be referenced
%           together in one group.
%           For Example: If there are 64 channels and you want channels 1
%           to 8 referenced together and the rest referenced together, then
%           channelGroups will be {[1:8],[9:64]}
%           Default: All channels considered for common reference.
% 
%           bad_channels:  Array with channels to be excluded from CAR even if they
%           are included in channel_groups. Default =[]
%
% Output:
%      -- carData = Common Average Referenced data  (in form  channels x samples)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin<1
    error('data is a required input.')
end

numChannels= size(data,1);

if nargin<2 || isempty(channel_groups)
    channel_groups = {[1:numChannels]};
end

if nargin<3
    bad_channels  =[];
end


carData      = zeros(size(data));

% for each channel, get the channel's group.
numGroups   =    length(channel_groups);
grId    =   zeros(numChannels,1);
for i= 1: numGroups
    grId(channel_groups{i}) = i;
end
% channels which are not in any group all under group zero

for ichan=1:numChannels
%     keyboard
    cInd =find(grId== grId(ichan));   % all the channels belonging to same group as channel: ichan
    bchan = intersect(cInd,bad_channels); % see if there is any bad channel in that group
    % if there is a bad channel in that group do not use it in CAR
    % remember not to use bad channels in any further analysis.
    goodchan = setdiff(cInd,bchan);
    if isempty(goodchan)  % if all channels in the group are bad, dont do anything
        carData(ichan,:) = data(ichan,:) ;
    else
        groupAvg = mean(data(goodchan,:));
        carData(ichan,:) = data(ichan,:) - groupAvg;
    end
end

end

