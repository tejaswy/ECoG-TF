function zscoredData = zscoreActivity(data,baseInd)
% zscore relative to the mean and standard deviation of baseline data for each channel
% inputs:
% --required
% data is ecog in chan x time
% baseline is ecog indices of baseline activity
% outputs:
% zscoredData

if nargin<2  
    error('data and baseline are required inputs');
end
if isempty(baseInd)
    fprintf('\nNo baseline activity. Returning input data\n')
    zscoredData = data;
else

% fprintf('\nZscoring with Baseline data...\n')
[~,time]   =  size(data);
baseData     = data(:,baseInd);
baseMean     = mean(baseData,2);
basestd      = std(baseData,1,2);
zscoredData  = (data -repmat(baseMean,1,time))./repmat(basestd,1,time);

end

