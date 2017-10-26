function plotHist(Trials,Dims)
% plots histograms of feature dist
% inputs:
% Trials: struct with features and Labels
% Dims:  Dimensions of the feature vector , default=1
% output
% f= figure handle

if nargin<2
    Dims=1;
end

features = [Trials.features]';
Labels   = [Trials.labels];

numTrials = length(Labels);

numClasses = length(unique(Labels));

colorStr = 'rgbkm';

for i=Dims
    for j=1:numClasses
    hist(features(i,Labels==j));
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',colorStr(j),'EdgeColor','k','facealpha',0.5)
    hold on;
    end
end



end
