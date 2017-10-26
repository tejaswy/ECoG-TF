function featRanks =  getFeatRanks(Trials, method)
% Rank features based on some criterion
% method ='ANOVA', 'LDA','r2'
% inputs:
% Trials is a structure of size number of Epochs,
% with fields :  features (numChannels* bins)
%             :  label
% method: 1. 'anova':  p value from anova
%           2.  'lda': LDA weights
%
% output: rank order of fetaures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

features       =     [PitchTrials.features]';
% features should be numEpochs x numDimensions
Labels         =     [Trials.label];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(method,'anova')
    featRanks = getAnovaRanks(features,Labels);
elseif strcmp(method,'lda')
    featRanks = getLdaRanks(features,Labels);
else
    error('invalid method')
end




end

function featRanks = getLdaRanks(features,Labels)
cls   = fitcdiscr(features,Labels(:),'DiscrimType','Linear','Gamma',1);
W     = cls.Coeffs(1,2).Linear;
[~ , featRanks ]= sort(abs(W));
end



function featRanks = getAnovaRanks(features,labels)

for i = 1:numchannels
    X= features(i,:);
    p(i)= anova1(X,labels(:)','off');
end
[sorted_p,featRanks]=sort(p);

end
