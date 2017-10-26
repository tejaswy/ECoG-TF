function [features, classLabels] = restructureData(ecog, fs, Labels, params)

% ecog is numChannels x time
% Features would be 28x6xnumChannelsxnumTrials

[numChannels,numSamples]  = size(ecog);

freqRange  = params.freqRange;
numBands   = params.numBands;
numBins    = params.numBins;
binWidth   = params.binWidth; % in milliseconds



% tau is the window width of each trial
tau =    (numBins*binWidth)*fs/1000;
delay=150;
data =    zeros(numChannels,numSamples,numBands);
pad=1;    method='squared';
for i=1: numBands
    opts.freqrange= freqRange(i:i+1);
    [data(:,:,i) , ~] =  estimatePowerSeries(ecog, fs, opts,pad, method );
end

% mean power in bin is mean of amp^2


if params.toNorm ==1  % normalize by scoring each channel's frequency component
    meanPower = mean(data,2); stdPower = std(data,[],2);
    meanPower = repmat(meanPower,1,size(data,2),1);
    stdPower = repmat(stdPower,1,size(data,2),1);
    data = (data-meanPower)./stdPower;
    
end



Labels= Labels(2:end-1,[1 2 3]);
% leave out the first and last trial to avoid negative indices when centering with tau

if    params.toRand ==1   % to randomize originial trials
    randOrd   =   randperm(size(Labels,1));
    Labels    =   Labels(randOrd,:);
end

numEpochs= length(Labels);


% Split the original trials into tau ms trials


trialInfo=[];
if params.leaveTransition==0
    % Dont leave
    
    for i=1:numEpochs
        win      = Labels(i,1):Labels(i,2);
        tempLen =length(win) - mod(length(win),tau);
        subwin   = reshape(win(1:tempLen),tau,[]);
        tStart= subwin(1,:)';
        tStop=  subwin(end,:)';
        tLabels = repmat(Labels(i,3),length(tStart),1);
        trialInfo = [trialInfo; tStart,tStop,tLabels];
    end
    
elseif params.leaveTransition==1
    %leaving out first and last transition trials
    
    for i=1:numEpochs
        win      = Labels(i,1):Labels(i,2);
        tempLen =length(win) - mod(length(win),tau);
        subwin   = reshape(win(1:tempLen),tau,[]);
        tStart= subwin(1,2:end-1)';
        tStop=  subwin(end,2:end-1)';
        tLabels = repmat(Labels(i,3),length(tStart),1);
        trialInfo = [trialInfo; tStart,tStop,tLabels];
    end
    
else
    
    for i=1:numEpochs  %Use only first trials
        tStart= Labels(i,1)-tau/2 ;
        tStop= Labels(i,1)+tau/2;
        tLabels =  Labels(i,3);
        trialInfo = [trialInfo; tStart,tStop,tLabels];
    end
    
end

% from every starting point , take a tau ms window
numTrials = length(trialInfo);


features =   zeros(numChannels,numBands,numBins,numTrials );
classLabels = zeros(1,numTrials);

for i=1:numTrials
    %     midPoint= floor((trialInfo(i,1)+trialInfo(i,2))/2);
    %     win  =    midPoint-floor(delay)+1: midPoint+ceil(delay);
    win = trialInfo(i,1)-delay : trialInfo(i,2)-delay-1;
   
    tempLen =length(win) - mod(length(win),tau);
    temp =    data(:,win(1:tempLen),:);
    temp  =   reshape(temp,numChannels,[],numBins,numBands); % chxbinwidthxnumBinsxnumBands
    temp2(1:numChannels,:,:)  =   squeeze(mean(temp,2)); % chx numBins xnumBands
    
    features(:,:,:,i)= permute(temp2,[1,3,2]);
    classLabels(1,i) =    trialInfo(i,3);
end

% to remove line noise bands
if params.noLine
    freqCols = [1:13,17:28,32:34];
    features =   features(:,freqCols,:,:);
end
end

