function Trials = getEpochs(data, Labels,fs, tau,delay,numbins, groupDelay)
% extracts epochs from data and labels and puts them in Trials structure
%
% Inputs:

% Required:
% data : ecog data or usually band power time series of size ch x time
% fs   : sampling frequency
% Labels : Usually a matrix with #Epochs rows  x three columns :
% [startIndex, stopIdex, Label] , indices are always ecog time sample
% indices, not absolute time. if you have absolute time , use getIndices to
% get closest ecog sample index.


% Optional:
% We extract features from a specific window for each trial

% tau = window size centered around onset in ms   Default: trial width
% delay = window center from trial onset in ms : Default : mid point of trial,
%   negative delay means window center is after trial onset
%   zero delay is window center is at trial onset
% numbins = number of bins to break down the tau ms window: default=1
% groupDelay = if data is BP or HGP, the filters induce group delay and should be accounted for when computing features;
%               Default=0;

% Output:
% Trials: Structure of size number of Epochs with fields: Features (numchannels * numbins vector), label

if nargin<3    error('Not enough info')   ;  end
if nargin<4  tau= []; end
if nargin<5 delay= []; end
if nargin<6 numbins=1; end
if nargin<7 groupDelay=0; end

if ~isempty(tau)  tauPts=tau*fs/1000;  end
if ~isempty(delay) delayPts=delay*fs/1000; end

tstart= Labels(:,1);
tstop= Labels(:,2);
numEpochs = length(Labels);
numChannels = size(data,1);
for i=1:numEpochs
    Trials{i}.label = Labels(i,3);
    
    if isempty(tau) tauPts = tstop(i)-tstart(i)+1; end
    if isempty(delay)  delayPts = ceil((tstop(i)+tstart(i))/2) ;  end
    win      = tstart(i)-delayPts-floor(tauPts/2)+groupDelay:tstart(i)-delayPts+floor(tauPts/2)-1+groupDelay;
    subwin   = reshape(win,[],numbins);
    temp = zeros(numbins,numChannels);
    for k=1:numbins
        temp(k,1)=   mean(data(:,subwin(:,k)));
        
    end
    Trials{i}.features = reshape(temp,[],1);
end


end