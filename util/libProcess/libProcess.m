function libProcess
% This library has functions required for data pre-processing of ecog data. Uses chronux
% and eeglab packages.
% Certain variables are used repeatedly in this library.
%
% data:  always a matrix with dimensions numChannels x Samples
% Fs:     sampling frequency
%
% freqrange: band pass frequencies in the form [fmin fmax] where fmin >=0 and
%        fmax<=Fs/2. optional (default [0 Fs/2])
%
% Other parameters are discussed in individual routines as and when they
% are used.
%
% Here is the list of functions:
%
% carData = CAR(data, channel_groups, bad_channels)
% [cleanData, Sorig, Sclean, f, amps, freqs, g] =  LineNoiseFilter(ecog,fs,opts)
% zscoredData = zscoreActivity(data,baseline)
% detrendedData = detrendData(data, fs, movingwin)
% [badChannels,measure] = findBadChannels(data,fs,elec,threshold,measure,norm,precomp,freqrange)
% [b,filtorder] = filterDesign(fs,locutoff,hicutoff,filtorder,revfilt,plotfreqz,minphase,wtype)
% [estPowers, gd] =  estimatePowerSeries(data, fs, opts,pad, method )
% filtSignal =  movAvgFilter(data, fs, windowWidth,method)
% filtData = firfiltecog(data, fs,b)
%
% Author : Tejaswy Pailla,
% Date   : 5/28/2016



%%%%%%%%%%%%%% test script testing all the above functions %%%
%%%% test data is finger flexion training data from bp
clc
clear all
close all
load('tempEcog.mat');
ecog = tempEcog;
fs   = 1000;

%% FOR ALL SCRIPTS
%rows of data are channels and columns are time samples, if not, do transpose
if size(ecog,1)>size(ecog,2)
    ecog=ecog';
end

%% Line noise removal
% opts can  have fields 'linefreqs', 'scanforlines','p'or'alpha','bandwidth','winsize',
% 'winstep,'tau','pad','computepower','normSpectrum','verb','plotfigures'


opts.p=0.005; opts.bandwidth= 10; opts.tau=50;
[cleanData, Sorig, Sclean, f, amps, freqs, g] =  LineNoiseFilter(ecog,fs,opts);

%Uncomment to  See how the spectra look like
% figure
% for ch =1 % plot for one channel
% plot(f, Sorig(ch,:)); hold on; plot(f, Sclean(ch,:),'r');
% legend('Orig','Clean');
% title('before and after line noise removal')
% xlabel('freq'); ylabel('10*log10(pxx)')
% pause;  hold off;
% end

% note that peaks at 180 and 300 Hz are hard to remove
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Remove bad channels using cleanData generated after removing line noise

data = cleanData;
elec=[]; threshold=[]; norm=[]; precomp=[];
measure='prob'; freqrange=[2 200];
[badChannelsProb,measure] = findBadChannels(data,fs,elec,threshold,measure,norm,precomp,freqrange);
fprintf('Bad channels with prob measure:\n');
fprintf(' %d \n', badChannelsProb)
measure='kurt';
[badChannelsKurt,measure] = findBadChannels(data,fs,elec,threshold,measure,norm,precomp,freqrange);
fprintf('Bad channels with kurt measure:\n');
fprintf(' %d \n', badChannelsKurt)

fprintf('common channels from both methods\n');
fprintf('%d\n',intersect(badChannelsKurt, badChannelsProb))

% I would go with prob measure , one channel might be super responsive to
% task and be picked out as bad one by kurt.. just a theory though
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Common Average Referencing
%%%%%%%%%%%%%%%%%%%%%%% Now do CAR on cleanData with bad channels left out
bad_channels = badChannelsProb;
channel_groups=[];
carData = CAR(data, channel_groups, bad_channels);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%zscore with baseline data if available%%%
% baseline is when the subject is at rest , you can do this even with ISI
% epochs
baselineIndices=[];
zscoredData = zscoreActivity(carData,baselineIndices);

% % uncomment to see the plot
% ch=1;
% plot(carData(ch,:)); hold on; plot(zscoredData(ch,:),'r');
% legend('orig','zscored')
% %%

%% detrending

movingwin=[5 1]; % in seconds
detrendedData = detrendData(carData, fs, movingwin);

%  % uncomment to see the plot
% ch=1;
% smCar =  movAvgFilter(carData(ch,:), fs, 5,[]);
% smDet =  movAvgFilter(detrendedData(ch,:), fs, 5,[]);
%
% plot(smCar(ch,:)); hold on; plot(smDet(ch,:),'r');
% legend('orig','detrended')

%% TESTING FILTER DESIGN

revfilt=0; plotfreqz=0;minphase=0; wtype='hamming'; filtorder=[];
locutoff = 15; hicutoff=40; 
[b,filtOrder] = filterDesign(fs,locutoff,hicutoff,filtorder,revfilt,plotfreqz,minphase,wtype);
freqz(b,1,[],fs); 
suptitle('Filter amp and phase response')
[gd,w]=grpdelay(b,1);
gd= mode(gd); % this is the group delay introduced by the filter;


ch=1;
temp= carData(ch,:);
filtTemp = filter(b,1,temp);%firfiltecog(temp,b);
[pxxOrig,f] = pwelch(temp,fs);
[pxxFilt,f] = pwelch(filtTemp,fs);
f= linspace(0,fs/2,length(f));
plot(f,10*log10(pxxOrig));
hold on
plot(f,10*log10(pxxFilt),'r'); 
legend('orig','Filtered')
 xlabel('freq'); ylabel('10*log10(pxx)');
 title('before and after filtering');



%% testing filtering function firfiltecog that accounts for group delay

revfilt=0; plotfreqz=0;minphase=0; wtype='hamming'; filtorder=[];
locutoff = 15; hicutoff=40; 
[b,filtOrder] = filterDesign(fs,locutoff,hicutoff,filtorder,revfilt,plotfreqz,minphase,wtype);

filtTemp = firfiltecog(temp,fs,b);

 
%% Get powers from data by filtering and  squaring amp 

temp= carData(1,:);

opts.freqrange= [15 40]; pad=1;
method='squared'; % only analytic amp
[estPowers, gd] =  estimatePowerSeries(temp, fs, opts,pad, method );
smLow= movAvgFilter(estPowers, fs, 2,'tsmovavg');

opts.freqrange= [75 150];pad=1;
method='squared'; % only analytic amp
[estPowers, gd] =  estimatePowerSeries(temp, fs, opts,pad, method );
smHigh= movAvgFilter(estPowers, fs, 2,'tsmovavg');

figure

plot(smLow','r');hold on
plot(smHigh,'k'); 
legend('low power','high power')


%% testing smoothing function

% filtSignal =  movAvgFilter(temp, fs, windowWidth,method)










