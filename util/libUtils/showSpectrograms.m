function f=showSpectrograms(data,fs,timeWin,params,toZscore)
% show   spectrograms from ecog data
% inputs:
% data: time series data from one channel
% fs: sampling freq
% timeWin: window for which spectrogram to be shown: default(entire time)
%params: struct for computing spectrogram with fields
%    movingwin=[winwidth, winStep] in sec
%    Fs=sampling frequency
%    fpass=  % frequencies of interest
%    tapers
%toZscore=1 to show zscored S , 0 otherwise(Default:1)



if nargin<2
    error('data and fs are required')
end
if nargin<3
    timeWin=[0 2]; % show 2 second spectrogram by default
end
if nargin<4
    params=[];
end
if nargin<5
    toZscore=1;
end

params =getdefaultParams(fs,params);

ecogwin = timeWin(1)*fs+1:  timeWin(2)*fs;

[S,t,f] = mtspecgramc((data(1,ecogwin))',params.movingwin,params);

f=plot_matrix(zscore(S),t,f,'n'); grid on
colorbar('off')
colorbar('SouthOutside')
xlabel('Time [s]')
ylabel('Freq [Hz]')
title('Spectrogram')


end


function params=getdefaultParams(fs,params)
% parameters for spectrogram
if ~isfield(opts,'movingwin')
    params.movingwin         = [0.5 0.05]; % moving window size in [s]
end
if ~isfield(opts,'Fs')
    params.Fs         = fs;
end
if ~isfield(opts,'fpass')
    params.fpass      = [0 150];
end
if ~isfield(opts,'tapers')
    params.tapers     =  [4 7];
end

end