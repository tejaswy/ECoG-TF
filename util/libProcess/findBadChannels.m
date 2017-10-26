function [badChannels,measureVal] = findBadChannels(data,fs,elec,threshold,measure,norm,precomp,freqrange)
%%% find bad channels using  pop_rejchan() -
% reject artifacts channels in an EEG dataset using joint  probability of the recorded electrode.
% I use pop_rejchan() to first find measure and then determine thresh as
% mean+3*stddev of measures of all channels
%  Inputs:
%--required
%  data         (data as a matrix ch x times or a single vector)
%  fs           (sampling frequency)
%-- Optional inputs:
%   'elec'     - [n1 n2 ...] electrode number(s) to take into
%                consideration for rejection
%   'threshold' - [max] absolute thresold or activity probability
%                 limit(s) (in std. dev.) if norm is 'on'.
%   'measure'  - ['prob'|'kurt'] compute probability 'prob', kurtosis 'kurt'
%                or spectrum 'spec' for each channel. Default is 'kurt'.
%   'norm'     - ['on'|'off'] normalize measure above (using trimmed
%                normalization as described in the function jointprob()
%                and rejkurt(). Default is 'off'.
%   'precomp'  - [float array] use this array instead of computing the 'prob'
%                or 'kurt' measures.
%   'freqrange' - [min max] frequency range for spectrum computation.
%                Default is 1 to sampling rate divided by 2. The average
%                of the log spectral power is computed over the frequency
%                range of interest.
%
% Outputs:
%
%  badChannels  - indices of rejected electrodes
%   measure   - measure value for each electrode
%

%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if nargin<2, error('data and fs are required'); end
if nargin<3 || isempty(elec), elec= 1:size(data,1); end
if nargin<4,  threshold=[]; end
if nargin<5,  measure=[]; end
if nargin<6, norm=[]; end
if nargin<7, precomp=[]; end
if nargin<8, freqrange=[]; end

if isempty(norm), norm='off'; end 

[measureVal,~]=   computeMeasure(data,fs,elec,threshold,measure,norm,precomp,freqrange);

meanSpec = mean(measureVal);
stdSpec  = std( measureVal);

if strcmp(norm,'off'),
    thresh = meanSpec+(stdSpec*3);
else
    thresh =  stdSpec*3;
end
% if norm is on, then the measureVal are zscored, so use 3*std for threshold
ind = find(abs(measureVal)> thresh);
badChannels = elec(ind);



end