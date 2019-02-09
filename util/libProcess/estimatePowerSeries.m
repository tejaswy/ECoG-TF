function  [estPowers, group_delay] =  estimatePowerSeries(data, fs, opts,pad, method )
% Estimate power series by band pass filtering and squaring or using
% hilbert method,
%
% Inputs:
% --required
%  data         (data as a matrix ch x times or a single vector)
%  fs           (sampling frequency)
% -- optional
% opts          struct with details abt filter
%                @Params: freq range: default [75,150] or
%                           filtCoeff  : filter design coeffcients
% pad = 0 , dont account for group delay (use filter)
%       =1 ; use firfiltecog , account for group delay
%
% Output:
% estimated power time series and group delay (gd) introduced by filters
if nargin<2
    error('\n Insufficient inputs')
end
if nargin<3
    opts.freqrange=[75,150];
end

if nargin<4
    method = 'squared';
end

if isfield(opts, 'filtCoeff')
    filtCoeff = opts.filtCoeff;
elseif isfield(opts, 'freqrange')
    % do filter design
    filtCoeff = filterDesign(fs,opts.freqrange(1),opts.freqrange(2));
    % syntax:     filtcoeff = filterDesign(fs,locutoff,hicutoff,filtorder,revfilt,plotfreqz,minphase,wtype)
    
end



if pad==1
    filtSignal = firfiltecog(data,fs, filtCoeff);    
    group_delay=0;
else
    filtSignal = filter(filtCoeff,1,data')';
    
    group_delay=mode(grpdelay(filtCoeff));
end


if strcmp(method,'hilbert')
    estPowers   = (abs(hilbert(filtSignal')))';
elseif strcmp(method,'squared')
    estPowers = (filtSignal).^2;
else
    error('Invalid method');
end


end
