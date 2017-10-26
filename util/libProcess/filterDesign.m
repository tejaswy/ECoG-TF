function [b,filtorder] = filterDesign(fs,locutoff,hicutoff,filtorder,revfilt,plotfreqz,minphase,wtype)
% filter design script using modified pop_eegfiltnew(),pop_firws
% Designs windowed sinc type I linear phase FIR filter
%%%% NOTE: minphase is always 0
% Inputs:
% --required:
% fs: Sampling Frequency
% %   locutoff  - lower edge of the frequency pass band (Hz)
%               {[]/0 -> lowpass}
% -- optional
%   hicutoff  - higher edge of the frequency pass band (Hz)
%               {[]/0 -> highpass}

%   filtorder - filter order (filter length - 1). Mandatory even
%   revfilt   - [0|1] invert filter (from bandpass to notch filter)
%               {default 0 (bandpass)}
%   plotfreqz - [0|1] plot filter's frequency and phase response
%               {default 0}
%   minphase  - scalar boolean minimum-phase converted causal filter
%               {default false}
%   'wtype'       - char array window type. 'rectangular', 'bartlett',
%                   'hann', 'hamming', 'blackman', or 'kaiser' {default
%                   'blackman'}
% Outputs:
%  b  - filter coefficients to use with firfilt(b,1,x)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%   Window based filters' transition band width is defined by filter
%   order and window type/parameters. Stopband attenuation equals
%   passband ripple and is defined by the window type/parameters. Refer
%   to table below for typical parameters. (Windowed sinc) symmetric FIR
%   filters have linear phase and can be made zero phase (non-causal) by
%   shifting the data by the filters group delay (what firfilt does by
%   default). Pi phase jumps noticable in the phase reponse reflect a
%   negative frequency response and only occur in the stopband. pop_firws
%   also allows causal filtering with minimum-phase (non-linear!) converted
%   filter coefficients with similar properties. Non-linear causal
%   filtering is NOT recommended for most use cases.
%
%               Beta    Max stopband    Max passband    Max passband    Transition width    Mainlobe width
%                       attenuation     deviation       ripple (dB)     (normalized freq)   (normalized rad freq)
%                       (dB)
%   Rectangular         -21             0.0891          1.552           0.9 / m*             4 * pi / m
%   Bartlett            -25             0.0562          0.977                                8 * pi / m
%   Hann                -44             0.0063          0.109           3.1 / m              8 * pi / m
%   Hamming             -53             0.0022          0.038           3.3 / m              8 * pi / m
%   Blackman            -74             0.0002          0.003           5.5 / m             12 * pi / m
%   Kaiser      5.653   -60             0.001           0.017           3.6 / m
%   Kaiser      7.857   -80             0.0001          0.002           5.0 / m
%   * m = filter order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin<2
    error('fs and freq must be provided')
end

if nargin < 3 || isempty( hicutoff)
    hicutoff = [];
end
if nargin < 4  || isempty( filtorder)
    filtorder = [];
end
if nargin < 5 || isempty(revfilt)
    revfilt = 0;
end

if nargin < 6 || isempty(plotfreqz)
    plotfreqz = 0;
end
if nargin < 7 || isempty(minphase)
    minphase = 0;
end
if nargin<8
    wtype ='hamming';
end

% Constants
TRANSWIDTHRATIO = 0.25;
fNyquist = fs / 2;


% Check arguments
if locutoff == 0, locutoff = []; end
if hicutoff == 0, hicutoff = []; end
if isempty(hicutoff) % Convert highpass to inverted lowpass
    hicutoff = locutoff;
    locutoff = [];
    revfilt = ~revfilt;
end

edgeArray = sort([locutoff hicutoff]);

if isempty(edgeArray)
    error('Not enough input arguments.');
end
if any(edgeArray < 0 | edgeArray >= fNyquist)
    error('Cutoff frequency out of range');
end

if ~isempty(filtorder) && (filtorder < 2 || mod(filtorder, 2) ~= 0)
    error('Filter order must be a real, even, positive integer.')
end

% Max stop-band width
maxTBWArray = edgeArray; % Band-/highpass
if revfilt == 0 % Band-/lowpass
    maxTBWArray(end) = fNyquist - edgeArray(end);
elseif length(edgeArray) == 2 % Bandstop
    maxTBWArray = diff(edgeArray) / 2;
end
maxDf = min(maxTBWArray);



% Transition band width and filter order
if isempty(filtorder)
    
    % Default filter order heuristic
    if revfilt == 1 % Highpass and bandstop
        df = min([max([maxDf * TRANSWIDTHRATIO 2]) maxDf]);
    else % Lowpass and bandpass
        df = min([max([edgeArray(1) * TRANSWIDTHRATIO 2]) maxDf]);
    end
    filtorder = pop_firwsord(wtype, fs, df);
    % pop_firwsord Inputs:
%   wtype - char array window type. 'rectangular', 'bartlett', 'hann',
%           'hamming', {'blackman'}, or 'kaiser'
%   fs    - scalar sampling frequency {default 2}
%   df    - scalar requested transition band width
%   dev   - scalar maximum passband deviation/ripple (Kaiser window
%           only)

%     filtorder = 3.3 / (df / fs); % Hamming window
%     filtorder = ceil(filtorder / 2) * 2; % Filter order must be even.
%     
else
    
  filtorder  = ceil(m / 2) * 2; % Make filter order even (type 1)
    
end

filterTypeArray = {'lowpass', 'bandpass'; 'highpass', 'bandstop (notch)'};

% Passband edge to cutoff (transition band center; -6 dB)
dfArray = {df, [-df, df]; -df, [df, -df]};
cutoffArray = edgeArray + dfArray{revfilt + 1, length(edgeArray)} / 2;

% Window
winArray = windows(wtype, filtorder + 1);

% Filter coefficients
if revfilt == 1
    filterTypeArray = {'high', 'stop'};
    b = firws(filtorder, cutoffArray / fNyquist, filterTypeArray{length(cutoffArray)}, winArray);
else
    b = firws(filtorder, cutoffArray / fNyquist, winArray);
end

 causal = minphase;
% Plot frequency response
if plotfreqz
     plotfresp(b, 1, [], fs, causal);
end



end

