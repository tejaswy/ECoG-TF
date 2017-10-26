function data = firfiltecog(data,fs, b)
%Uses method in  firfilt() - Pad data with DC constant, filter data with FIR filter,
%             and shift data by the filter's group delay
%
% Usage:
%   >> EEG = firfilt(EEG, b, nFrames);
%
% Inputs:
% -- required
%   EEG           - EEGLAB EEG structure
%   b             - vector of filter coefficients
%
% --Optional inputs:
%   nFrames       - number of frames to filter per block {default 4*fs}
%
% Outputs:
%   EEG           - EEGLAB EEG structure
%
% Note:
%   Higher values for nFrames increase speed and working memory
%   requirements.
%


if nargin < 2
    error('Not enough input arguments.');
end

nFrames=4*fs;
% Filter's group delay
if mod(length(b), 2) ~= 1
    error('Filter order is not even.');
end
groupDelay = (length(b) - 1) / 2;

dcArray = [1 size(data,2)+ 1];
%  keyboard

for iDc = 1:(length(dcArray) - 1)
    % Pad beginning of data with DC constant and get initial conditions
        ziDataDur = min(groupDelay, dcArray(iDc + 1) - dcArray(iDc));
        [temp, zi] = filter(b, 1, double([data(:, ones(1, groupDelay) * dcArray(iDc)) ...
                                  data(:, dcArray(iDc):(dcArray(iDc) + ziDataDur - 1))]), [], 2);

        blockArray = [(dcArray(iDc) + groupDelay):nFrames:(dcArray(iDc + 1) - 1) dcArray(iDc + 1)];
        for iBlock = 1:(length(blockArray) - 1)

            % Filter the data
            [data(:, (blockArray(iBlock) - groupDelay):(blockArray(iBlock + 1) - groupDelay - 1)), zi] = ...
                filter(b, 1, double(data(:, blockArray(iBlock):(blockArray(iBlock + 1) - 1))), zi, 2);

            % Update progress indicator
%             [step, strLength] = mywaitbar((blockArray(iBlock + 1) - groupDelay - 1), size(data, 2), step, nSteps, strLength);
        end

        % Pad end of data with DC constant
        temp = filter(b, 1, double(data(:, ones(1, groupDelay) * (dcArray(iDc + 1) - 1))), zi, 2);
        data(:, (dcArray(iDc + 1) - ziDataDur):(dcArray(iDc + 1) - 1)) = ...
            temp(:, (end - ziDataDur + 1):end);

        % Update progress indicator
%         [step, strLength] = mywaitbar((dcArray(iDc + 1) - 1), size(data, 2), step, nSteps, strLength);

    
    
end

end
