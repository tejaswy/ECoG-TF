function  [measure,indelec]=   computeMeasure(data,fs,elec,threshold,measure,norm,precomp,freqrange)
%computes the measure required to identify bad electrodes


numChannels = size(data,1);


if nargin <2                   error('data and fs are required'); end
if nargin<3 || isempty(elec)        elec=1:numChannels; end
if nargin<4 || isempty(threshold)   threshold= 0; end
if nargin<5 || isempty(measure)     measure='kurt'; end
if nargin<6 || isempty(norm)        norm='off'; end
if nargin<7|| isempty(precomp)      precomp=[]; end
if nargin<8 || isempty(freqrange)   freqrange=[1, fs/2]; end


opt.norm    = norm;
opt.measure = measure;
opt.precomp = precomp;
opt.freqrange = freqrange;
opt.elec      = elec;
opt.threshold = threshold;


indelec = [];
measure = [];


% compute the joint probability
% -----------------------------
if strcmpi(opt.norm, 'on')
    normval = 2;
else
    normval = 0;
end;
if strcmpi(opt.measure, 'prob')
    fprintf('Computing probability for channels...\n');
    [ measure, indelec ] = jointprob( reshape(data(opt.elec,:,:), length(opt.elec), size(data,2)*size(data,3)), opt.threshold, opt.precomp, normval);
elseif strcmpi(opt.measure, 'kurt')
    fprintf('Computing kurtosis for channels...\n');
    [ measure, indelec ] = rejkurt( reshape(data(opt.elec,:,:), length(opt.elec), size(data,2)*size(data,3)), opt.threshold, opt.precomp, normval);
else
%     fprintf('Computing spectrum for channels...\n');
    error('I dunno how this works. Will get back to it later')
    %     [measure freq] = pop_spectopo(EEG, 1, [], 'EEG' , 'plot','off');
end
% colors = cell(1,length(opt.elec)); colors(:) = { 'k' };
% colors(find(indelec)) = { 'r' }; colors = colors(end:-1:1);
fprintf('%d electrodes labeled for rejection\n', length(find(indelec)));

% output variables
indelec     = find(indelec)';

if isempty(indelec), return; end;

end

