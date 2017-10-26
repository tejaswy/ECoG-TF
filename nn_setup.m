function [net, info,imdb] = nn_setup(AEWeights,inputs, fflabels,expDir,nH_l1,nH_l2,numClasses,fixNumLayers)
if nargin<3
    expDir='tempDir';
    fpath = fullfile(pwd,expDir);
    if isequal(exist(fpath, 'dir'),7) % 7 = directory
        rmdir('tempDir','s');
    end
end
input_h= size(inputs,1);
input_w= size(inputs,2);
input_d= size(inputs,3);

opts.batchNormalization = true;
opts.networkType = 'simplenn' ;

sfx = opts.networkType ;
if opts.batchNormalization, sfx = [sfx '-bnorm'] ; end

opts.train = struct() ;


imdb  = getFFImdb(inputs, fflabels,input_h,input_w,input_d,numClasses);

net= initializeNetwork(AEWeights,input_h, input_w, input_d,numClasses,opts,nH_l1,nH_l2,fixNumLayers);


net.meta.classes.name = arrayfun(@(x)sprintf('%d',x),1:numClasses,'UniformOutput',false) ;

% --------------------------------------------------------------------
%                                                                Train
% --------------------------------------------------------------------

switch opts.networkType
    case 'simplenn', trainfn = @cnn_train ;
    case 'dagnn', trainfn = @cnn_train_dag ;
end


opts.expDir = expDir;

[net, info] = trainfn(net, imdb, getBatch(opts), ...
    'expDir', opts.expDir, ...
    net.meta.trainOpts, ...
    opts.train, ...
    'val', find(imdb.images.set == 2));%,...
end
% --------------------------------------------------------------------
function fn = getBatch(opts)
% --------------------------------------------------------------------
switch lower(opts.networkType)
    case 'simplenn'
        fn = @(x,y) getSimpleNNBatch(x,y) ;
    case 'dagnn'
        bopts = struct('numGpus', numel(opts.train.gpus)) ;
        fn = @(x,y) getDagNNBatch(bopts,x,y) ;
end
end
% --------------------------------------------------------------------
function [images, labels] = getSimpleNNBatch(imdb, batch)
% --------------------------------------------------------------------

images = imdb.images.data(:,:,:,batch) ;
labels = imdb.images.labels(1,batch) ;
end
% --------------------------------------------------------------------
function inputs = getDagNNBatch(opts, imdb, batch)
% --------------------------------------------------------------------
images = imdb.images.data(:,:,:,batch) ;
labels = imdb.images.labels(1,batch) ;
if opts.numGpus > 0
    images = gpuArray(images) ;
end
inputs = {'input', images, 'label', labels} ;
end

function net= initializeNetwork(AEWeights,input_h, input_w, input_d,numClasses,opts,nH_l1,nH_l2,fixNumLayers)

opts.batchNormalization = true ;

net.layers = {} ;

h1   = 1; w1=input_w;  d1= input_d; n1= nH_l1;
h2   = 1 ; w2=1; d2=nH_l1;  n2=nH_l2;
h3 = input_h; w3=1;d3=n2; n3= numClasses;


weights_l1  =  AEWeights.weights_init_l1;
biases_l1  =  AEWeights.biases_init_l1;

net.layers{end+1} = struct('type', 'conv', ...
    'weights', {{single(weights_l1),single(biases_l1)}}, ...
    'stride', [h1,w1], ...
    'pad', 0) ;

net.layers{end+1} = struct('type', 'relu') ;



sig= n2;
net.layers{end+1} = struct('type', 'conv', ...
    'weights', {{single(normrnd(0,sqrt(2/sig), [h2,w2,d2,n2])), zeros(1,n2,'single')}}, ...
    'stride', [h2,w2], ...
    'pad', 0) ;

net.layers{end+1} = struct('type', 'relu') ;
% 
sig= input_h*n2;

net.layers{end+1} = struct('type', 'conv', ...
    'weights', {{single(normrnd(0,sqrt(2/sig), [h3,w3,d3,n3])), zeros(1,n3,'single')}}, ...
    'stride', [h3,w3], ...
    'pad', 0) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

net.layers{end+1} = struct('type', 'softmaxloss','class',1:numClasses) ;


if fixNumLayers==1
    net.layers{1}.learningRate=[0 0];
elseif  fixNumLayers==2 
    net.layers{1}.learningRate=[0 0];
     net.layers{3}.learningRate=[0 0];
end

% Meta parameters
net.meta.inputSize = [input_h,input_w, input_d] ;

net.meta.trainOpts.learningRate = linspace(0.025,0.001,25);

net.meta.trainOpts.batchSize = 10;
net.meta.trainOpts.numEpochs = numel(net.meta.trainOpts.learningRate) ;



% Fill in default values
net = vl_simplenn_tidy(net) ;

% Switch to DagNN if requested
switch lower(opts.networkType)
    case 'simplenn'
        % done
    case 'dagnn'
        net = dagnn.DagNN.fromSimpleNN(net, 'canonicalNames', true) ;
        net.addLayer('top1err', dagnn.Loss('loss', 'classerror'), ...
            {'prediction', 'label'}, 'error') ;
        net.addLayer('top5err', dagnn.Loss('loss', 'topkerror', ...
            'opts', {'topk', 5}), {'prediction', 'label'}, 'top5err') ;
    otherwise
        assert(false) ;
end
end


function net = insertBnorm(net, l)
% --------------------------------------------------------------------
assert(isfield(net.layers{l}, 'weights'));
ndim = size(net.layers{l}.weights{1}, 4);
layer = struct('type', 'bnorm', ...
    'weights', {{ones(ndim, 1, 'single'), zeros(ndim, 1, 'single')}}, ...
    'learningRate', [1 1 0.05], ...
    'weightDecay', [0 0]) ;
net.layers{l}.biases = [] ;
net.layers = horzcat(net.layers(1:l), layer, net.layers(l+1:end)) ;
end



function imdb  = getFFImdb(data, labels,input_h,input_w,input_d,numClasses)
l1= floor(numel(labels)*0.8);

l2= length(labels)-l1;

set = [ones(1,l1), 2*ones(1,l2)];
if ndims(data)==4
    data= single(data);
else
    data = single(reshape(data,input_h,input_w,input_d,[]));
end

dataMean = mean(data(:,:,:,set == 1), 4);
data = bsxfun(@minus, data, dataMean) ;

imdb.images.data = data;
imdb.images.data_mean = dataMean; 
imdb.images.labels = labels;
imdb.images.set = set ;
imdb.meta.sets = {'train', 'val', 'test'} ;
imdb.meta.classes = arrayfun(@(x)sprintf('%d',x),1:numClasses,'uniformoutput',false) ;

end