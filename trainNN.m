
function  [net_minObj, confusionMatrix_minObj,info,net_minErr, confusionMatrix_minErr,best]=trainNN(AEWeights,numClasses,trainFeatures, testFeatures,trainLabels,testLabels,nH_l1,nH_l2,fixNumLayers)

datetime=datestr(now);
datetime=strrep(datetime,':','_'); %Replace colon with underscore
datetime=strrep(datetime,'-','_');%Replace minus sign with underscore
datetime=strrep(datetime,' ','_');%Replace space with underscore

expDir= datetime; 

[net, info,imdb] = nn_setup(AEWeights,trainFeatures, trainLabels,expDir,nH_l1,nH_l2,numClasses,fixNumLayers);


% % % % %  selecting network with min validation set objective
[val, minValObjEpoch]=  min([info.val.objective]);
ep              = num2str(minValObjEpoch);
netAddress      = fullfile(pwd,expDir,char(strcat('net-epoch-',ep,'.mat')));
net             = load(netAddress);
net              = net.net;

net              = cnn_ff_deploy(net);

dataMean = imdb.images.data_mean;
testFeatures= bsxfun(@minus, testFeatures,dataMean) ;

% show the classification result

error=0;

best = zeros(size(testLabels));
for i=1:length(testLabels)
    temp = single(testFeatures(:,:,:,i));
    res = vl_simplenn(net,temp);    
    scores = squeeze(gather(res(end).x)) ;    
    [bestScore, best(i)] = max(scores) ;
    
end
confusionMatrix= confusionmat(testLabels, best, 'order',1:5);

net_minObj= net;
confusionMatrix_minObj = confusionMatrix;


% % % % %  selecting network with min validation set Error

[val, minValerrEpoch]=  min([info.val.top1err]);
ep              = num2str(minValerrEpoch);
netAddress      = fullfile(pwd,expDir,char(strcat('net-epoch-',ep,'.mat')));
net             = load(netAddress);
net             = net.net;
net             = cnn_ff_deploy(net);

dataMean = imdb.images.data_mean;
testFeatures = bsxfun(@minus, testFeatures,dataMean) ;

% show the classification result

error=0;

best = zeros(size(testLabels));
for i=1:length(testLabels)
    temp = single(testFeatures(:,:,:,i));
    res = vl_simplenn(net,temp);
    scores = squeeze(gather(res(end).x)) ;
    [bestScore, best(i)] = max(scores) ;
    
end
confusionMatrix= confusionmat(testLabels, best, 'order',1:5);

net_minErr= net;
confusionMatrix_minErr = confusionMatrix;
rmdir(expDir,'s');
end

