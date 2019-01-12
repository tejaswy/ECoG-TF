
function trainingAE(features,nH_l1,activationFunction,sparsityProportion,useGPU,trainedWeightsFile_l1)


[numChannels,nF,nB,numTrials] = size(features);

count=1;
for trInd=1:numTrials
    for chInd=1:numChannels
        xTrain{count}= squeeze(features(chInd,:,:,trInd));
       
        count=count+1;
    end
end


autoenc1 = trainAutoencoder(xTrain,nH_l1,'EncoderTransferFunction',activationFunction,...
    'DecoderTransferFunction',activationFunction,'SparsityProportion',sparsityProportion,...
    'ScaleData', false, 'ShowProgressWindow',false,'UseGPU',useGPU,'MaxEpochs',200);

weights   =      autoenc1.EncoderWeights  ; %size  numKernelsx168;
biases    =      autoenc1.EncoderBiases;
weights   =      reshape(weights',nF, nB, nH_l1);

save(trainedWeightsFile_l1,'autoenc1','weights','biases');

end
