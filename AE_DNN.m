%% Script to show steps for decoding finger flexion from ECoG data using autoencoder initialized CNNs
% Tejaswy Pailla
% 10/25/2017

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Loads data and label files for Finger Flexion data set
% --- Set what channels to use here in UseElectrodes

% Performs preprocessing if necessary (Line noise removal and CAR)

% Data restructuring : Creates trials in fxt format , i.e. #Freq Bins x #Time Bins
% --- Set the parameters for data restructuring using params

% Performs N fold CV by first randomizing trials, splitting into 5 parts, train on 4 , test on 1.
% Train Autoencoders using each channel as a trial using training features
% only and Uses the Autoencoder weights to initialize or fix NN layers

% Trains the Neural Network by using CNN implementation from MatConvNet and saves the net
% which minimizes Validation set Cost function

% This is not optimal as we have to give learning rate
% train CNN with 4:1 ratio of train and validation sets and using  gradient descent with batchsize=10.
% and batch normalization 
% we are subtracting mean from train and test features , so that network will have inputs around mean 0.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(0)
dataFolder = './data/' ; % folder with finger flexion data from Kail Miller's stanford repo
addpath(genpath('./util')); % folder with utility files

% set up matconvnet 
matconvpath= './matconvnet-1.0-beta19';  % or add path to matconvnet folder
run(fullfile(matconvpath, 'matlab', 'vl_setupnn.m')) ;

% only the following 5 subjects are used from Kai's dg dataset (finger flexion task)
subjects      = {'bp','wc','cc','jc','zt'};
arrLocation   = {'l','l','r','l','l'}; % their arrayLocations from table in documentation
SubCode       = {'A','B','C','D','E'};  % Use SubCodes for figures when publishing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% parameters to set


% Preprocessing steps setting
doLineNoiseRemoval    = 0; % not necessary since we remove frequency bins with line noises
doCAR                 = 1; % common average referencing
UseElectrodes         = [1 3 4];  %SMC  % Read doc to know the labels of regions


%
fs             =  1000;  % sampling Frequency
freqRange      =   4:4:140;
numfreqBandBins= length(freqRange)-1;

numCVFolds     =   5;


params.freqRange  = freqRange;
params.noLine     = 1;
params.numBands   = length(freqRange)-1;
params.numBins    = 6;
params.binWidth   = 100; % in ms
params.toRand     = 1;
params.toNorm     = 1;     % Set 1 if each frequency filtered channel is to be zscored
% transition trials are the first and last windows in each flexion
params.leaveTransition  = 0; %transition trials are the beginning and end trials in each flexion 0 is dont leave , 1 is leave transition trials, 2 is use only transition trials
numClasses        = 5;


% For autoencoder

activationFunction    = 'satlin';
sparsityProportion    = 0.05;   % 0 is highly sparse, range [0-1]
useGPU                = true; 

fixNumLayers = 1;

% nH_li is the number of hidden nodes in layer i

for subIndex= 1:5
    
    % paths to store trained autoencoder weights and CNNs 
    SubweightsDir    = ['AE_weights/SMC/',subjects{subIndex}];
    SubnnDir         = ['NN/SMC/',subjects{subIndex}];
    
    
    
    
    %%%%%%%%%%%%%%%% Load raw data  file
    
    load([dataFolder,subjects{subIndex},'/',subjects{subIndex},'_fingerflex.mat'])
    
    %%%%%%%%%%%%%%% Load Labels
    % Discrete labels are obtained from processed finger flexion traces
    % after thresholding and I made the FFLabels in the format numTrials x 4 where columns are :
    % [startIndex, stopIndex, finger labels (1 for rest and 2 to 6 for
    % fingers), movement/noMovement (1 or 2)]
    
    
    load(['./FFLabels/Labels_',subjects{subIndex},'.mat']);
    actualLabels =Labels(:,1:3);
    if numClasses==5 % remove rest trials
        ind  = find(actualLabels(:,3)~=1);
        Labels = actualLabels(ind,:);
        Labels(:,3)= Labels(:,3)-1;
    end
    
    %%%%%%%%%%% Format Ecog data:
    % raw data is in timesx ch, I use ch x time  in my util scripts. So do transpose
    
    if isempty(UseElectrodes)    %     Use all elec
        ecog  =  data';
        channelNumbers = 1:size(ecog,1);
    else
        ch   =  find(ismember(elec_regions,UseElectrodes));
        ecog = data(:,ch)';
        channelNumbers = ch;
    end
    
    %%%%%%%%%%%%%%%% Preprocess data
    if doLineNoiseRemoval
        opts.p=0.005; opts.bandwidth= 10; opts.tau=50;
        [ecog, ~, ~,~,~, ~, ~] =  LineNoiseFilter(ecog,fs,opts);
    end
    if doCAR
        ecog = CAR(ecog);
    end
    
    
    %%%%%%%%%%%%%%%% Restructure data for Autoencoder into time-freq
    %%%%%%%%%%%%%%%% format
    [features, classLabels] = restructureData(ecog, fs, Labels, params);
    % these features are nCh x numBands xnumBins xnumTrials
    
    
    [numChannels,nF,nB,numTrials] = size(features);
      
      % Make CV folds and train
      % Do not use cvpartition function to get indices for training and
      % testing. We donot want to have trials in the same flexion epoch split in
      % training and test sets .. since ecog data might be correlated, So do CV folds partition manually
    finalnumEpochs =   numTrials-mod(numTrials,numCVFolds);
    Ind            = reshape(1:finalnumEpochs,[], numCVFolds);
    
    for cv =1: numCVFolds
        
     
        testInd = Ind(:,cv);
          
        trainInd = setdiff(1:finalnumEpochs,testInd);
        
        testFeatures = features(:,:,:,testInd);
        testLabels = classLabels(testInd);
     
        trainFeatures = features(:,:,:,trainInd);
        trainLabels = classLabels(trainInd);
        
        for nH_l1 = 10:10:100
            
            
            if ~(isequal(exist( SubweightsDir, 'dir'),7)) % to save AE  weights
                mkdir(SubweightsDir);
            end
            
            
            trainedWeightsFile_l1 = [SubweightsDir,'/cv_',num2str(cv),'_n1_',num2str(nH_l1),'.mat'];
            trainingAE(trainFeatures,nH_l1,activationFunction,sparsityProportion,useGPU,trainedWeightsFile_l1)
            
            for nH_l2=10:10:50
                SubNNDir =[SubnnDir,'/n1_',num2str(nH_l1),'/n2_',num2str(nH_l2)]; % this will save the net
                if ~(isequal(exist(SubNNDir, 'dir'),7))
                    mkdir(SubNNDir);
                end
                
                % initialize weights layer 1
                weights_init_l1     = zeros(1,nF,nB,nH_l1);    biases_init_l1    = zeros(1,nH_l1);
                load(trainedWeightsFile_l1,'weights','biases');
                weights_init_l1(:,:,:,1:nH_l1) =   weights;         biases_init_l1(:,1:nH_l1)      =  biases;
                
                % initialize weights layer 2
                weights_init_l2     = [];%zeros(1,1,nH_l1,nH_l2);
                biases_init_l2   =[]; % zeros(1,nH_l2);
                
                AEWeights.weights_init_l1  = weights_init_l1;
                AEWeights.weights_init_l2  = weights_init_l2;
                AEWeights.biases_init_l1  = biases_init_l1;
                AEWeights.biases_init_l2  = biases_init_l2;
                
                
                
                [net_minObj, confusionMatrix_minObj,info,net_minErr, confusionMatrix_minErr,predLabels]= ...
                    trainNN(AEWeights,numClasses,trainFeatures, testFeatures,trainLabels,testLabels,...
                    nH_l1,nH_l2,fixNumLayers);
                
                save([SubNNDir,'/cv',num2str(cv),'.mat'],'info','confusionMatrix_minErr','net_minErr',...
                    'confusionMatrix_minObj','net_minObj');
            end
        end
        end
        
         %
    clear features classLabels ecog actualLabels Labels data classLabels
    
    end
   
    
