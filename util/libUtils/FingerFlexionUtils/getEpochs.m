function Labels = getEpochs(folderName,filename)

%output Labels will have size numTrials x 4; columns are startIndex,
%stopIndex, finger labels[1:6], no movement/ movemnet[1,2];

% Total Subjects =9
% cc_ddg.csv  jp_ddg.csv  wm_ddg.csv  ht_ddg.csv  mv_ddg.csv  zt_ddg.csv
% bp_ddg.csv  jc_ddg.csv  wc_ddg.csv
% subjects ht,jp,wm have bad data
if nargin<2
    filename   = 'cc_ddg.csv';
end

if nargin<1
    folderName = '/net/derData/FingerFlexion/DiscreteDataGlove';  
end



M  = csvread(fullfile(folderName,filename));

M  = reshape(M,[],5);



[samples, numFingers]= size(M);


% curFlag is the current value in the M matrix
% prevFlag is the value at prev time Pt

tempMatrix  =   zeros(size(M));

for fIdx= 1: numFingers
    prevFlag = M(1,fIdx);
    prevStop=1;
    for i= 1:samples
        curStart =i;
        curFlag  =  M(i,fIdx);
        
        
        if prevFlag~=curFlag
            if curFlag==1 && prevFlag==0
                
                tempMatrix(prevStop:curStart-1,fIdx)=0;
                prevStop= curStart;
            elseif prevFlag>1
                
                tempMatrix(prevStop:curStart-1,fIdx)= fIdx;
                prevStop= curStart;
            end
            
        end
        
        prevFlag=curFlag;
    end
end

scLabels =zeros(samples,1);
for i= 1:numFingers
    scLabels(tempMatrix(:,i)~=0)=i;
end


prevLabel= scLabels(1);
tstart = 1;
Labels =[];
for i=1:samples
    curLabel= scLabels(i);
    if prevLabel~=curLabel;
        Labels = [Labels; tstart,i-1,prevLabel];
        prevLabel= curLabel;
        tstart=i;
    end
    
end
% get movement and no movement periods
r= find(Labels(:,3)~=0);
Labels(r,4)=1;


% to print number of trials for each finger
% for i=1:5
% a= find(gtLabels(:,3)==i);
% length(a)
% end

Labels(:,3)= Labels(:,3)+1;
Labels(:,4)= Labels(:,4)+1;
end





