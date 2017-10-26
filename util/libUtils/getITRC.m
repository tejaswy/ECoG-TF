function [Pc,I] = getITRC(Pyc,Pc)
%calculate ITRC
%inputs:
% Pyc= transpose of confusion matrix  (vertical Y, horizontal X)
% Pc = input probability distribution: default uniform
% outputs:
% I = ITRC i..e. Blahut arimito
% Pc = maximizing input prob distribution



[N,M]=size(Pyc);% Pyc is transpose of conf matrix
if nargin<2 || isempty(Pc)
    Pc=(1/M).*ones(M,1);
end
% Pc should be column vector
Pc= Pc(:);


%normalize
numCat = length(Pc);
Pyc =Pyc./repmat(sum(Pyc,1),numCat,1);
Pc=Pc./sum(Pc);

t= find(Pyc==0); Pyc(t)= 0.00000001;

count=0;
e=0.00001;
Iu=1;Il=0;

Qy=Pyc*Pc;
F=zeros(1,M);
while (Iu-Il)>e
     count=count+1;
    for j=1:M
        temp=0;
        for k=1:N
            temp=temp+((Pyc(k,j))*log((Pyc(k,j))/Qy(k)));
        end
        F(j)=exp(temp);
    end
    x=F*Pc;
    Il=log2(x);
    Iu=log2(max(F));
    
    if (Iu-Il)<e
       
        Cc= {Il};
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % to calculate mutual information
        temp2=0;
        for j=1:M
            temp1=0;
            for k=1:N
                temp1=temp1+((Pyc(k,j))*log2((Pyc(k,j))/Qy(k)));
            end
            temp2=temp2+Pc(j)*temp1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        I=temp2; %mutual information
        %Pc={Pc};
        break;
    else
        if count>10000
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % to calculate mutual information
            temp2=0;
            for j=1:M
                temp1=0;
                for k=1:N
                    temp1=temp1+((Pyc(k,j))*log2((Pyc(k,j))/Qy(k)));
                end
                temp2=temp2+Pc(j)*temp1;
               
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            I=temp2; %mutual information
           % Pc={Pc};
             break
        end
            
            Pc=(1/x).*F'.*Pc;
            Qy=Pyc*Pc;
 
    end
end



end