function I = calcMI(Pyc,Pc)
% Calculate Mutual Information I(X;Y) = H(X)-H(X/Y)
%inputs:
% Pyc= transpose of confusion matrix  (vertical Y, horizontal X)
% Pc = input probability distribution: default uniform
% outputs:
% I = ITR i..e. mutual info

[N,M]=size(Pyc);

if nargin<2 || isempty(Pc)
    Pc=(1/M).*ones(M,1);
end
% Pc should be column vector
Pc= Pc(:);
%Normalize Pyc
numCat = length(Pc);
Pyc =Pyc./repmat(sum(Pyc,1),numCat,1);
Pc=Pc./sum(Pc);
t= find(Pyc==0); Pyc(t)= 0.00000001;

% to calculate mutual information
Qy=Pyc*Pc;

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

end