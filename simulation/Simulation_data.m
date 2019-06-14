clear;
clc;

%% Generating beta
beta=zeros(1,2000);
beta(1)=1.5;
beta(2)=-1.2;
beta(3)=1.8;
beta(4)=-2;
beta(5)=2.5;
beta(6)=-1.2;
beta(7)=1;
beta(8)=-1.5;
beta(9)=2;
beta(10)=-1.6;
beta_t=beta';

%% Generating simulation data 
train_size=100;     % size = 100, 200
X1 = normrnd(0, 1, train_size, size(beta,2));
X2 = normrnd(0, 1, train_size, size(beta,2));
X3 = normrnd(0, 1, train_size, size(beta,2));
[n,p]=size(X1);

% Setting correlation
correlation = [0, 0.2, 0.4, 0.6, 0.8];
cor = correlation(1);             % correlation =0, 0.2, 0.4, 0.6, 0.8
for i=1:n
    for j=2:p
        X1(i,j) =  X1(1,1) * cor + X1(i,j) * (1-cor);
        X2(i,j) =  X2(1,1) * cor + X2(i,j) * (1-cor);
        X3(i,j) =  X3(1,1) * cor + X3(i,j) * (1-cor);
    end
end

noise = [0, 0.2, 0.4, 0.6, 0.8];
noise1 = [0.1, 0.2, 0.3];
sigm = noise(3);         % noise = 0, 0.2, 0.4, 0.6, 0.8
sigm1 = noise1(1);
sigm2 = noise1(2);
sigm3 = noise1(3);

l1 = X1 * beta' + sigm * normrnd(0, 1, n, 1);
prob1=exp(l1)./(1 + exp(l1));

l2 = X2 * beta' + sigm * normrnd(0, 1, n, 1);
prob2=exp(l2)./(1 + exp(l2));

l3 = X3 * beta' + sigm * normrnd(0, 1, n, 1);
prob3=exp(l3)./(1 + exp(l3));

for i=1:train_size
    if prob1(i)>0.5
        Y1(i)=1;
    else
        Y1(i)=0;
    end
    if prob2(i)>0.5
        Y2(i)=1;
    else
        Y2(i)=0;
    end
    if prob3(i)>0.5
        Y3(i)=1;
    else
        Y3(i)=0;
    end
end

Y1 = [Y1;ones(1,length(Y1))]';
Y2 = [Y2;2*ones(1,length(Y2))]';
Y3 = [Y3;3*ones(1,length(Y3))]';
label = [Y1;Y2;Y3];

dlmwrite('Dataset1_matrix.txt',X1,"\t");
dlmwrite('Dataset2_matrix.txt',X2,"\t");
dlmwrite('Dataset3_matrix.txt',X3,"\t");
dlmwrite('Dataset_label.txt',label,"\t");
