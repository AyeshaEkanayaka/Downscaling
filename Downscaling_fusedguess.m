clc
clear all

%Import data from .csv files
%Contains a total of 185 months from 1914 to 2099
GCM_flat=table2array(readtable('Interp_GCMs_flat.csv'));
GCM_weighted=table2array(readtable('Interp_GCMs_weighted.csv'));

Trend=table2array(readtable('Trend_weighted.csv'));

GCM_Trend_flat=GCM_flat-Trend;
GCM_Trend_weighted=GCM_weighted-Trend;

%Create an array
Current_months=89:105;%Time period MUR available (20003 to 2019)
Basis_months=1:137;%GCMs Going from past to future from 1914 to 2050

Z1=GCM_Trend_weighted(:,Current_months);
Z1_basis=GCM_Trend_weighted(:,1:Basis_months(end));
Z2=GCM_Trend_flat(:,Current_months);

%Standardize using pixel-wise means and stds
mu1=mean(Z1_basis,2);
sigma1=std(Z1_basis,[],2);

mu2=mean(Z2,2);
sigma2=std(Z2,[],2);

% Introduction and Setup Data

% read in your data as a p x n x m 'data_array'
% number of variables x number of locations x number of realizations

data_array=NaN([2 height(GCM_Trend_flat) width(Current_months)]);

data_array(1,:,:)=(Z1-mu1)./sigma1;
data_array(2,:,:)=(Z2-mu2)./sigma2;


p = size(data_array,1);
n = size(data_array,2);
m = size(data_array,3);


addpath(genpath('src'));
%% Basis Setup

% construct your orthogonal basis matrix 'smallPhi' here
% we use EOFs 
% be sure to save it if your basis setup is costly like this SVD
% smallPhi dimension is n x l where l = 2000 is the number of basis functions
% need smallPhi' * smallPhi = eye(l)

reshaped_data_array = zeros([n p*m]);
for i = 1:n
    reshaped_data_array(i,:) = vec(data_array(:,i,:));
end

%Include more months for basis calculation
reshaped_data_array=horzcat((Z1_basis-mu1)./sigma1,reshaped_data_array);

%Choose l such that the total variance explained is 99.5%
[~,S,~] = svd(reshaped_data_array,'econ');

Cum_variance=NaN([height(S) 1]);
for j=1:height(S)
    S_diags=diag(S)/sum(diag(S));
    Cum_variance(j,1)=sum(S_diags(1:j));
end    
plot(1:height(S),Cum_variance);
Var99_explained=find(Cum_variance>=0.995);
xline(Var99_explained(1),'--r');
title('Percentage of variance explained');
xlabel('EOF');
hold on
scatter(Var99_explained(1),Cum_variance(Var99_explained(1)),'r')
hold off
text(Var99_explained(1)+0.5,Cum_variance(Var99_explained(1))-0.01,['(' num2str(Var99_explained(1)) ',', num2str(Cum_variance(Var99_explained(1))) ')']);


[U,~,~] = svd(reshaped_data_array,'econ');

clear reshaped_data_array

smallPhi= U(:,1:Var99_explained(1));
l = size(smallPhi,2);

%% nugget estimate

% constant diagonal Q

x0 = [1,1];
lb = [0 0];
ub = [50 50];
nuggetestimates = zeros([p 1]);
smoothnessestimates = zeros([p 1]);

for i=1:p
    tt = smallPhi'* squeeze(data_array(i,:,:));
smallPhi_S_Phi = tt*tt'/m;
    trS = norm(squeeze(data_array(i,:,:)),"fro")^2/m;
    [xval, fval]= fmincon(@(x) nuggetlikelihood_orthog(x,smallPhi_S_Phi,trS,n),x0,[], [], [], [], lb, ub);
    nuggetestimates(i) = xval(1);
    smoothnessestimates(i) = xval(2);
end

% exponentially varying diagonal Q for each level
% 
% x0 = [1,0.05,1];
% lb = [0 0 0];
% ub = [50 0.1 50];
% 
% nuggetestimates_exp = zeros([p 1]);
% smoothnessestimates_exp = zeros([p 1]);
% marginalvarianceestimates_exp = zeros([p 1]);
% 
% for i=1:p
%     tt = smallPhi'* squeeze(data_array(i,:,:));
%     smallPhi_S_Phi = tt*tt'/m;
%     trS = norm(squeeze(data_array(i,:,:)),"fro")^2/m;
%     [xval, fval]= fmincon(@(x) nuggetlikelihood_expdiag_orthog(x,smallPhi_S_Phi,trS,n),x0,[], [], [], [], lb, ub);
%     nuggetestimates_exp(i) = xval(1);
%     smoothnessestimates_exp(i) = xval(2);
%     marginalvarianceestimates_exp(i) = xval(3);
% end

% see supplementary material for more explanation here
% nugget estimates were the same out to several decimal places
% for both constant diagonal and exponentially varying Q
% in our analysis

inv_error_variances = 1./nuggetestimates;

%% projecting data and creating matrix

% this is creating the matrix
% Phi' * Dinv * Y_datamatrix
% with multiplications dimensions (pl x pn) * (pn x pn) * (pn x m)
% efficiently using a Kronecker product trick

projected_data = zeros(p*l,m);

for i=1:m
projected_data(:,i)= vec((inv_error_variances .* data_array(:,:,i)) * smallPhi);
end


%% main algorithms

% no fusion
lambda=0.01;
graphguess = BGLblockdiag_orthogonal(lambda,inv_error_variances,projected_data);

% fusion
rho=1;
fusedguess = FBGL_orthogonal(lambda,rho,inv_error_variances,projected_data,graphguess);

% mle, no penalty
mleguess = BGLblockdiag_nopenalty(inv_error_variances,projected_data);


%Load Y1
%MUR location information 
%(column 3 and 4 are longitudes and latitudes respectively)
MUR_attribs=table2array(readtable('MUR_attribs.csv'));

%Total number of future months
t=width(GCM_Trend_weighted)-Current_months(end);

%Save predictions and predictive stds
%First 3 columns to store longitude,latitude and coral mask
Predicted_temp=NaN([n t+3]);
Predicted_temp(:,1:3)=MUR_attribs(:,3:5);
Predicted_sd=NaN([n t+3]);
Predicted_sd(:,1:3)=MUR_attribs(:,3:5);

%First need to obtain GLS estimate for W1=(W11,W12,....,W1l)
%Matrix A1 required for linear transformation of W to obtain W1. i.e. A1*W=W1
A1_start=eye(l);
A1= NaN([l p*l]);
for r = 1:l
col1=2*r-1;
col2=2*r;
A1(1:l,col1:col2) = horzcat(A1_start(1:l,r),zeros(l,1));
end

%Q matrix
Q=num2cell(fusedguess,[1,2]);
Q=blkdiag(Q{:});

%Calculate Q_inv
fusedguess_inv=NaN([p p l]);

for r = 1:l
fusedguess_inv(:,:,r) = inv(fusedguess(:,:,r));
end

Q_inv=num2cell(fusedguess_inv,[1,2]);
Q_inv=blkdiag(Q_inv{:});

%Begin downscaling
Start=cputime;
tic
for month = 1:t

%New Z1 and the corresponding Trend
Z1_new=(GCM_Trend_weighted(:,month+Current_months(end))-mu1)./sigma1; %Standardize
Trend_new=Trend(:,month+Current_months(end));

%W1_hat
Inv=Q+transpose(smallPhi*A1)*smallPhi*A1/nuggetestimates(1);
Phi_Sigma_Z1_inv_Z1=transpose(smallPhi)*Z1_new/nuggetestimates(1)-transpose(smallPhi)*smallPhi*A1*inv(Inv)*transpose(smallPhi*A1)*Z1_new/nuggetestimates(1)^2;
Phi_Sigma_Z1_inv_Phi=transpose(smallPhi)*smallPhi/nuggetestimates(1)-transpose(smallPhi)*smallPhi*A1*inv(Inv)*transpose(smallPhi*A1)*smallPhi/nuggetestimates(1)^2;

W1_hat=inv(Phi_Sigma_Z1_inv_Phi)*Phi_Sigma_Z1_inv_Z1;

%E[Z2|Z1]=Phi*E[W2|Z1]
%Var[Z2|Z1]=tau^2*I+Phi'*Var[W2|Z1]*Phi

Exp_W2_given_Z1=NaN([l 1]);
Var_W2_given_Z1=NaN([l 1]);

for r = 1:l
Ql_inv=fusedguess_inv(:,:,r);
Exp_W2_given_Z1(r,1) = Ql_inv(1,2)*W1_hat(r,1)/Ql_inv(2,2);
Var_W2_given_Z1(r,1) = Ql_inv(1,1)-Ql_inv(1,2)*Ql_inv(2,1)/Ql_inv(2,2);
end
  
Predicted_temp(:,3+month)=Trend_new(:,1)+mu2+((smallPhi*Exp_W2_given_Z1).*sigma2);

for loc = 1:n
Predicted_sd(loc,3+month)=sqrt(smallPhi(loc,:).^2*Var_W2_given_Z1+nuggetestimates(2)) ;
end

end
toc 
csvwrite('Predicted_temp.fusedguess.csv',Predicted_temp);
csvwrite('Predicted_sd.fusedguess.csv',Predicted_sd);