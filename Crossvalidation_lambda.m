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


lambda_tune=0:1:10;
MSE_lambda_tune=NaN([1 width(lambda_tune)]);
MSE_q=NaN([1 width(Z1)]);


for tune = 1:width(MSE_lambda_tune)

for q = 1:width(MSE_q)

Z1_q=Z1;
Z1_q(:,q)=[];
Z2_q=Z2;
Z2_q(:,q)=[];

mu1_q=mean(horzcat(Z1_basis,Z1_q),2);
sigma1_q=std(horzcat(Z1_basis,Z1_q),[],2);
mu2_q=mean(Z2_q,2);
sigma2_q=std(Z2_q,[],2);

data_array=NaN([2 height(GCM_Trend_flat) width(Z1_q)]);

data_array(1,:,:)=(Z1_q-mu1_q)./sigma1_q;
data_array(2,:,:)=(Z2_q-mu2_q)./sigma2_q;

addpath(genpath('src'));

p = size(data_array,1);
n = size(data_array,2);
m = size(data_array,3);

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

[~,S,~] = svd(reshaped_data_array,'econ');

Cum_variance=NaN([height(S) 1]);
for j=1:height(S)
    S_diags=diag(S)/sum(diag(S));
    Cum_variance(j,1)=sum(S_diags(1:j));
end    

Var99_explained=find(Cum_variance>=0.995);

[U,~,~] = svd(reshaped_data_array,'econ');


smallPhi= U(:,1:Var99_explained(1));
l = size(smallPhi,2);

%% nugget estimate

% constant diagonal Q

x0 = [1,1];
lb = [0 0];
ub = [50 50];
nuggetestimates= zeros([p 1]);
smoothnessestimates = zeros([p 1]);

for i=1:p
    tt = smallPhi'* squeeze(data_array(i,:,:));
    smallPhi_S_Phi = tt*tt'/m;
    trS = norm(squeeze(data_array(i,:,:)),"fro")^2/m;
    [xval, fval]= fmincon(@(x) nuggetlikelihood_orthog(x,smallPhi_S_Phi,trS,n),x0,[], [], [], [], lb, ub);
    nuggetestimates(i) = xval(1);
    smoothnessestimates(i) = xval(2);
end


inv_error_variances = 1./nuggetestimates;


projected_data = zeros(p*l,m);

for i=1:m
projected_data(:,i)= vec((inv_error_variances .* data_array(:,:,i)) * smallPhi);
end


%% main algorithms

% no fusion
lambda=lambda_tune(1,tune);
graphguess = BGLblockdiag_orthogonal(lambda,inv_error_variances,projected_data);

% fusion
rho=0;
fusedguess = FBGL_orthogonal(lambda,rho,inv_error_variances,projected_data,graphguess);

%Y1

Y1=(Z1(:,q)-mu1_q)./sigma1_q;
Trend2=Trend(:,q+Basis_months(end));
% Obtain GLS estimate for W1=(W11,W12,....,W1l)


%A1*W=W1
A1_start=eye(l);
A1= NaN([l p*l]);
for r = 1:l
col1=2*r-1;
col2=2*r;
A1(1:l,col1:col2) = horzcat(A1_start(1:l,r),zeros(l,1));
end

Q=num2cell(fusedguess,[1,2]);
Q=blkdiag(Q{:});

%calculate Q_inv
fusedguess_inv=NaN([p p l]);

for r = 1:l
fusedguess_inv(:,:,r) = inv(fusedguess(:,:,r));
end

Q_inv=num2cell(fusedguess_inv,[1,2]);
Q_inv=blkdiag(Q_inv{:});

Inv=Q+transpose(smallPhi*A1)*smallPhi*A1/nuggetestimates(1);
Phi_Sigma_Y1_inv_Y1=transpose(smallPhi)*Y1/nuggetestimates(1)-transpose(smallPhi)*smallPhi*A1*inv(Inv)*transpose(smallPhi*A1)*Y1/nuggetestimates(1)^2;
Phi_Sigma_Y1_inv_Phi=transpose(smallPhi)*smallPhi/nuggetestimates(1)-transpose(smallPhi)*smallPhi*A1*inv(Inv)*transpose(smallPhi*A1)*smallPhi/nuggetestimates(1)^2;

W1_hat=inv(Phi_Sigma_Y1_inv_Phi)*Phi_Sigma_Y1_inv_Y1;

Exp_W2_given_Y1=NaN([l 1]);

for r = 1:l
Ql_inv=fusedguess_inv(:,:,r);
Exp_W2_given_Y1(r,1) = Ql_inv(1,2)*W1_hat(r,1)/Ql_inv(2,2);
end

Predicted_Z2=mu2+((smallPhi*Exp_W2_given_Y1).*sigma2);

MSE_q(:,q)=mean((Z2(:,q)-Predicted_Z2).^2);

end
MSE_lambda_tune(:,tune)=mean(MSE_q);

end

plot(lambda_tune,MSE_lambda_tune);
xlabel('lambda');
ylabel('MSE');
title('Tuning-lambda');