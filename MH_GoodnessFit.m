function  [Pearson,NRMSE,R2, adjR2, VIP]= my_GoodnessFit(Y,Yhat,k,W)
% Goodness of Fit function
% INPUT:  Y: time*dimension (true trajectory)
%         Yhat: time*dim  (predicted trajectory)
%         k: number of features
%         W: samples*dim (W = [W_x, W_y, W_z, W_vx , and so on])  
%         k = the number of predictors
% this function calculates
% pearson correlation
% NRMSE
%the coefficient of determination. R2
% adjusted R2
% set W=0 if the method does not reqiure W


n = size(Y,1); % number of observations
Ybar = mean(Y);
SSE = sum((Y - Yhat).^2);
SST = sum((Y - Ybar).^2);

NRMSE = sqrt(SSE/n)./(max(Y)-min(Y));
Pearson = diag(corr(Y,Yhat))';

if W~=0
[J, F] = size(W);
J = J-1;
ExplainedSST = sum((Yhat - Ybar).^2); % for each variable
VIP = sqrt(sum(W(2:end,:).^2.*ExplainedSST,2)*(J/(F*sum(ExplainedSST))));
else
   VIP = 0; 
end
R2 = 1 - (SSE./SST);
beta = (n-1)/(n-k-1);
adjR2 = 1 - beta*(SSE./SST);


end