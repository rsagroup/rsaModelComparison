function [R2,mCorr,Ypred,h]=encode_crossval(Y,Z,partition,method,varargin); 
% function [R2,corr,Ypred,h]=encode_crossval(Y,Z,partition,method,varargin)
% crossvalidation of an encoding model across the partitions
%
% INPUT: 
%    Y: NxP matrix Data 
%    Z: NxQ matrix for random effects 
%    partions: Nx1 vector of partitions
%    method: Method for pattern estimation 
%            'linregress':     Normal, unregularised linear regression 
%            'pcm_EM':         Bayesian (Ridge) regression, using EM to estimate constrained G matrix (see pcm_EM) 
%            'pcm_EM_free':    Bayesian (Ridge) regression, completely unstrained variances and covariances
%            'pcm_NR':         Bayesian (Ridge) regression, uses NR to estimate constraiend G matrix (see pcm_NR) 
%            'pcm_NR_diagfree: Bayesian (Ridge) regression, diagnonal matrix, but unqiue elements 
%            'ridge_fixed'     Bayesian (Ridge) regression with a fixed regularisation parameter (inv(G)*sigma2) 
%            'GCV'             Generalized cross-validation (golub et al., 1979) 
% VARARGIN: 
%    'X', X     : design matrix for fixed effects removal (also passed to
%                 function for ReML estimation) 
%    'Gc',Gc    : Variance components for pcm_NR
%    'Ac',Ac    : Variance components for pcm_EM
%    
% OUTPUT: 
%    R2         : the portion of correctly predicted variance of the betas
%    corr       : correlation between predicted and observed values 
%    Ypred      : 1*n vector with predicted classes
%    h          : Hyper-parameters (meaning depends on the pcm version used) 
% 
% (c) 2016 Joern Diedrichsen 

[N,P] = size(Y);      % 
meanS = 0; 
X     = []; 
Gc    = []; 
Gd    = []; 
G     = [];
sigma2 = []; 
Ac    = []; 
numIter = [];
h     = []; 
evalCrit = 'RSS'; 
vararginoptions(varargin,{'X','meanS','numIter','Gc','Ac','Gd','G','sigma2','evalCrit'});

% REmove X matrix from Y and Z
% -------------------------------------------------
if (~isempty(X))
    Y=Y-X*pinv(X)*Y; % remove X matrix 
    Z=Z-X*pinv(X)*Z; % remove X matrix 
end; 


% Remove the mean from the Y-data (across voxels) 
% -------------------------------------------------
if meanS==1
    a  =  pinv(Z)*sum(Y,2)/P;
    Y  =  bsxfun(@minus,Y,Z*a);
end; 

% Loop over partitions 
% -------------------------------------------------
part=unique(partition);  % Get unique partitions
for i=1:length(part)
    trainI = partition~=part(i);
    testI  = partition==part(i);
               
    Ytrain = Y(trainI,:);
    Ztrain = Z(trainI,:);    
    Ztest  = Z(testI,:);
    
    if (~isempty(X)) 
        Xtrain = X(trainI,:);
        Xtrain = Xtrain(:,find(sum(Xtrain.^2,1)~=0)); 
    else 
        Xtrain = []; 
    end; 
    
    % Determine the predicted patterns. For linear regression this is easy,
    % For Bayesian (ridge) regression, we determine the regularisation
    % parameter using PCM. 
    % -------------------------------------------------
    switch(method) 
        case 'linregress' 
            u = pinv(Ztrain)*Ytrain; 
        case 'pcm_EM' 
            [G,H,u]=pcm_EM(Ytrain,Ztrain,'Ac',Ac,'X',Xtrain); 
            h(i,:)=H(:,end)';
            h(i,1:end-1)=h(i,1:end-1).^2; % Square them for comparision purposes 
        case 'pcm_EM_free' 
            [G,H,u]=mvpattern_covcomp_free(Ytrain,Ztrain); 
            h(i,:)=H(:,end)';
        case 'pcm_NR' 
            [G,h(i,:),u]=pcm_NR(Ytrain,Ztrain,'Gc',Gc,'X',Xtrain);
        case 'pcm_NR_diagfree' 
            [G,h(i,:),u]=pcm_NR_diagfree(Ytrain,Ztrain,'X',Xtrain);
        case 'pcm_NR_diag' 
            [G,h(i,:),u]=pcm_NR_diag(Ytrain,Ztrain,'Gc',Gc,'Gd',Gd,'X',Xtrain);  
        case 'ridgeFixed' 
            % u1 = G*Ztrain'*((Ztrain*G*Ztrain'+eye(size(Ztrain,1))*sigma2)\Ytrain);
            u = (Ztrain'*Ztrain+pinv(G)*sigma2)\(Ztrain'*Ytrain); 
        case 'GCV' 
            [u,h(i,1)] = encode_gcv(Ytrain,Ztrain); 
    end; 
    
    % Generate prediction 
    % -------------------------------------------------
    Ypred(testI,:)=Ztest*u; 
end;

% Evaluate prediction by calculatin R2 
% -------------------------------------------------
SST=sum(sum(Y.*Y));
res=Y-Ypred; 
SSR=sum(sum(res.*res)); 
R2=1-SSR/SST;
mCorr = sum(sum(Y.*Ypred))/...
      sqrt(sum(sum(Y.*Y))*sum(sum(Ypred.*Ypred))); 
