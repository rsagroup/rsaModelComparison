function [nlmlSum,dnlmlSum,wN,VN,nlml] = rsa_marglRidgeIndividEB(logtheta, X, Y, Sigma);
% function [nlmlSum,dnlmlSum,wN,VN,nlml] = rsa_marglRidgeIndividEB(logtheta, X, Y, Sigma);
% Marginal likelihood of the multiple regression model with known
% multivariate noise covariance (Sigma) 
% Inetgrates over different subejcts (observations) which may have
% different design matrices and different Sigma's 
% 
% y_n = X_n*beta_n + noise_n; 
% 
% The prior variance of the regression coefficients is V0 = diag(exp(logtheta))
% 
% INPUT: 
%   logtheta : log prior variance on the betas 
%   X        : Design matrix (N x Q (x numSubj)  
%   Y        : Data matrix (N x numSubj) 
%   Sigma    : NxN noise Matrix (x numSubj) 
% OUTPUT 
%   nlmlSum  : the negative log marginal likelihood (summed over subjects)              
%   dnlmlSum :  a (column) vector of partial derivatives of the negative
%                 log marginal likelihood wrt each log hyperparameter
%   wN       : the posterior mean of the regression coefficients
%   VN       : the posterior Variance of the regression coefficients
%   nlml     : negative log likelihood for each Subject (for model comp) 
% Joern Diedrichsen

[N, numReg,depthX] = size(X);
[N, numSubj]       = size(Y);

% Deal with possible different design matrices for different subjects 
if (depthX ==1) 
    X=repmat(X,1,1,numSubj); 
elseif (depthX~=numSubj)
    error('X must needs to be a matrix or have a size in the 3rd dimension of numSubj'); 
end;

% Prior covariance of the regression coefficients
V0 = diag(exp(logtheta));              

for s=1:numSubj
    
    % Precompute outer products of the X-vectors 
    XX = zeros(N,N,numReg); 
    for i=1:numReg
        XX(:,:,i)=X(:,i,s)*X(:,i,s)'; 
    end; 

    S = Sigma(:,:,s) + X(:,:,s)*V0*X(:,:,s)' ;       % compute training set covariance matrix
    L = chol(S)';               % cholesky factorization of the covariance
    alpha = solve_chol(L',Y(:,s));   % Convenience function
    
    nlml(s)  = 0.5*sum(sum(alpha.*Y(:,s),2)) + sum(log(diag(L))) + 0.5*N*log(2*pi);  % Negative log-likihood
    if (nargout>1)
        invS     = (L'\(L\eye(N)));
        W        = alpha*alpha'-invS;           % this is (alpha*alpha' - inv(S))
        for i=1:numReg
            thetai = exp(logtheta(i));
            dnlml(i,s) = -1/2*sum(sum(W.*(XX(:,:,i))))*thetai;              
        end; 
    end;
    if (nargout>2)
        wN(:,s)  = V0*X(:,:,s)'*alpha;                            % regression coefficients by Matrix inversion
    end;
    if (nargout>3)
        VN(:,:,s)  = V0 - V0*X(:,:,s)'*invS*X(:,:,s)*V0;                   % Posterior variance
        VN(:,:,s) = (VN(:,:,s)+VN(:,:,s)')/2; % Prevent asymmrety through roundoff.
    end;
end;
% Sum marginal likelihoods for optimization 
nlmlSum  = sum(nlml); 
if (nargout>1)
    dnlmlSum = sum(dnlml,2); 
end; 