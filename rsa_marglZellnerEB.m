function [nlmlSum,dnlmlSum,wN,VN,nlml] = rsa_marglZellnerEB(logtheta, X, Y, Sigma);
% function [nlmlSum,dnlmlSum,wN,VN,nlml] = rsa_marglRidgeEB(logtheta, X, Y, Sigma);
% Marginal likelihood of the Baysean multiple regression model with known
% multivariate noise covariance (Sigma) 
% Integrates over different subjects (observations) which may have
% different design matrices and different Sigma's 
% 
% y_n = X_n*beta_n + noise_n; 
% 
% The prior variance of the regression coefficients is V0 = inv(X'*X) * exp(logtheta)
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
% Joern Diedrichsen%
% Joern Diedrichsen

[N, numReg,depthX] = size(X);
[N, numSubj] = size(Y);
if (depthX ==1) 
    X=repmat(X,1,1,numSubj); 
elseif (depthX~=numSubj)
    error('X must needs to be a matrix or have a size in the 3rd dimension of numSubj'); 
end;


for s=1:numSubj
    V0X = (X(:,:,s)'*X(:,:,s))\X(:,:,s)'*exp(logtheta);              % Prior covariance of the regression coefficients
    S = Sigma(:,:,s) + X(:,:,s)*V0X ;       % compute training set covariance matrix
    L = chol(S)';                      % cholesky factorization of the covariance
    alpha = solve_chol(L',Y(:,s));     % Convenience function
    
    nlml(s)  = 0.5*sum(sum(alpha.*Y(:,s),2)) + sum(log(diag(L))) + 0.5*N*log(2*pi);  % Negative log-likihood
    if (nargout>1)
        invS     = (L'\(L\eye(N)));
        W        = alpha*alpha'-invS;           % this is (alpha*alpha' - inv(S))
        dnlml(s) = -1/2*sum(sum(W.*(X(:,:,s)*V0X)));            % Derivative of L
    end;
    if (nargout>2)
        wN(:,s)  = V0X'*alpha;                            % regression coefficients by Matrix inversion
    end;
    if (nargout>3)
        VN(:,:,s)  = V0 - V0X'*invS*V0X';                   % Posterior variance
        VN(:,:,s) = (VN(:,:,s)+VN(:,:,s)')/2; % Prevent asymmrety through roundoff.
    end;
end;
% Sum marginal likelihoods for optimization 
nlmlSum  = sum(nlml); 
if (nargout>1)
    dnlmlSum = sum(dnlml);
end; 
