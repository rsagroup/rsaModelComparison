function [nlmlSum,dnlmlSum,wN,VN,nlml] = rsa_marglNonlin(theta,Model, Y, Sigma);
% function [nlmlSum,dnlmlSum,wN,VN,nlml] = rsa_marglNonlin(theta,Model, Y, Sigma);
% returns marginal likelihood (and derivatives) for nonlinear models with 2 levels of parameters: 
% 1.level are omega parameters that are subject-specific. These are
%    marginalised out 
% 2.level are group "hyper"-parameters that
%    The first Model.numComp are the prior variance on the omegas 
%    The next Model.numNonlin are nonlinear components 
% INPUT: 
% -theta: Hyper parameters consist of liniear and nonlinear parts
%
% -Model: Model structure determing 
%   Model.numComp       : Number of linearly seperable components (at least 1)
%   Model.numPrior      : Number of prior parameters on the component coefficients 
%   Model.numNonlin     : Number of nonlinear parameters 
%   Model.nonlinP0      : Starting value of nonlinear(mixing) parameters 
%   Model.constantParams: Cell array of additional parameters to function 
%   Model.fcn           : Function returning RDM and derivatives 
% 
% -Y: Data (i.e., distance) for all subjects 
%       Size: (number of distance pairs)x(number of subjects)
% -Sigma: Covariance matrix of distance for all subjects
%       Size: (number of distance pairs)x(number of distance pairs)x(number of subjects)
% 
% OUTPUT: 
%   nlml     is the returned value of the negative log marginal likelihood
%   dnlml    is a (column) vector of partial derivatives of the negative
%                 log marginal likelihood wrt each log hyperparameter
%   wN       is the posterior Model of the regression coefficients
%   VN       is the posterior Variance of the regression coefficients
%
% Joern Diedrichsen


% Prior covariance of the regression coefficients
logtheta = theta(1:Model.numPrior);  % Get the number precesion parameter for regressors 
nonlinP  = theta([1:Model.numNonlin]+Model.numPrior); 

% Prior variance-covariance matrix  
V0 = diag(exp(logtheta));              

% Get the linear design matrix and the derviates of the design matrix 
% in respect to the non-linear parameters 
M = feval(Model.fcn,nonlinP,Model.constantParams{:}); 
X = permute(M.RDM,[2 1 3]); 
dX = permute(M.dRDMdTheta,[2 1 3 4]); % ay
[N, numSubj] = size(Y);
[N, numReg,depthX] = size(X);
if (depthX ==1) 
    X=repmat(X,1,1,numSubj); % use the same model for all subjects
    dX = repmat(dX,1,1,1,numSubj); % ay
elseif (depthX~=numSubj)
    error('X must needs to be a matrix or have a size in the 3rd dimension of numSubj'); 
end;

% nlml    = zeros(numSubj);
% dnlml   = zeros(numReg+Model.numNonlin,numSubj);

for s=1:numSubj
    
    % Precompute outer products of the X-vectors 
    XX = zeros(N,N,numReg); 
    for i=1:numReg
        XX(:,:,i)=X(:,i,s)*X(:,i,s)'; 
    end; 

    S = Sigma(:,:,s) + X(:,:,s)*V0*X(:,:,s)' ;                              % Compute training set covariance matrix
    L = chol(S)';                                                           % Hack: Cholesky factorization of the covariance
    alpha = solve_chol(L',Y(:,s));                                          % Hack: Convenience function (=S^-1*Y)
    
    % Negative log-likihood
    nlml(s)  = 0.5*sum(sum(alpha.*Y(:,s),2)) + sum(log(diag(L))) + 0.5*N*log(2*pi);  
    
    % Derivatives
    if (nargout>1)
        invS     = (L'\(L\eye(N)));
        W        = alpha*alpha'-invS;                                       % This is 2 times derivative of negative log-likelihood wrt S (=2*dl/dS)
        
        % Derivative of negative log-likelihood wrt log prior parameters (=(dl/dS)x(dS/dV0)xV0)
        for i=1:numReg
            thetai = exp(logtheta(i));
            dnlml(i,s) = -1/2*sum(sum(W.*(XX(:,:,i))))*thetai;              
        end; 
        % Derivative of negative log-likelihood wrt nonlinear parameters
        % (ay) may need some hack to speed up
        dSdTheta = zeros(N,N);
        for i=1:Model.numNonlin
           for d=1:N
              for dd=1:N                  
                  dSdTheta(d,dd) = dX(d,:,i,s)*V0*X(dd,:,s)' + X(d,:,s)*V0*dX(dd,:,i,s)';           
              end
           end
           thetai = 1;%exp(nonlinP(i)); no need for log-transform here
           dnlml(Model.numPrior+i,s) = -1/2*sum(sum(W.*dSdTheta))*thetai;
        end 
    end;
    
    % Regression coefficients
    if (nargout>2)
        wN(:,s)  = V0*X(:,:,s)'*alpha;                                      % Regression coefficients by Matrix inversion
    end;
    if (nargout>3)
        VN(:,:,s)  = V0 - V0*X(:,:,s)'*invS*X(:,:,s)*V0;                    % Posterior variance
        VN(:,:,s) = (VN(:,:,s)+VN(:,:,s)')/2;                               % Prevent asymmrety through roundoff.
    end;
    %keyboard; 
end;
% Sum marginal likelihoods for optimization 
nlmlSum  = sum(nlml); 
if (nargout>1)
    dnlmlSum = sum(dnlml,2); 
end; 