function [omega,varOmega,loglike,sigma,loglike2]=rsa_fitModelIRLS(Model,d,Sig,numPart,numVox,varargin)
% function omega=fitModelIRLS(Model,d,Sig,numPart,numVox,varargin)
% Fits the distance matrix using iteratively-reweighted-least-squares 
% INPUT 
%    d          : Crossvalidated Mahalanobis Distances, in a flat N x nDist vectorized representation
%                 with nDist = nCond*(nCond-1)/2 
%    X          : Component design matrix of model nDist x H matrix 
%    Sig        : Estimated noise-covariance matrix of the betas. Should be 
%                       a) a Nx1 vector (assuming independence of betas Sig <- eye(K)*Sig
%                       b) a Nx(KxK) matrix in flattened Sigma in matrix notation 
%                       c) vectorized covariance matrix with K*(K-1)/2 columns 
%                 partition (KxK), or a single scale for 
%    nPart      : Number of partitions (runs), used to calucalte d
%    nVox       : Number of voxels 
%  VARARGIN (USEROPTIONS) 
%    tol        : Tolerance for convergence (default 1e-05) 
%    maxIter    : Maximal number of iterations (default 100) 
%    priorVar   : Prior Variance of regression coefficients (defaults to inf)
% OUTPUT
%    omega      : Estimated maximum a-posteriori component weight coefficients N x H matrix 
%    varOmega   : Vectorised posterior variance of omega 
%    loglike    : Data Log-likelihood of the distances under the model &
%                 omega (not marginal log likelihood!) 
% (c) 2015 Joern Diedrichsen & Atsushi Yokoi 
% joern.diedrichsen@googlemail.com 
% --------------------------------
Opt.tol     = 1e-04;    % tolerance for updating
Opt.maxIter = 100;      % maximum iteration
Opt.priorVar = inf;    % Variance of the prior 
Opt.assumeVoxelInd = 1; % Assume that the overall size of V is exact? ? If not, it uses the residuals from the distances to get scaling factor sigma 

Opt = rsa_getUserOptions(varargin,Opt);

% Check the input size of distances
[numSubj,numDist]   = size(d); 
numCond   =  ceil(sqrt(2*numDist)); 
if numCond*(numCond-1)/2 ~= numDist
    error('Distance needs to have proper length for pdist (nCond*(nCond-1)/2)');
end

% Check input size of the noise-covariance 
[n,lengthSig]   = size(Sig); 
if (n~=numSubj) 
    error('Sig (noise covariance) needs as many rows as distances');
end; 
if (lengthSig~=1 && lengthSig~=numCond*(numCond-1)/2+numCond && lengthSig~=numCond*numCond)
    error('Sig (noise covariance) needs to have either 1 or nCond*(nCond-1)/2+nCond or nCond x nCond columns');
end; 

% Check input size of component design matrix 
if (isstruct(Model))
    X = vertcat(Model(:).RDM)'; 
else
    X = Model; 
end;

[m,numReg]   = size(X); 
if (m~=numDist) 
    error(message('Component design matrix needs to have M rows')); 
end; 

% Check input of nVox 
if (length(numVox)==1) 
    numVox=ones(numSubj,1)*numVox; 
elseif (length(numVox)~=numSubj) 
    error('nVox must have the same length as the number of RDMs or be 1'); 
end; 

% Generate distance contrast matrix for the distances 
C    = rsa_indicatorMatrix('allpairs',[1:numCond]); 

% Regularisation matrix 
L    = 1./Opt.priorVar*eye(numReg);

% Transfer to single for speed 
d= single(d); 
X= single(X); 
C= single(C); 
L= single(L); 
Sig = single(Sig); 

% Preallocate results for speed 
omega     = zeros(numSubj,numReg);       % Resulting estimated weights 
omegaIter = zeros(numReg,Opt.maxIter); 
% Start iterative fitting 
for n=1:numSubj 
    
    % Get the current noise variance matrix into the right form 
    if (lengthSig==1)   % Scalar noise variance 
        thisSig = eye(numCond)*Sig(n,1); 
    elseif (lengthSig ==numCond*(numCond-1)/2+numCond)
        thisSig                  = zeros(numCond);                % make matrix 
        thisSig(tril(true(numCond),0)) = Sig(n,:);                % Assign the lower triag
        transSig                 = thisSig';                % Now make symmetric           
        thisSig(thisSig==0)      = transSig(thisSig==0);    
    else 
        thisSig    = reshape(Sig(n,:),numCond,numCond); 
    end; 
    
    % Initial covariance matrix 
    omegaIter(:,1)=(X'*X)\(X'*d(n,:)');
    V   = rsa_varianceLDC(X*omegaIter(:,1),C,thisSig,numPart,numVox(n)); 

    % Iterate until convergence 
    for i=2:Opt.maxIter
        % General least square estimator. The following is faster than 
        % iV = pinv(V); 
        % omegaIter(:,i)=pinv(X'*iV*X + L)'*X'*iV*d(n,:)';                
        % omegaIter(:,i)=pinv(X'*(V\X) + L)'*X'*(V\d(n,:)');   
        iVX = V\X; 
        omegaIter(:,i)=(X'*iVX+L)\iVX'*d(n,:)'; 
        V = rsa_varianceLDC(X*omegaIter(:,i),C,thisSig,numPart,numVox(n));         % Update variance estimate 
    
        % check convergence: maximal deviation in proportion change 
        if i>1
            if (max(abs(1-(omegaIter(:,i)./omegaIter(:,i-1))))<Opt.tol)
                break;
            end
        end
        if i>=Opt.maxIter
            warning('iteration not converged.');
        end; 
    end;
    omega(n,:)    = omegaIter(:,i)'; 
    varOmega(n,1) = 0; 
    % Calculate Data loglikelihood
    prederr       = d(n,:)'-X*omega(n,:)';   % Prediction error of the model 
    if (Opt.assumeVoxelInd)   % Assume independence of voxels (or known dependence) 
        loglike(n,1)  = -numDist/2*log(2*pi)-1/2*logdet(V)-1/2*prederr'*(V\prederr); % Log likelihood of the data at maximum likelihood 
    else 
        sigma(n,1)    = prederr'*(V\prederr)/numDist;       % ML estimate 
        loglike(n,1)  = -numDist/2*log(2*pi)-numDist/2*log(sigma(n,1))-1/2*logdet(V)-1/(2*numDist); % Log likelihood of the data at maximum likelihood 
    end; 
    % fprintf('.'); 
end; 

    % fprintf('\n'); 


