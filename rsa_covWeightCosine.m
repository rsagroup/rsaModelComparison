function r=rsa_covWeightCosine(A,B,removeInt) % Covariance-Weighted cosine similarity measure 
% function r=rsa_covWEightCosine(A,B,removeInt); 
% Covariance-Weighted cosine similarity measure between distances 
% INPUT: 
%   A: N1xD matrix, with each row consistuting a vectorized RDM 
%   B: N2xD matrix, with each row constituting a vectorized RDM 
%   removeInt: Remove intercept 0: cosine Angle 1: correlation 
% Set defaults for removed Intercept 
if (nargin<3) 
    removeInt=0; 
end; 
persistent sq; % Define a persistent sq, to prevent repeated precomputation 
D=size(A,2); 
if (size(B,2)~=D) 
    error('A and B need to have the same number of columns'); 
end;

% If no fitting prewhitening matrix for the distances can be found in
% memory - make one 
if (size(sq,1)~=D)
    K = ceil(sqrt(2*D)); % (1 + sqrt(1+8*D))/2, but works for large n
    if K*(K-1)/2 ~= D
        error('bad input size'); 
    end; 
    C=pcm_indicatorMatrix('allpairs',[1:K]); 
    CC = C*C';
    Var= CC.*CC; 
    [V,L]=eig(Var); 
    l=real(diag(L));
    sq = V*bsxfun(@rdivide,V',sqrt(l)); % Slightly faster than sq = V*diag(1./sqrt(l))*V';
end; 
    
wA=A*sq; 
wB=B*sq; 
if (removeInt) 
    wA=bsxfun(@minus,wA,mean(wA,2));
    wB=bsxfun(@minus,wB,mean(wB,2));
end; 
wA=bsxfun(@rdivide,wA,sqrt(sum(wA.^2,2)));
wB=bsxfun(@rdivide,wB,sqrt(sum(wB.^2,2)));
r=wA*wB';
