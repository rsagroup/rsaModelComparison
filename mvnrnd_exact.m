function U=mvnrnd_exact(G,P); 
% function U=mvnrnd_exact(G,P); 
% Generates multivariate normal random data with 
% a sample second moment matrix G=(U*U')/P, which is EXACTLY G
K = size(G,1); 
U = normrnd(0,1,K,P);
E = (U*U'); 
Z = E^(-0.5)*U;   % Make random orthonormal vectors 
A = cholcov(G); 
if (size(A,1)>P)
    error('not enough columns to represent G'); 
end; 
U = A'*Z(1:size(A,1),:)*sqrt(P); 
