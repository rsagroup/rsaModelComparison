function dist_pw = rdm_prewhitenDist(dist,V,varargin)
%function dist_pw = rdm_prewhitenDist(dist,V)
% prewhitens the distances using the covariance matrix of the vector
% distances (V)
% 
% INPUT:  - dist:   vector containing distances (1 x d) vectorized form of rdm - upper triangular, use rsa_vectorizeRDM)
%         - V:      covariance matrix of the vector distances (d x d)
%
% OUTPUT: - dist_pw: vector containing prewhitened distances (1 x d)
%
% VARARGIN: - type: 'multi' or 'uni' - multivariate or univariate prewhitening
%                   default multi
% Note: allows multiple distance matrices (e.g. k) prewhitened at the same time
%       - dist (k x d)
%       - V    k cells of size d x d
%       - output k x d
%
type = 'multi';
vararginoptions(varargin,{'type'});
% check the type and dimension of inputs 
if iscell(V)
    nRDM = size(V,2);
else
    nRDM = 1;
    tmp{1} = V;
    clear V; V = tmp;
end
if size(dist,2) ~= size(V{1},1) || (size(V{1},1)~=size(V{1},2))
    error('Incorrect dimensions of input arguments!');
end 

dist_pw=zeros(nRDM,size(V{1},2));
% do the calculation
for i=1:nRDM
    switch type
        case 'multi'
            dist_pw(i,:) = dist(i,:)*V{i}^(-1/2);
        case 'uni'
            dist_pw(i,:) = dist(i,:).*(diag(V{i}).^(-1/2))';
    end
end
