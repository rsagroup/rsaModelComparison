function [V,delta,ksi]=rsa_covarianceDist(partVec,condVec,varargin)
%function [V,delta,ksi]=covariance_dist(partVec,condVec)
%calculates the covariance of distance estimates V, together with the
%signal component (delta) and noise component (ksi)
%
% INPUT:
%       - partVec:  partition vector (K x 1)
%       - condVec:  condition vector (K x 1)
%
% OUTPUT:
%       - V:        covariance of distance estimates
%       - delta:    signal component 
%       - ksi:      error component
%
% VARARGIN:
%       - distType: 'crossval' or 'ncv' (whether distances are
%                   calculated with crossvalidation or not) - changes V
%                   formula
%       - G:        true second moment matrix
%       - sigma:    true noise matrix
%       - data:     beta estimates (K x P; K:cond/run, P:voxels)
%       - nVox:     number of voxels (if data not given, otherwise
%                   estimated from data)
%
distType = 'crossval'; % crossval or ncv
data = [];
G = [];
sigma = [];
nVox = [];
vararginoptions(varargin,{'distType','G','sigma','data','nVox'});

X = indicatorMatrix('identity_p',condVec);          % per condition
C = indicatorMatrix('allpairs',unique(condVec)');   % condition pairs contrast
nPart = numel(unique(partVec));

switch distType % correct denominator in the V equation
case 'crossval'
    den = nPart * (nPart-1);
case 'ncv'
    den = nPart * nPart;
end

% if G not provided - estimate from data
if isempty(G)
    % check if data given
    if isempty(data)
        error('Provide data or G matrix');
    else
        nVox = size(data,2);
        U = pinv(X)*data; % calculate mean pattern across runs
        G = U*U';
    end
end
delta = C*G*C'; % signal component
% if sigma not given, estimate from data
if isempty(sigma)
    sigma=0;
    for i = 1:nPart
        sigma = sigma + (data(partVec==i,:)-delta)*(data(partVec==i,:)-delta)'/(nPart-1)*nVox;
    end
end
ksi = C*sigma*C'; % noise component
V = 4*(delta.*ksi)/nPart + 2*(ksi.^2)/den;
