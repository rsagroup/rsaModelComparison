function [omega,logEvidence,theta,logEvidenceSplit]=rsa_fitModelHierarchEB(Model,dist,Sigma,numPart,numVox,varargin);
% function [omega,logEvidence,logtheta]=rsa_fitModelRidgeIndividEB(X,dist,sigma,numPart,numVox,varargin);
% Does a linear fir (component reweighting) or nonlinear fit (remixing)
% fit.
% For linear fit, set the numNonlin field to 0.
% INPUT:
%  Model   : Model structure with optional fields
%    Model.RDM :  numReg  x numDist (x numSubj): Design matrix,
%             if different across subject, 3-d structure with one slice per subject
%    Model.prior: type of Bayesian prior ('Ridge','RidgeIndivid','Zellner','ZellnerIndivid')
%    Model.numComp       : Number of linearly seperable components (at least 1)
%    Model.numPrior      : Number of prior parameters on the component coefficients
%    Model.numNonlin     : Number of nonlinear parameters (>0)
%    Model.nonlinP0      : Starting value of nonlinear(mixing) parameters
%    Model.constantParams: Cell array of additional parameters to function
%    Model.fcn           : Function returning RDM and design matrix the
%                           partial derviates of the design matrix
%    Model.ridgeparam: List of the regularisation parameter for all linear variables
%                      which parameter refers to which column (defaults to [1:numCol])
%
%  dist    : numSubj x numDist          : Data, with subjects stacked in rows
%  Sigma   : numCond x numCond x numSubj : Estimated variance matrix of
%  numPart : Number of partitions
%  numVox  : numSubj x 1 Effective number of voxels
% OPTIONAL
%  'minimizeLength': Length of line search parameter (see help minimize)
% OUTPUT:
%   omega:       estimates of the regression coefficients, 1 row per subject
%   logEvidence: Marginal likelihood of the model
%   logtheta:    logvariance of prior, as determined by empirical Bayes
import rsa.util.*;
import rsa.*;

% Options: Not used right now
Opt = [];
Opt.minimizeLength = 10;
Opt=rsa.getUserOptions(varargin,Opt);


[numSubj,numDist]=size(dist);
numCond = ceil(sqrt(numDist*2));
if (numCond*(numCond-1)/2~=numDist)
    error('bad size of distances');
end;

Y=dist';
C= indicatorMatrix('allpairs',[1:numCond]);

% Make sigma: is this a scalar, or a matrix per subject
for s=1:numSubj
    if (size(Sigma,2)==1)
        S = eye(numCond)*Sigma(s);
    elseif (size(Sigma,3)==numSubj) %
        S = Sigma(:,:,s);
    elseif (size(Sigma,1)==numSubj & size(Sigma,2)==numCond*numCond)
        S = reshape(Sigma(s,:),numCond,numCond);
    elseif  (size(Sigma,1)==numSubj & size(Sigma,2)==numCond*(numCond-1)/2)
        S = rsa_squareIPM(Sigma(s,:));
    else
        error('Sigma does not have legal size.');
    end;
    SigmaDist(:,:,s) = rsa_varianceLDC(zeros(1,numDist),C,S,numPart,numVox(s));
end;

% Simple linear component model
if (Model.numNonlin==0)
    X=permute(Model.RDM,[2 1 3]);
    numReg = size(X,2);
    switch (Model.prior)
        case 'Ridge'
            theta0 = 0;                                               % variance parameters: number of regressors
            theta  = minimize(theta0,@rsa_marglRidgeEB,Opt.minimizeLength,X,Y,SigmaDist);     % Minmize the ridge parameter for the group
            [~,~,omega,~,nlml]=rsa_marglRidgeEB(theta, X, Y, SigmaDist);     % Estimate the regression coefficients seperately
        case 'RidgeIndivid'
            theta0 = zeros(numReg,1);                                               % variance parameters: number of regressors
            theta  = minimize(theta0,@rsa_marglRidgeIndividEB,Opt.minimizeLength,X,Y,SigmaDist);     % Minmize the ridge parameter for the group
            [~,~,omega,~,nlml]=rsa_marglRidgeIndividEB(theta, X, Y, SigmaDist);     % Estimate the regression coefficients seperately
        case 'Zellner'
            theta0 = 0;                                               % variance parameters: number of regressors
            theta  = minimize(theta0,@rsa_marglZellnerEB,Opt.minimizeLength,X,Y,SigmaDist);     % Minmize the ridge parameter for the group
            [~,~,omega,~,nlml]=rsa_marglZellnerEB(theta, X, Y, SigmaDist);     % Estimate the regression coefficients seperately
        case 'ZellnerIndivid'
            theta0 = zeros(numReg,1);                                               % variance parameters: number of regressors
            theta  = minimize(theta0,@rsa_marglZellnerIndividEB,Opt.minimizeLength,X,Y,SigmaDist);     % Minmize the ridge parameter for the group
            [~,~,omega,~,nlml]=rsa_marglZellnerIndividEB(theta, X, Y, SigmaDist);     % Estimate the regression coefficients seperately
        otherwise
            error('Prior needs to be Ridge, RidgeIndivid,Zellner,ZellnerIndivid');
    end;
    % Otherwise nonlinear fit, right now using only a Inidivdual Ridge prior
else
    theta0              = [zeros(Model.numPrior,1);Model.nonlinP0'];                % Add the variance parameters
    [theta, nlmls]      = minimize(theta0,@rsa_marglNonlin,Opt.minimizeLength,Model,Y,SigmaDist);    % Minmize the ridge parameter for the group
    [~,~,omega,~,nlml]  = rsa_marglNonlin(theta, Model, Y, SigmaDist);              % Estimate the regression coefficients seperately
end;
omega       =  omega';
logEvidence =  -nlml';
theta    = theta';

% Also return the logEvidence split by different regressors (remove at some
% point?)
if nargout>3
    X=permute(Model.RDM,[2 1 3]);
    V0 = diag(exp(theta(1:Model.numPrior)));
    [mL1,mL2]=rsa_marglRidgeSplit(V0, X, Y, SigmaDist);                % Estimate the regression coefficients seperately
    logEvidenceSplit=[mL1 mL2];
end;


