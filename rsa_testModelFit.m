function varargout=rsa_testModelFit(what,varargin)
% Testing for contstructng and Fitting component models
%
%
% Define optional parameter
Opt.rootPath='/Users/joern/Desktop/rsaworkshop/rsatoolbox.v.2a/Demos/Demo_component';
Opt.rootPath='/Users/joern/Talks/2015/02_Demo_component';
% Opt.rootPath='/Users/jdiedrichsen/Talks/2015/02_Demo_component';
tendigitDir='/Users/joern/Talks/2015/07_COSMO/practical';

import rsa.*;
import rsa.util.*;
import rsa.stat.*;
import rsa.rdm.*;

%  Make the representational model matrices from features
switch (what)
    case 'yokoiModel'               % returns Yokoi model of sh1
        load(fullfile(Opt.rootPath,'Yokoi_features.mat'));
        M=rsa.stat.features2ModelRDMs({Features(1).F,'sequence',...
            Features(2).F,'2nd chunk',...
            Features(3).F,'chunk',...
            Features(4).F,'finger'});
        M=rsa.rdm.vectorizeRDMs(M);
        M=rsa_vectorizeIPM(M);
        varargout={M};
    case 'categoryModel'               % returns Yokoi model of sh1
        numStim = 5;    % Stimuli per category
        features{1}=kron(eye(2),ones(numStim,1));
        features{2}=[eye(numStim);zeros(numStim)];
        features{3}=[zeros(numStim);eye(numStim)];
        
        
        M=rsa.stat.features2ModelRDMs({features{1},'between_cat',...
            features{2},'stim1',...
            features{3},'stim2'});
        M=vectorizeRDMs(M);
        M=rsa_vectorizeIPM(M);
        varargout={M};
    case 'randomModel'                  % Random model of rewighting matrix
        numStim = 20;                               % Stimuli in total
        
        F=normrnd(0,1,numStim,numStim-1);
        C = rsa.util.indicatorMatrix('allpairs',[1:numStim]);
        for i=1:numStim-1
            M(i).name=sprintf('Dimension %d',i);
            M(i).IPM = F(:,i)*F(:,i)';
            M(i).RDM = diag(C*M(i).IPM*C')';
        end;
        M=rsa_vectorizeIPM(M);
        varargout={M};
    case 'simulate_data'                % Generate data from model structure: Generate activities and then distances
        % 1. model distance by d = X*omega
        % 2. get G from d by G = (H(-0.5*d)H'), H: centering matrix (I-1/N)
        % 3. get trueU (u=cholcov(G)*normrand(0,1))
        % 4. get Y (= Za*u + noise)
        % 5. get d_hat and sigma_hat by crossvalidation
        %
        % Note that G=UU'= P * sum(omega*Gc);
        % Noise*Noise' = P * var_e
        % This for a SNR of 1:1 omega needs to be same size as var_e
        
        Model  = varargin{1};  % Input Model RDM
        
        % Bring into Model Data frame structure if not already
        if length(Model)>1
            Model = struct2dataframe(Model);
        end;
        
        % Get info of RDMs
        D.numComp    = size(Model.RDM,1);                       % Number of components
        D.numCond    = size(squareform(Model.RDM(1,:)),1);      % Number of conditions
        
        % (Simulated) Experimental parameters: default
        D.numPart       = 9;                % Number of runs
        D.numVox        = 160;              % Number of voxels
        D.omega         = zeros(1,D.numComp);   % hyperparamter on distance model
        D.var_e         = 1;                % Noise variance: Like omega, this is cumulative over voxels
        D.numSim        = 1000;             % Number of simulation
        
        %- allow to get user options
        D               = getUserOptions({varargin{2:end}},D);
        D.N             = D.numPart*D.numCond;          % Number of trials
        part            = kron([1:D.numPart]',ones(D.numCond,1));            % Partitions
        conditions      = kron(ones(D.numPart,1),[1:D.numCond]');            % Conditions
        
        % Get design matrix
        Za  = kron(ones(D.numPart,1),eye(D.numCond));
        
        X           = Model.RDM';
        if (size(D.omega,1)==1)
            D.omega=repmat(D.omega,D.numSim,1);
        end;
        S = [];
        
        for n = 1:D.numSim
            
            % Make true distances
            D.d_true    = D.omega(n,:)*X';
            % Now make true covariance matrix: infer from distances if not
            % given
            if (isfield(Model,'IPM'));   % See if there is a inner product matrix
                G  = D.omega(n,:)*Model.IPM;
                G  = rsa_squareIPM(G);
            else
                H  = eye(D.numCond)-ones(D.numCond)/D.numCond;
                G  = -0.5 * H * squareform(D.d_true) * H';      % mean subtracted
            end;
            
            % Generate true pattern from predicted distances
            trueU = mvnrnd_exact(G,D.numVox)*sqrt(D.numVox);
            % Generate artificial data
            Noise     = sqrt(D.var_e) * randn(D.N,D.numVox);  % Strech the variance of the noise over voxels, just as omega
            Y(:,:,n)  = Za*trueU + Noise;
            
            % Calc cross-validated distance and noise estimate
            [d_hat,Sig_hat] = rsa_distanceLDC(Y(:,:,n),part,conditions);
            
            S.RDM(n,:) = d_hat';
            S.Sig_hat(n,:)= Sig_hat(tril(true(D.numCond),0))';  % Get vectorized version of the variance-covariance matrix
            S.sig_hat(n,1)= mean(diag(Sig_hat));               % Noise variance estimate
        end;
        
        varargout = {S, Y, D};
    case 'simulate_direct'              % Simulates distances directly for speed
        Model  = varargin{1};  % Input Model RDM
        
        % Bring into Model Data frame structure if not already
        if length(Model)>1
            Model = struct2dataframe(Model);
        end;
        
        % Get info of RDMs
        D.numComp    = size(Model.RDM,1);                       % Number of components
        D.numCond    = size(squareform(Model.RDM(1,:)),1);      % Number of conditions
        
        % (Simulated) Experimental parameters: default
        D.numPart       = 9;                % Number of runs
        D.numVox        = 160;              % Number of voxels
        D.omega         = zeros(1,D.numComp);   % hyperparamter on distance model
        D.var_e         = 1;                % Noise variance: Like omega, this is cumulative over voxels
        D.numSim        = 1000;             % Number of simulation
        
        %- allow to get user options
        D               = getUserOptions({varargin{2:end}},D);
        
        % Make true distances
        X           = Model.RDM';
        D.d_true    = D.omega*X';
        
        
        Sigma     = eye(D.numCond)*D.var_e;
        D.numDist = D.numCond*(D.numCond-1)/2;
        C = indicatorMatrix('allpairs',[1:D.numCond]);
        V=rsa.stat.varianceLDC(D.d_true,C,Sigma,D.numPart,D.numVox);
        
        % Generate data
        cholV=cholcov(V);
        
        S.RDM = normrnd(0,1,D.numSim,D.numDist)*cholV;
        S.RDM = bsxfun(@plus,S.RDM,D.d_true);
        
        varargout={S,D};
    case 'test_varianceLDC'             % Tests if our Normal approximation to the Distances is reasonable
        M1=rsa_testModelFit('yokoiModel');
        
        D.omega   = [0 0 0 1];    % component weights
        D.var_e   = 2;    % Noise variance
        D.numPart = 2;
        D.numVox  = 10;
        [S,Y,D]=rsa_testModelFit('simulate_data',M1,D);     % Get simulated data
        
        X=vertcat(M1.RDM);                                                            % Make design matrix
        C = rsa.util.indicatorMatrix('allpairs',[1:D.numCond]);
        dpred = D.omega*X;
        Sigma = eye(D.numCond)*mean(S.sig_hat); % Average, regularized Sigma estimate
        Vpred = rsa_varianceLDC(dpred(1,:),C,Sigma,D.numPart,D.numVox);
        dhat  = mean(S.RDM);
        Vhat  = cov(S.RDM);
        
        scale= max([Vpred(:);Vhat(:)]);
        figure(3);
        subplot(2,2,1);
        imagesc(Vpred,[-scale scale]);
        title('predicted covariance');
        xlabel('Distance'); ylabel('Distance');
        subplot(2,2,2);
        imagesc(Vhat,[-scale scale]);
        title('measured covariance');
        xlabel('Distance'); ylabel('Distance');
        subplot(2,2,3);
        x=[1:length(dpred(1,:))];
        plot(x,dpred(1,:),'k:',x,dhat,'r');
        subplot(2,2,4);
        imagesc(Vhat,[-scale scale]);
        cla;
        colorbar;
        figure(2);
        histplot(S.RDM(:,1:5),'numcat',100);
        
        varargout={S};
    case 'simple_modelfit'
        M(1)=load('Muscle_model.mat');
        M(2)=load('Naturalstats_model.mat');
        D.numSim    = 5000;
        D.var_e   = 1;    % Noise variance
        D.numPart = 8;
        D.numVox  = 80;
        D.omega   = 0.2; % sum(randn(D.numSim,2).^2,2)/10;    % component weights: chi2(2) distributed
        L=[];
        for om=[0.05:0.1:0.8];
             D.omega   = om;
            T=[];
            for m=1:2
                [S,Y,D]=rsa_testModelFit('simulate_data',M(m),D);
                S.truemodel=ones(D.numSim,1)*m;
                % S.omega=D.omega;
                T=addstruct(T,S);
            end;
            T.kendallM1 = corr(T.RDM',M(1).RDM','type','Kendall');
            T.kendallM2 = corr(T.RDM',M(2).RDM','type','Kendall');
            T.spearmanM1 = corr(T.RDM',M(1).RDM','type','Spearman');
            T.spearmanM2 = corr(T.RDM',M(2).RDM','type','Spearman');
            
            T.pearsonM1 = corr(T.RDM',M(1).RDM');
            T.pearsonM2 = corr(T.RDM',M(2).RDM');
            
            nRDM=bsxfun(@rdivide,T.RDM,sqrt(sum(T.RDM.^2,2)));
            T.fixedM1 = nRDM*M(1).RDM';
            T.fixedM2 = nRDM*M(2).RDM';
            
            [T.loglikNull]=rsa.stat.fitModelNull(T.RDM,T.sig_hat,8,D.numVox);
            [T.weightM1,~,T.loglikM1]=...
                rsa.stat.fitModelIRLS(M(1).RDM',T.RDM,T.sig_hat,8,D.numVox);
            [T.weightM2,~,T.loglikM2]=...
                rsa.stat.fitModelIRLS(M(2).RDM',T.RDM,T.sig_hat,8,D.numVox);
            %         subplot(1,2,1);
            %         barplot(T.truemodel,[T.pearsonM1 T.pearsonM2]);
            %         subplot(1,2,2);
            %         barplot(T.truemodel,[T.loglikM1-T.loglikNull T.loglikM2-T.loglikNull]);
            %         s=(T.truemodel-1)*2-1;
            %         pivottable(T.truemodel,[],s.*(T.pearsonM2-T.pearsonM1)>0,'mean');
            %         pivottable(T.truemodel,[],s.*(T.loglikM2-T.loglikM1)>0,'mean');
            s=(T.truemodel-1)*2-1;
            K.propCorr(1,1)=mean(s.*(T.kendallM2-T.kendallM1)>0);
            K.propCorr(2,1)=mean(s.*(T.spearmanM2-T.spearmanM1)>0);
            K.propCorr(3,1)=mean(s.*(T.pearsonM2-T.pearsonM1)>0);
            K.propCorr(4,1)=mean(s.*(T.fixedM2-T.fixedM1)>0);
            K.propCorr(5,1)=mean(s.*(T.loglikM2-T.loglikM1)>0);
            K.omega=ones(5,1)*om;
            K.method=[1:5]';
            L=addstruct(L,K); 
        end;
        lineplot(D.omega,D.propCorr,'split',D.method,'style_thickline','leg',{'spearman','pearson','cosineAng','logLike'},'subset',D.method>1);
        ylabel('Proportion correct');

        varargout={L,T};
    case 'test_bayesRegress'            % Test the Bayesian regression model comparisons
        % Determines the regularisation and sigma from a group of
        % numParticipants measures
        M1 = rsa_testModelFit('chordingModel');
        M2 = M1(1:3);   % lacking the third term
        D.numExp  =   10;           % 100 Experiments
        D.numSubj = 10;     % 12 Partitipants
        D.numSim  = D.numExp * D.numSubj;
        D.omega   = zeros(1,length(M1));
        % D.omega = [0 0 0 0.01];
        D.omega(1)= 0.03;
        D.omega(10)= 0.03;
        D.var_e   = 80;
        D.numPart = 8;
        
        % figure(1);
        % rsa.fig.imageRDMs(M1);
        [S,Y,D]=rsa_testModelFit('simulate_data',M1,D);
        S.exp   = kron([1:D.numExp]',ones(D.numSubj,1));
        S.subj  = kron(ones(D.numExp,1),[1:D.numSubj]');
        for i=1:D.numExp
            indx = find(S.exp==i);
            % Null Model: assuming all distances are zero
            S.logE0(indx,1)=...
                rsa.stat.fitModelNull(S.RDM(indx,:),S.sig_hat(indx,:),D.numPart,D.numVox);
            % GLS fitting: Assuming the noise covariance based on all distances zero
            S.omegaGLS(indx,:) =...
                rsa.stat.fitModelGLS(M1,S.RDM(indx,:),S.sig_hat(indx,:),D.numPart,D.numVox);
            % IRLS fitting: Taking into account the noise covariance under
            % the current model fit
            S.omegaIRLS(indx,:) =...
                rsa.stat.fitModelIRLS(M1,S.RDM(indx,:),S.sig_hat(indx,:),D.numPart,D.numVox);
            % Using overall ridge estimated by empirical Bayes
            [S.omega1(indx,:),S.logE1(indx,1),S.logtheta1(indx,1)]=...
                rsa.stat.fitModelRidgeEB(M1,S.RDM(indx,:),S.sig_hat(indx,:),D.numPart,D.numVox);
            % Using overall ridge and empirical Bayes: Model 2 for
            % comparision
            [S.omega2(indx,:),S.logE2(indx,1),S.logtheta2(indx,1)]=...
                rsa.stat.fitModelRidgeEB(M2,S.RDM(indx,:),S.sig_hat(indx,:),D.numPart,D.numVox);
            % Using Individual ridge paramtersestimated by empirical Bayes
            [S.omegaIn(indx,:),S.logEIn(indx,1),lT,S.logEInSplit(indx,:)]=...
                rsa.stat.fitModelRidgeIndividEB(M1,S.RDM(indx,:),S.sig_hat(indx,:),D.numPart,D.numVox);
            S.logthetaIn(indx,:)=repmat(lT,length(indx),1);
            % Ceiling model
            res=bsxfun(@minus,S.RDM(indx,:),mean(S.RDM(indx,:)));
            S.logEC(indx,1)=...
                rsa.stat.fitModelNull(res,S.sig_hat(indx,:),D.numPart,D.numVox);
        end;
        subplot(3,2,1);
        barplot([],S.omegaGLS)
        subplot(3,2,2);
        barplot([],S.omegaIRLS)
        subplot(3,2,3);
        barplot([],S.omega1);
        subplot(3,2,4);
        barplot([],[S.logE0 S.logE1 S.logE2 S.logEC]);
        subplot(3,2,5);
        barplot([],S.omegaIn);
        subplot(3,2,6);
        barplot([],[S.logEInSplit]);
        
        varargout={S};
end;



