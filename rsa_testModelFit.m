function varargout=rsa_testModelFit(what,varargin)
% Testing for Constructing and fitting component models
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
            try
                [d_hat,Sig_hat] = rsa_distanceLDC(Y(:,:,n),part,conditions);
            catch
                tmp1 = indicatorMatrix('identity',conditions);
                tmp2 = indicatorMatrix('allpairs',unique(conditions)');
                [d_hat,Sig_hat] = distance_ldc_sigma(Y(:,:,n),tmp1,tmp2,part);
            end
            S.RDM(n,:) = d_hat';
            S.Sig_hat(n,:)= Sig_hat(tril(true(D.numCond),0))';  % Get vectorized version of the variance-covariance matrix
            S.sig_hat(n,1)= mean(diag(Sig_hat));               % Noise variance estimate
            Sigma(:,:,n) = Sig_hat;
        end;
        S.numVox = ones(size(S.RDM,1),1)*D.numVox;
        varargout = {S, Y, D, Sigma};
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
        % This compares a number of different Ways of applying the prior 
        %  - Single Ridge 
        %  - Individual Ridge 
        %  - Single Zellner 
        %  - Individual Zellner 
        M1 = sh1_getRDMmodelTau1([-1 0 0 0],1,[1 2 4 6],'sqEuclidean');
        M1.X = M1.RDM;      % Design matrix 
        D.numExp  =   80;           % 100 Experiments
        D.numSubj = 10;     % 12 Partitipants
        D.numSim  = D.numExp * D.numSubj;
        % D.omega   = zeros(1,size(M1.RDM,1));
        D.omega = [1 0 0 0]; 
        D.var_e   = 50;
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
                rsa_fitModelNull(S.RDM(indx,:),S.sig_hat(indx,:),D.numPart,S.numVox(indx));
            % IRLS fitting: Taking into account the noise covariance under
            % the current model fit
            S.omegaIRLS(indx,:) =...
                rsa_fitModelIRLS(M1.X',S.RDM(indx,:),S.sig_hat(indx,:),D.numPart,S.numVox(indx));
            
            % Using overall ridge estimated by empirical Bayes
            M1.prior='Ridge'; 
            [S.omega1(indx,:),S.logE1(indx,1),S.logtheta1(indx,1)]=...
                rsa_fitModelHierarchEB(M1,S.RDM(indx,:),S.sig_hat(indx,:),D.numPart,S.numVox(indx));
            
            % Using Individual ridge paramtersestimated by empirical Bayes
            M1.prior='RidgeIndivid'; 
            [S.omega2(indx,:),S.logE2(indx,1),lT]=...
                rsa_fitModelHierarchEB(M1,S.RDM(indx,:),S.sig_hat(indx,:),D.numPart,S.numVox(indx));
            S.logtheta2(indx,:)=repmat(lT,length(indx),1);
            
            % Using overall Zellner g-prior 
            M1.prior='Zellner'; 
            [S.omega3(indx,:),S.logE3(indx,1),S.logtheta3(indx,1)]=...
                rsa_fitModelHierarchEB(M1,S.RDM(indx,:),S.sig_hat(indx,:),D.numPart,S.numVox(indx));
            
            % Using individual Zellner g-prior
            M1.prior='ZellnerIndivid'; 
            [S.omega4(indx,:),S.logE4(indx,1),lT]=...
                rsa_fitModelHierarchEB(M1,S.RDM(indx,:),S.sig_hat(indx,:),D.numPart,S.numVox(indx));
            S.logtheta4(indx,:)=repmat(lT,length(indx),1);

            
            % Ceiling model
            res=bsxfun(@minus,S.RDM(indx,:),mean(S.RDM(indx,:)));
            S.logEC(indx,1)=...
                rsa_fitModelNull(res,S.sig_hat(indx,:),D.numPart,S.numVox(indx));
        end;
        subplot(3,2,1);
        barplot([],S.omegaIRLS)
        title('IRLS'); 
        subplot(3,2,2);
        barplot([],S.omega1);
        title('Ridge'); 
        subplot(3,2,3);
        barplot([],S.omega2);
        title('RidgeIndivid'); 
        subplot(3,2,4);
        barplot([],S.omega3);
        title('Zellner'); 
        subplot(3,2,5);
        barplot([],S.omega4);
        title('ZellnerIndivid'); 
        
        
        varargout={S};
    case 'yokoiModel_synth'               % simulate Yokoi model of chunking experiment         
        % True hyper parameters
        modelTerms  = varargin{1};
        tau         = varargin{2};        % temporal decay of features 
        v           = varargin{3};            % variance of true omegas
        w0          = varargin{4};            % mean of true omegas;        
        
        % Experimental parameters
        D.numExp    = 1;     % 100 Experiments
        D.numSubj   = 14;     % 12 Partitipants
        D.var_e     = 50;
        D.numPart   = 8;        
        D           = getUserOptions({varargin{5:end}},D);
        D.numSim    = D.numExp * D.numSubj;        
        D.omega     = repmat(w0,D.numSim,1); 
        % generate true omega for individual participant: constant 
        % mvnrnd(w0,diag(v),D.numSim);
        % D.omega     = ssqrt(D.omega.*D.omega);
        
        % Define model structure (using one-digit, two-digit, chunk, and
        % sequence models)
        logtau          = log(tau);             % log of tau
        constantParams  = {1,modelTerms,'sqEuclidean'}; %[1 2 4 6]
        Model           = sh1_getRDMmodelTau1(logtau,constantParams{:});
        Model.name_orig = Model.name;
        for m=1:numel(Model.name)
            Model.name{m} = sprintf('%s (logtau=%2.0d)',Model.name{m},logtau(m));
        end
        % figure(1);
        % rsa.fig.imageRDMs(Model);        
        Model.constantParams    = constantParams;
        Model.numComp           = numel(modelTerms);
        Model.numPrior          = numel(v);
        Model.numNonlin         = numel(tau);
        Model.nonlinP0          = zeros(1,Model.numNonlin);
        Model.constantParams    = constantParams;
        Model.fcn               = @sh1_getRDMmodelTau1;        
        
        % Simulate data
        [S,~,D,Sigma] = rsa_testModelFit('simulate_data',Model,D);
        S.exp   = kron([1:D.numExp]',ones(D.numSubj,1));
        S.subj  = kron(ones(D.numExp,1),[1:D.numSubj]');
        
        % adjust true value based on mean of X
        S.omega_true  = D.omega;             
        normX = sqrt(mean(Model.RDM.^2,2));  
        S.omega_trueNorm  = bsxfun(@times,S.omega_true,normX');
        varargout = {S,D,Sigma,Model};                
    case 'test_HierarchEB_checkgrad'            % check if derivatives of logmerginallikelihood are correct
        ptb = varargin{1};
        
        % True hyper parameters
        modelTerms  = [1 2 4 6];
        tau         = [0.05 13 5 1000];        % temporal decay of features 
        v           = [1 1 1 1];            % variance of true omegas
        w0          = [3 0 0 0];            % mean of true omegas;        
                
        % Synthesise data with given hyper parameters
        [S,D,Sigma,Model] = rsa_testModelFit('yokoiModel_synth',modelTerms,tau,v,w0);                
        
        % Check gradients for safety
        numDist = size(S.RDM,2);
        C       = indicatorMatrix('allpairs',[1:D.numCond]);
        for iter=1:10
            % Make sigma
            for s=1:D.numSim
                SigmaDist(:,:,s) = rsa_varianceLDC(zeros(1,numDist),C,Sigma(:,:,s),D.numPart,D.numVox);
            end;
            
            theta0 = 10*rand(size(Theta));
            diff(iter)=checkgrad('rsa_marglNonlin',theta0',ptb,Model,S.RDM',SigmaDist);
            fprintf('theta=[');
            fprintf('%3.0f, ',theta0);
            fprintf(']\n');
        end
    case 'test_HierarchEB_estimate_yokoiModel'  % test model by estimating 
        % True hyper parameters
        modelTerms  = [1 2 4 6];
        tau         = [0.05 13 5 1000];         % temporal decay of features 
        v           = [1 1 1 1];                % variance of true omegas
        w0          = [0 3 0 0];                % mean of true omegas;        

        % Synthesise data with given hyper parameters
        [S,D,Sigma,Model] = rsa_testModelFit('yokoiModel_synth',modelTerms,tau,v,w0);
               
        % Estimate both group hyper parameters using all data and then
        %   estimate omega
        
        % Redefine model for fit
        constantParams          = Model.constantParams;
        reducedModelTerms       = [1 2 4 6]; 
        constantParams{2}       = reducedModelTerms;
        Model                   = sh1_getRDMmodelTau1(zeros(length(reducedModelTerms)),constantParams{:});
        Model.constantParams    = {1,reducedModelTerms,'sqEuclidean'}; %[1 2 4 6]
        Model.numComp           = numel(reducedModelTerms);
        Model.numPrior          = numel(reducedModelTerms);
        Model.numNonlin         = numel(reducedModelTerms);
        Model.nonlinP0          = zeros(1,Model.numNonlin);
        Model.fcn               = @sh1_getRDMmodelTau1; 
        
        for i=1:D.numExp
            indx = find(S.exp==i);
            [S.omega_hat(indx,:),S.logEvidence(indx,:),theta,S.logEvidenceSplit(indx,:)] = ...
                rsa_fitModelHierarchEB(Model,S.RDM(indx,:),Sigma(:,:,indx),D.numPart,repmat(D.numVox,D.numSim,1));
            
            % Summary result
            S.logvarV(indx,:) = kron(ones(D.numSubj,1),theta(1:Model.numPrior));
            S.logtau(indx,:)  = kron(ones(D.numSubj,1),theta(Model.numPrior+1:end));
            
            % post-process of omega to adjust by mean of regressor
            X = sh1_getRDMmodelTau1(theta(Model.numPrior+1:end),Model.constantParams{:});
            normX = sqrt(mean(X.RDM.^2,2));  
            S.omega_hatn(indx,:) = bsxfun(@times,S.omega_hat(indx,:),normX');
        end;
        
        % Plot result: omegas
        figure;
        subplot(2,2,1);
        xpos = myboxplot([],S.logEvidenceSplit(:,1+Model.numComp:end));
        title('log-evidence split');
        set(gca,'XTick',xpos,'XTicklabel',Model.name_orig);
        rotateXLabels(gca,45);
        
        subplot(2,2,2);
        xpos = myboxplot([],S.omega_hatn); hold on
        xrange = xpos(2)-xpos(1);
        drawline(0,'dir','horz');
        for reg=1:length(reducedModelTerms)
            drawline(S.omega_true(1,ismember(modelTerms,reducedModelTerms(reg))),...
                'dir','horz','linewidth',2,'lim',...
                [xpos(reg) xpos(reg)]+[-0.5 0.5]*xrange,'color',[1 0 0])
        end
        title('true and estimated \omega');
        set(gca,'XTick',xpos,'XTicklabel',Model.name_orig);
        rotateXLabels(gca,45);
        
        % Prior parameters
        T = tapply(S,{'exp'},{S.omega_hat,'nanmean','name','omega_hat'},...
            {S.logvarV,'nanmean','name','logvarV'},...
            {S.logtau,'nanmean','name','logtau'},...
            {S.omega_true,'nanmean','name','omega_true'});
        
        subplot(2,2,3);
        barplot([],T.logvarV);
        title('Group-wise hyper parameters log prior var','fontsize',12);
        
        % Nonlinear parameters
        subplot(2,2,4)
        barplot([],T.logtau);
        title('Group-wise hyper parameters log tau','fontsize',12);
        
        varargout={S,T};                
    case 'test_multiple_minima'  % test how to get around the problem of multiple local minima in the nonlinear fit 
        % Simulation 1: Possible confusion between single finger and chunk?
        scenario = varargin{1}; 
        switch (scenario) 
            case 1 
                v           = [0 0 0 0];                % variance of true omegas
                modelTerms  = [1 2 4 6];
                logtau      = [-2 0 2 0];         % temporal decay of features 
                w0          = [1 0 0 0];                % mean of true omegas;        
                nonlinP0    = [-2 0 2 0;...             % Start from true values (1)
                               2 0 -2 0];               % reverse roles          (2)
            case 2 
                v           = [0 0 0 0];                % variance of true omegas
                modelTerms  = [1 2 4 6];
                logtau      = [0 -2 2 0];         % temporal decay of features 
                w0          = [0 3 3 0];                % mean of true omegas;        
                nonlinP0    = [0 -2 2 0;...             % Start from true values 
                               0 2 -2 0];               % reverse roles 
        end; 
                   
        % Synthesise data with given hyper parameters
        D.numExp =10; 
        D.numSubj=10; 
        T=[]; 
        [S,D,Sigma,Model] = rsa_testModelFit('yokoiModel_synth',modelTerms,exp(logtau),v,w0,D);
                       
        % Redefine model for fit
        reducedModelTerms       = [1 2 4 6]; 
        constantParams{2}       = reducedModelTerms;
        Model.numComp           = numel(reducedModelTerms);
        Model.numPrior          = numel(reducedModelTerms);
        Model.numNonlin         = numel(reducedModelTerms);
        Model.fcn               = @sh1_getRDMmodelTau1; 

        for k=1:size(nonlinP0,1) 
            Model.nonlinP0          = nonlinP0(k,:);
            for i=1:D.numExp
                indx = find(S.exp==i);

                [S.omega_hat(indx,:),S.logEvidence(indx,:),theta] = ...
                    rsa_fitModelHierarchEB(Model,S.RDM(indx,:),Sigma(:,:,indx),...
                    D.numPart,repmat(D.numVox,D.numSim,1),'minimizeLength',1000);
            
                % Fitted group parameters 
                S.logvarV(indx,:) = kron(ones(D.numSubj,1),theta(1:Model.numPrior));
                S.logtau(indx,:)  = kron(ones(D.numSubj,1),theta(Model.numPrior+1:end));
            
                % post-process of omega to adjust by norm of regressor
                X = sh1_getRDMmodelTau1(theta(Model.numPrior+1:end),Model.constantParams{:});
                normX = sqrt(mean(X.RDM.^2,2));  
                S.omega_hatNorm(indx,:) = bsxfun(@times,S.omega_hat(indx,:),normX');
            end;
            S.startingVal = ones(size(S.logEvidence))*k; 
            T=addstruct(T,S); 
        end; 
        subplot(2,2,1); 
        myboxplot(T.startingVal,T.omega_hatNorm); 
        title('omega estimates'); 
        xlabel('starting Val'); 
        drawline(0,'dir','horz'); 
        
        for i=1:3 
            subplot(2,2,1+i); 
            scatterplot(T.logtau(T.startingVal==1,i),T.logtau(T.startingVal==2,i),'identity'); 
            title(sprintf('logtau (%d) estimates',i)); 
            xlabel('startingval 1'); 
            ylabel('startingval 2'); 
        end; 
        
        varargout={T};        
end;



