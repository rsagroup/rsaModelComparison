function varargout=rsa_testModelCompare(what,varargin)
% Testing for Constructing and fitting component models
% Define optional parameter
% pt.rootPath='/Users/joern/Desktop/rsaworkshop/rsatoolbox.v.2a/Demos/Demo_component';
% Opt.rootPath='/Users/joern/Talks/2015/02_Demo_component';
% Opt.rootPath='/Users/jdiedrichsen/Talks/2015/02_Demo_component';
% tendigitDir='/Users/joern/Talks/2015/07_COSMO/practical';
% chordDir='/Users/joern/Projects/ChordPatternDist/analysis';
baseDir = '/Users/jdiedrichsen/Dropbox (Diedrichsenlab)/Projects/modelCompare';

% Use develop branch of the RSA toolbox
import rsa.*;
import rsa.util.*;
import rsa.stat.*;
import rsa.rdm.*;

%  Make the representational model matrices from features
switch (what)
    % Prepare 3 different models and bring them into the correct format
    case 'prep5FingerModel'
        M(1)=load(fullfile(baseDir,'Muscle_model.mat'));
        M(2)=load(fullfile(baseDir,'Naturalstats_model.mat'));
        save(fullfile(baseDir,'Model_fiveFinger.mat'),'M');
        for i=1:length(M)
            M(i).RDM=M(i).RDM./sqrt(sum(M(i).RDM.^2));
        end;
        varargout={M};
    case 'prepChord'
        % Standardisation, but no absolute value
        D=load(fullfile(chordDir,'distance_emg.mat'));
        M(1).RDM=mean(D.dist(D.method==6,:));
        M(1).name='Musclemodel';
        
        D=naturalglove_analyze('naturalstats_distance_direct','Velocity_Cov_MCP_chords_01.mat');
        M(2).RDM=mean(D.dist);
        M(2).name='Naturalstats';
        
        for i=1:length(M)
            M(i).RDM=M(i).RDM./sqrt(sum(M(i).RDM.^2));
        end;
        save(fullfile(baseDir,'Model_chords.mat'),'M');
        varargout={M};
    case 'prepStartA'
        load(fullfile(baseDir,'START_B_visualiseRDMs','modelRDMs_A2.mat'));
        M=modelRDMs_struct;
        indx=[];
        for i=1:length(M)
            if ~isempty(strfind(M(i).name,'compl'))
                indx=[indx i];
            end;
        end;
        M=M(indx(1:end-1));
        for i=1:length(M)
            M(i).RDM=M(i).RDM./sqrt(sum(M(i).RDM.^2));
        end;
        save(fullfile(baseDir,'START_compl.mat'),'M');
    case 'evaluateSubspace' % Looks at the subspace overlap of models in terms of their principal components
        M=varargin{1};
        
        for m=1:2
            % Prep variance components
            RDM{m}=rsa.rdm.squareRDM(M{m}.RDM);
            numCat=size(RDM{m},1);
            H = eye(numCat)-ones(numCat)/numCat;  % Centering matrix
            G{m}=-0.5*H*RDM{m}*H;
            % Get the design matrix for encoding models
            [V,Lam]=eig(G{m});
            [lambda{m},i]=sort(diag(Lam),'descend');
            X{m}=bsxfun(@times,V(:,i),sqrt(lambda{m}'));
        end;
        for k=1:numCat-1
            R2(k)=subspace_overlap(X{1}(:,1:k),X{2}(:,1:k));
        end
        subplot(2,1,1); % Plot the amount of variance captured
        cs1 =cumsum(lambda{1})./sum(lambda{1});
        cs2= cumsum(lambda{2})./sum(lambda{2});
        plot([1:numCat]',cs1,'b',[1:numCat]',cs2,'r');
        subplot(2,1,2);
        plot([1:k],R2);
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
        
        X           = double(Model.RDM');
        % Make true distances
        D.d_true    = double(D.omega*X');
        
        % Now make true covariance matrix: infer from distances if not
        % given
        if (isfield(Model,'IPM'));   % See if there is a inner product matrix
            G  = D.omega*Model.IPM;
            G  = rsa_squareIPM(G);
        else
            H  = eye(D.numCond)-ones(D.numCond)/D.numCond;
            G  = -0.5 * H * squareform(D.d_true) * H';      % mean subtracted
        end;

        
        S = [];
        
        for n = 1:D.numSim
            
            
            % Scale if necessary
            if (isfield(D,'phi'))
                G=G.*D.phi(n);
            end;
            
            % Generate true pattern from predicted distances
            trueU = mvnrnd_exact(G,D.numVox);
            % Generate artificial data
            Noise     = sqrt(D.var_e) * randn(D.N,D.numVox);  % Strech the variance of the noise over voxels, just as omega
            Y(:,:,n)  = Za*trueU + Noise;
            
            % Calc cross-validated distance and noise estimate
            % try
            [d_hat,Sig_hat] = rsa.distanceLDC(Y(:,:,n),part,conditions);
            Xcon = indicatorMatrix('identity',conditions);
            d_1 = pdist(pinv(Xcon)*Y(:,:,n),'euclidean');
            % catch
            %
            %     tmp2 = indicatorMatrix('allpairs',unique(conditions)');
            %     [d_hat,Sig_hat] = distance_ldc_sigma(Y(:,:,n),tmp1,tmp2,part);
            % end
            S.RDM(n,:) = d_hat';
            S.RDMn(n,:)= (d_1.*d_1)';
            S.Sig_hat(n,:)= Sig_hat(tril(true(D.numCond),0))';  % Get vectorized version of the variance-covariance matrix
            S.sig_hat(n,1)= mean(diag(Sig_hat));               % Noise variance estimate
            Sigma(:,:,n) = Sig_hat;
        end;
        S.numVox = ones(size(S.RDM,1),1)*D.numVox;
        varargout = {S, Y, D, Sigma,conditions,part};
    case 'modelCompare'
        % RSA_methods={'spearman','pearson','fixed','loglikIRLS'};
        % RSA simulations
        % rsa_testModelCompare('modelCompare','model','Model_fiveFinger.mat','numSim',1000,'outfile','sim_rsa_Exp1.mat','methods',RSA_methods,'Omega',[0:0.1:0.8]);
        % rsa_testModelCompare('modelCompare','model','Model_chords.mat','numSim',100,'outfile','sim_rsa_Exp2.mat','methods',RSA_methods,'Omega',[0:0.05:0.3]);
        % rsa_testModelCompare('modelCompare','model','Model_START.mat','numSim',1,'outfile','sim_rsa_Exp3.mat','methods',RSA_methods,'Omega',[0:0.1:0.8]);
        % NEW RSA simulation
        % RSA_methods={'spearman','pearson','pearsonNc','pearsonSq','pearsonNcSq','cosine','loglikPCM'};
        % rsa_testModelCompare('modelCompare','model','Model_fiveFinger.mat','numSim',1000,'outfile','sim_rsan_Exp1.mat','methods',RSA_methods,'Omega',[0:0.1:0.8]);
        % rsa_testModelCompare('modelCompare','model','Model_chords.mat','numSim',1000,'outfile','sim_rsan_Exp2a.mat','methods',RSA_methods,'Omega',[0:0.05:0.3]);
        % rsa_testModelCompare('modelCompare','model','Model_START.mat','numSim',200,'outfile','sim_rsan_Exp3a.mat','methods',RSA_methods,'Omega',[0:0.1:0.8]);
        % rsa_testModelCompare('modelCompare','model','Model_START.mat','numSim',200,'outfile','sim_rsan_Exp3b.mat','methods',RSA_methods,'Omega',[0:0.1:0.8]);   
        % 
        % Probabilistic RSA simlation 
        % RSA_methods={'pearson','cosine','cosineWNull','cosineWData','loglikIRLS'};
        % rsa_testModelCompare('modelCompare','model','Model_fiveFinger.mat','numSim',1000,'outfile','sim_rsaProb_Exp1.mat','methods',RSA_methods,'Omega',[0:0.1:0.8]);
        % rsa_testModelCompare('modelCompare','model','Model_chords.mat','numSim',1000,'outfile','sim_rsaProb_Exp2a.mat','methods',RSA_methods,'Omega',[0:0.05:0.3]);
        % rsa_testModelCompare('modelCompare','model','Model_START.mat','numSim',50,'outfile','sim_rsaProb_Exp3a.mat','methods',RSA_methods,'Omega',[0:0.1:0.6]);
        % 
        % 
        % PCM simulations
        %         rsa_testModelCompare('modelCompare','model','Model_fiveFinger.mat','numSim',500,'outfile','sim_pcm_Exp1.mat','methods',{'loglikPCM'},'Omega',[0:0.1:0.8]);
        %         rsa_testModelCompare('modelCompare','model','Model_chords.mat','numSim',100,'outfile','sim_pcm_Exp2.mat','methods',{'loglikPCM'},'Omega',[0:0.05:0.3]);
        %         rsa_testModelCompare('modelCompare','model','Model_START.mat','numSim',100,'outfile','sim_pcm_Exp3.mat','methods',{'loglikPCM'},'Omega',[0:0.1:0.8]);
        %
        % OPT_methods={'encodePCM','loglikIRLS','loglikPCM'};
        % rsa_testModelCompare('modelCompare','model','Model_fiveFinger.mat','numSim',1000,'outfile','sim_opt_Exp1a.mat','methods',OPT_methods,'Omega',0.3);
        % rsa_testModelCompare('modelCompare','model','Model_chords.mat','numSim',1000,'outfile','sim_opt_Exp2a.mat','methods',OPT_methods,'Omega',0.15);
        % rsa_testModelCompare('modelCompare','model','Model_START.mat','numSim',500,'outfile','sim_opt_Exp3a.mat','methods',OPT_methods,'Omega',0.5);
        %
        % Encode simulations
        % rsa_testModelCompare('modelCompare','model','Model_fiveFinger.mat','numSim',1000,'outfile','sim_encPCMCorr_Exp1.mat','methods',{'encodePCMCorr'},'Omega',[0:0.1:0.8]);
        % rsa_testModelCompare('modelCompare','model','Model_fiveFinger.mat','numSim',1000,'outfile','sim_encRidgeCorr_Exp1.mat','methods',{'encodeRidgeCorr'},'Omega',[0:0.1:0.8]);
        % rsa_testModelCompare('modelCompare','model','Model_chords.mat','numSim',1000,'outfile','sim_encRidgeCorr_Exp2.mat','methods',{'encodeRidgeCorr'},'Omega',[0:0.05:0.3]);
        % rsa_testModelCompare('modelCompare','model','START_compl.mat','numSim',500,'outfile','sim_encRidgeCorr_Exp3.mat','methods',{'encodeRidgeCorr'},'Omega',[0:0.1:0.8]);
        Opt.methods ={'kendall','spearman','pearson','fixed','loglikIRLS','loglikPCM','encodeReg','encodePCM'};
        Opt.methods ={'encodePCMcorr'};
        Opt.model   = 'Model_fiveFinger.mat';
        Opt.numSim  = 1000;
        Opt.numPart = 8;
        Opt.numVox  = 160;
        Opt.outfile = [];
        Opt.Omega   = [0 0.3 0.6 0.9];
        
        Opt=rsa.getUserOptions(varargin,Opt);
        
        L=[];   % Summary stats
        U=[];   % details on each simulation
        
        numMethods = length(Opt.methods);
        if (ischar(Opt.model))
            load(fullfile(baseDir,Opt.model));
        else
            M=Opt.model;
        end;
        
        numModels = length(M);
        
        % Experimental constants
        test=rsa.rdm.squareRDM(M{1}.RDM);
        D.numCat = size(test,1);
        D.var_e   = 1;    % Noise variance
        D.numPart = Opt.numPart;
        D.numVox  = Opt.numVox;
        D.numSim  = Opt.numSim;
        
        % Prep some convenience matrices
        H = eye(D.numCat)-ones(D.numCat)/D.numCat;  % Centering matrix
        part  = kron([1:D.numPart]',ones(D.numCat,1));
        cond  = kron(ones(D.numPart,1),[1:D.numCat]');
        Xpart = indicatorMatrix('identity',part);
        Xcond = indicatorMatrix('identity',cond);                            
        C=pcm_indicatorMatrix('allpairs',[1:D.numCat]);  % Contrast vector for all pairs 

        
        % Prep the Variance components and regression matricies for each model
        for m=1:numModels
            % Prep variance components
            RDM{m}=rsa.rdm.squareRDM(M{m}.RDM);
            G{m}=-0.5*H*RDM{m}*H;
            % Get the design matrix for encoding models
            [V,Lam]=eig(G{m});
            [lambda{m},i]=sort(diag(Lam),'descend');
            X{m}=bsxfun(@times,V(:,i),sqrt(lambda{m}'));
            X{m}=X{m}(:,lambda{m}>eps);
            numReg(m) = size(X{m},2);
        end;
        
        % Loop over different signal levels
        for om=Opt.Omega;
            D.omega   = om;
            fprintf('%2.3f\n',om);
            T=[];
            % Generate data from the two models
            Y=[];
            for m=1:numModels
                [S,y,D]=rsa_testModelCompare('simulate_data',M{m},D);
                S.truemodel=ones(D.numSim,1)*m;
                T=addstruct(T,S);
                Y=cat(3,Y,y);
            end;
            T.omega = ones(D.numSim*numModels,1)*om;
            for m=1:numModels
                % Calculate the different forms of correlation,
                for meth = 1:numMethods
                    tic;
                    switch(Opt.methods{meth})
                        case 'kendall'
                            T.kendall(:,m)    = corr(T.RDMn',M{m}.RDM','type','Kendall');
                        case 'spearman'
                            T.spearman(:,m)   = corr(T.RDMn',M{m}.RDM','type','Spearman');
                        case 'pearson'
                            T.pearson(:,m)    = corr(T.RDM',M{m}.RDM');
                        case 'pearsonNc'
                            T.pearsonNc(:,m)    = corr(T.RDMn',M{m}.RDM');
                        case 'pearsonNcSq'
                            T.pearsonNcSq(:,m)    = corr(sqrt(T.RDMn'),sqrt(M{m}.RDM'));
                        case 'pearsonSq'
                            T.pearsonSq(:,m)    = corr(sqrt(T.RDM'),sqrt(M{m}.RDM'));
                        case 'cosine'
                            nRDM=bsxfun(@rdivide,T.RDM,sqrt(sum(T.RDM.^2,2)));
                            mRDM=(M{m}.RDM')/sqrt(sum(M{m}.RDM.^2));
                            T.cosine(:,m)   =   nRDM*mRDM;
                        case 'cosineWNull' % Cosine angle weighted by the covariance structure under the Null-hypothesis
                            varD = rsa_varianceLDC(zeros(D.numCat),C,mean(T.sig_hat),D.numPart,D.numVox); % Get the variance 
                            T.cosineWNull(:,m)   =   cosineW(T.RDM,M{m}.RDM,varD); 
                        case 'cosineWData' % Cosine angle weighted by the covariance structure under the Full hypothesis
                            for n=1:size(T.RDM,1)
                                G=H*squareform(-0.5*T.RDM(n,:))*H'; 
                                G=pcm_makePD(G); 
                                varD = rsa_varianceLDC(G,C,T.sig_hat(n,:),D.numPart,D.numVox);
                                T.cosineWData(n,m)   =  cosineW(T.RDM(n,:),M{m}.RDM,varD);
                            end; 
                        case 'fixed'
                            nRDM=bsxfun(@rdivide,T.RDM,sqrt(sum(T.RDM.^2,2)));
                            [T.weight(:,m),T.fixed(:,m),T.loglikeFixed(:,m)] = ...
                                rsa.stat.fitModelOLS(M{m}.RDM',T.RDM);
                        case 'loglikIRLS'
                            % Likelihood under normal approximation with assumed Sigma
                            [T.weight(:,m),~,T.loglikIRLS(:,m)]=...
                                rsa_fitModelIRLS(M{m}.RDM',T.RDM,T.sig_hat,D.numPart,D.numVox);
                        case 'loglikIRLSsig'
                            % Likelihood under normal approximation with
                            % inferred structure, but open sigma for
                            % scaling
                            [T.weight(:,m),~,T.loglikIRLSsig(:,m)]=...
                                rsa_fitModelIRLS(M{m}.RDM',T.RDM,T.sig_hat,D.numPart,D.numVox,'assumeVoxelInd',0);
                        case 'loglikPCM'
                            % PCM-based likelihood.
                            OPT.fitScale = 1;
                            OPT.runEffect = [];
                            for n=1:size(T.RDM,1)
                                YY = Y(:,:,n)*Y(:,:,n)';
                                fcn = @(x) pcm_likelihoodIndivid(x,YY,M{m},Xcond,Xpart,Opt.numVox,OPT);
                                [theta,T.loglikPCM(n,m)] = pcm_NR(zeros(2,1),fcn);
                                T.PCMs(n,m)=theta(1);
                                T.PCMe(n,m)=theta(2);
                            end;
                        case 'encodeReg'
                            % Encoding model without regularisation
                            for n=1:size(T.RDM,1)
                                [T.encodeReg(n,m),~,~,T.encodeRegCorr(n,m)] = encode_crossval(Y(:,:,n),Xcond*X{m}(:,1:2),part,'linregress','X',Xpart);
                            end;
                        case 'encodePCM'
                            % Encoding model with regularisation
                            for n=1:size(T.RDM,1)
                                [T.encodePCM(n,m)] = encode_crossval(Y(:,:,n),Xcond*X{m},part,'pcm_NR_diag','X',Xpart,'Gd',ones(numReg(m),1));
                            end;
                        case 'encodePCMcorr'
                            % Encoding model with regularisation
                            for n=1:size(T.RDM,1)
                                [~,T.encodePCMcorr(n,m)] = encode_crossval(Y(:,:,n),Xcond*X{m},part,'pcm_NR_diag','X',Xpart,'Gd',ones(numReg(m),1));
                            end;
                        case 'encodeRidge'
                            % Encoding model with regularisation
                            for n=1:size(T.RDM,1)
                                [T.encodeRidge(n,m)] = encode_crossval(Y(:,:,n),Xcond*X{m},part,'ridgeFixed','X',Xpart,'G',eye(numReg(m))*om,'sigma2',1);
                            end;
                        case 'encodeRidgeCorr'
                            % Encoding model with regularisation
                            for n=1:size(T.RDM,1)
                                [~,T.encodeRidgeCorr(n,m)] = encode_crossval(Y(:,:,n),Xcond*X{m},part,'ridgeFixed','X',Xpart,'G',eye(numReg(m))*om,'sigma2',1);
                            end;
                    end;
                    exetime(m,meth) = toc;
                end;
                fprintf('model %d done\n',m);
            end;
            % Now calculate the proportion of cases the correct model wins
            N = size(T.truemodel,1);
            for m=1:numMethods
                meancorrect = nan(N,1);
                for n=1:N
                    val = T.(Opt.methods{m})(n,:);
                    trueval = val(T.truemodel(n));
                    val(T.truemodel(n))=[];
                    meancorrect(n,:) = (sum(trueval>val) + ...
                        sum(trueval==val)*0.5)./(numModels-1);
                end;
                K.propCorr(m,1)  = mean(meancorrect);
                K.method(m,1)    = m;
                K.methodStr{m,1} = Opt.methods{m};
                K.avrgTime(m,1)    = mean(exetime(:,m));
            end;
            K.omega=ones(numMethods,1)*om;
            L=addstruct(L,K);
            U=addstruct(U,T);
        end;
        
        % Plot the results
        lineplot(L.omega,L.propCorr,'split',L.method,'style_thickline','leg',Opt.methods);
        ylabel('Proportion correct');
        drawline(0.5,'dir','horz');
        U.methods = Opt.methods;
        T=L;
        if (~isempty(Opt.outfile))
            save(fullfile(baseDir,Opt.outfile),'T','U');
        end;
        varargout={U,T};
    case 'calc_mean_correct'
        criterion = varargin{1};
        truemodel = varargin{2};
        [N,numModels]=size(criterion);
        meancorrect = nan(N,1);
        for n=1:N
            val = criterion(n,:);
            trueval = val (truemodel(n));
            val(truemodel(n))=[];
            meancorrect(n,:) = (sum(trueval>val) + ...
                sum(trueval==val)*0.5)./(numModels-1);
        end;
        varargout = {mean(meancorrect),meancorrect};
    case 'simple_modelfit_plot'
        doplot=[1 3 5 6 7 8];
        linecolor={'b','b','b','g','r','r'};
        linestyle={':','--','-','-',':','-'};
        D=load('simple_model.mat');
        methods = D.methods;
        D=rmfield(D,'methods');
        numMethods = length(methods);
        L=[];
        for om=unique(D.omega)'
            T=getrow(D,D.omega==om);
            s=(T.truemodel-1)*2-1;
            
            for m=1:length(methods)
                K.propCorr(m,1) = (sum((s.*diff(T.(methods{m}),[],2))>0)+...
                    sum(diff(T.(methods{m}),[],2)==0)*0.5) ...
                    ./size(T.omega,1);
                
            end;
            K.omega=ones(numMethods,1)*om;
            K.method=[1:numMethods]';
            L=addstruct(L,K);
        end;
        
        % Plot the results
        lineplot(L.omega,L.propCorr,'split',L.method,'style_thickline','leg',{methods{doplot}},...
            'linecolor',linecolor,'linestyle',linestyle,...
            'markercolor',linecolor,'markerfill',linecolor,...
            'markersize',5,'errorbars','none','subset',ismember(L.method,doplot));
        ylabel('Proportion correct');
        drawline(0.5,'dir','horz');
        varargout={L,D,methods};
    case 'test_scalefit'
        % This tests for different ways of doing a combined model fit
        % across a number of subjects, each having a different scale
        % of distances
        M1 = sh1_getRDMmodelTau1([-2 1.5 -0.5 0],1,[1 2 4 6],'sqEuclidean');
        M1.X = M1.RDM;      % Design matrix
        D.numExp  =  50;           % 100 Experiments
        D.numSubj = 14;     % Partitipants / experiment
        D.numSim  = D.numExp * D.numSubj;
        D.omega = [0 0 0.6 0];
        D.var_e   = 20;
        D.numPart = 8;
        % D.phi = chi2rnd(10,D.numSim,1)/10;
        models={[1 2 3 4],[2 3 4],[1 3 4],[1 2 3]};
        % figure(1);
        % rsa.fig.imageRDMs(M1);
        corrN(M1.X')
        [S,Y,D]=rsa_testModelFit('simulate_data',M1,D);
        
        S.exp   = kron([1:D.numExp]',ones(D.numSubj,1));
        S.subj  = kron(ones(D.numExp,1),[1:D.numSubj]');
        for i=1:D.numExp
            indx = find(S.exp==i);
            fprintf('%d\n',i);
            % for m=1:4
            %             S.R2(indx,1)=rsa_crossvalModel(getrow(M1,models{m}),S.RDM(indx,:),'method','OLS');
            S.R2(indx,1)=rsa_crossvalModel(M1,S.RDM(indx,:),'method','OLS');
            S.R2(indx,2)=rsa_crossvalModel(M1,S.RDM(indx,:),'method','NLS');
            % S.R2(indx,3)=rsa_crossvalModel(M1,S.RDM(indx,:),'method','ScaleOLS');
            % end;
        end;
        % S.R2d = bsxfun(@minus,S.R2(:,2:4),S.R2(:,1));
        barplot([],S.R2);
        varargout={S};
    case 'test_pcm_fixed' % Test PCM with a fixed block effect
        Opt.model   = 'Model_fiveFinger.mat';
        Opt.numSim  = 10;
        Opt.numPart = 3;
        Opt.numVox  = 40;
        Opt.outfile = [];
        Opt.Omega   = 0.2;
        U=[];
        
        Opt=rsa.getUserOptions(varargin,Opt);
        
        if (ischar(Opt.model))
            load(fullfile(baseDir,Opt.model))
        else
            M=Opt.model;
        end;
        numModels = length(M);
        
        % Experimental constants
        test=rsa.rdm.squareRDM(M(1).RDM);
        D.numCat = size(test,1);
        D.var_e   = 1;    % Noise variance
        D.numPart = Opt.numPart;
        D.numVox  = Opt.numVox;
        D.numSim  = Opt.numSim;
        
        % Prep some convenience matrices
        H = eye(D.numCat)-ones(D.numCat)/D.numCat;  % Centering matrix
        part  = kron([1:D.numPart]',ones(D.numCat,1));
        cond  = kron(ones(D.numPart,1),[1:D.numCat]');
        Xpart = indicatorMatrix('identity',part);
        Xcond = indicatorMatrix('identity',cond);
        
        % Prep the Variance components and regression matricies for each model
        for m=1:numModels
            % Prep variance components
            if isfield(M{m},'IPM')
                G{m}=rsa_squareIPM(M{m}.IPM);
            else
                RDM{m}=rsa.rdm.squareRDM(M{m}.RDM);
                G{m}=-0.5*H*RDM{m}*H;
            end;
            % Get the design matrix for encoding models
            [V,Lam]=eig(G{m});
            [lambda{m},i]=sort(diag(Lam),'descend');
            X{m}=bsxfun(@times,V(:,i),sqrt(lambda{m}'));
            X{m}=X{m}(:,lambda{m}>eps);
            numReg(m) = size(X{m},2);
        end;
        
        % Loop over different signal levels
        for om=Opt.Omega;
            D.omega   = om;
            fprintf('%2.3f\n',om);
            T=[];
            % Generate data from the two models
            Y=[];
            for m=1:numModels
                [S,y,D]=rsa_testModelCompare('simulate_data',M{m},D);
                S.truemodel=ones(D.numSim,1)*m;
                T=addstruct(T,S);
                Y=cat(3,Y,y);
            end;
            
            T.omega = ones(D.numSim*numModels,1)*om;
            for m=1:numModels
                for n=1:size(Y,3)
                    [G1,theta,~,l]=pcm_NR_diag(Y(:,:,n),Xcond*X{m},'Gd',ones(numReg(m),1),'X',Xpart,'hP',0);
                    T.loglik_NRdiag(n,m) = l(1);
                    T.signal_NRdiag(n,m)=theta(1);
                    T.noise_NRdiag(n,m)=theta(2);
                    [G3,theta2,~,l]=pcm_minimize(Y(:,:,n),Xcond,'Gc',{G{m}},'X',Xpart,'MaxIteration',1000);
                    T.loglik_min(n,m) = l(1);
                    T.signal_min(n,m)=theta2(1);
                    T.noise_min(n,m)=theta2(2);
                    [G2,theta1,~,l]=pcm_NR(Y(:,:,n),Xcond,'Gc',{G{m}},'X',Xpart,'hP',0);
                    T.loglik_NR(n,m) = l(1);
                    T.signal_NR(n,m)=theta1(1);
                    T.noise_NR(n,m)=theta1(2);
                end;
                fprintf('model %d done\n',m);
            end;
            % Now calculate the proportion of cases the correct model wins
            U=addstruct(U,T);
        end;
        varargout={U};
    case 'test_pcm'                     % Test PCM for comparisons
        M1 = sh1_getRDMmodelTau1([-2 1.5 -1 0],1,[1 2 4 6],'sqEuclidean');
        M1.X = M1.RDM;      % Design matrix
        D.numExp  =   1;           % 100 Experiments
        D.numSubj = 5;     % 12 Partitipants
        D.numSim  = D.numExp * D.numSubj;
        % D.omega   = zeros(1,size(M1.RDM,1));
        D.omega = [0 0 0.6 0];
        D.var_e   = 80;
        D.numPart = 8;
        ModelComp={1,2,3,4,[1 2 3 4],[2 3 4],[1 3 4],[1 2 4],[1 2 3]};
        
        % figure(1);
        % rsa.fig.imageRDMs(M1);
        corrN(M1.X');
        [S,Y,D,~,cond,part]=rsa_testModelFit('simulate_data',M1,D);
        
        S.exp   = kron([1:D.numExp]',ones(D.numSubj,1));
        S.subj  = kron(ones(D.numExp,1),[1:D.numSubj]');
        for i=1:D.numExp.*D.numSubj
            for m=1:length(ModelComp)
                M2=getrow(M1,ModelComp{m});
                [S.logEv(i,m)]=rsa_fitModelPCM(M2,Y(:,:,i),cond,part);
            end;
        end;
        keyboard;
        varargout={S};
    case 'encodeReg'
        % Test different dimensionalties of encoding models
        % rsa_testModelCompare('encodeReg','model','Model_fiveFinger.mat','numSim',1000,'outfile','sim_encodeReg_Exp1.mat');
        % rsa_testModelCompare('encodeReg','model','Model_fiveFinger.mat','numSim',50,'method','encodePCM','outfile','sim_encodePCM_Exp1.mat');
        % rsa_testModelCompare('encodeReg','model','Model_chords.mat','numSim',100,'outfile','sim_encodeReg_Exp2.mat','numReg',[1:8 10:2:20 25 30],'Omega',0.1);
        % rsa_testModelCompare('encodeReg','model','Model_chords.mat','numSim',10,'method','encodePCM','outfile','sim_encodePCM_Exp2.mat','numReg',[1:8 10:2:20 25 30],'Omega',0.1);
        % rsa_testModelCompare('encodeReg','model','Model_chords.mat','numSim',10,'method','encodeRidge','outfile','sim_encodeRidge_Exp2.mat','numReg',[1:8 10:2:20 25 30],'Omega',0.1);
        % rsa_testModelCompare('encodeReg','model','START_compl.mat','numSim',10,'outfile','sim_encodeReg_Exp3.mat','numReg',[1 5:5:30 40:10:90]);
        % rsa_testModelCompare('encodeReg','model','START_compl.mat','numSim',10,'method','encodePCM','outfile','sim_encodePCM_Exp3.mat','numReg',[1 5:5:30 40:10:90]);
        
        Opt.model   = 'Model_fiveFinger.mat';
        % Opt.outfile = 'sim_encodeReg_Exp1.mat';
        Opt.outfile = [];
        Opt.numSim  = 100;
        Opt.numPart = 8;
        Opt.numVox  = 160;
        Opt.Omega   = 0.3;
        Opt.numReg  = [];  % Number of regressors tested
        Opt.method  = 'encodeReg';
        
        Opt=rsa.getUserOptions(varargin,Opt);
        
        L=[];   % Summary stats
        U=[];   % details on each simulation
        
        if (ischar(Opt.model))
            load(fullfile(baseDir,Opt.model))
        else
            M=Opt.model;
        end;
        
        numModels = length(M);
        
        % Experimental constants
        test=rsa.rdm.squareRDM(M(1).RDM);
        D.numCat = size(test,1);
        D.var_e   = 1;    % Noise variance
        D.numPart = Opt.numPart;
        D.numVox  = Opt.numVox;
        D.numSim  = Opt.numSim;
        
        % Prep some convenience matrices
        H = eye(D.numCat)-ones(D.numCat)/D.numCat;  % Centering matrix
        part  = kron([1:D.numPart]',ones(D.numCat,1));
        cond  = kron(ones(D.numPart,1),[1:D.numCat]');
        Xpart = indicatorMatrix('identity',part);
        Xcond = indicatorMatrix('identity',cond);
        
        % Prep the Variance components and regression matricies for each model
        for m=1:numModels
            % Prep variance components
            RDM{m}=rsa.rdm.squareRDM(M{m}.RDM);
            G{m}=-0.5*H*RDM{m}*H;
            % Get the design matrix for encoding models
            [V,Lam]=eig(G{m});
            [lambda{m},i]=sort(diag(Lam),'descend');
            X{m}=bsxfun(@times,V(:,i),sqrt(lambda{m}'));
            X{m}=X{m}(:,lambda{m}>eps);
        end;
        if (isempty(Opt.numReg))
            Opt.numReg = [1:size(X{1},2)];   % Number of regressors
        end;
        % Loop over different signal levels
        for om=Opt.Omega;
            D.omega   = om;
            fprintf('%2.3f\n',om);
            T=[];
            % Generate data from the two models
            Y=[];
            for m=1:numModels
                [S,y,D]=rsa_testModelFit('simulate_data',M{m},D);
                S.truemodel=ones(D.numSim,1)*m;
                T=addstruct(T,S);
                Y=cat(3,Y,y);
            end;
            T.omega = ones(D.numSim*numModels,1)*om;
            for m=1:numModels
                % Calculate the different forms of correlation,
                % Encoding model without regularisation
                for n=1:size(T.RDM,1)
                    if (mod(n,ceil(Opt.numSim/20))==0);
                        fprintf('.');
                    end;
                    for v=Opt.numReg
                        switch(Opt.method)
                            case 'encodeReg'
                                [T.(sprintf('R2_%2.2d',v))(n,m),~,~,T.(sprintf('corr_%2.2d',v))(n,m)]= ...
                                    encode_crossval(Y(:,:,n),Xcond*X{m}(:,1:v),part,'linregress','X',Xpart);
                            case 'encodePCM'
                                [T.(sprintf('R2_%2.2d',v))(n,m),~,~,T.(sprintf('corr_%2.2d',v))(n,m)]= ...
                                    encode_crossval(Y(:,:,n),Xcond*X{m}(:,1:v),part,'pcm_NR_diag','X',Xpart,'Gd',ones(v,1));
                            case 'encodeRidge'
                                [T.(sprintf('R2_%2.2d',v))(n,m),~,~,T.(sprintf('corr_%2.2d',v))(n,m)]= ...
                                    encode_crossval(Y(:,:,n),Xcond*X{m}(:,1:v),part,'ridgeFixed','X',Xpart,'G',eye(v)*om,'sigma2',D.var_e);
                        end;
                    end;
                end;
                fprintf('model %d done\n',m);
            end;
            % Now calculate the proportion of cases the correct model wins
            N = size(T.truemodel,1);
            for i=1:length(Opt.numReg)
                name1=sprintf('R2_%2.2d',Opt.numReg(i));
                name2=sprintf('corr_%2.2d',Opt.numReg(i));
                for n=1:N
                    val1 = T.(name1)(n,:);
                    trueval1 = val1(T.truemodel(n));
                    val1(T.truemodel(n))=[];
                    meancorrect1(n,:) = (sum(trueval1>val1) + ...
                        sum(trueval1==val1)*0.5)./(numModels-1);
                    
                    val2 = T.(name2)(n,:);
                    trueval2 = val2(T.truemodel(n));
                    val2(T.truemodel(n))=[];
                    meancorrect2(n,:) = (sum(trueval2>val2) + ...
                        sum(trueval2==val2)*0.5)./(numModels-1);
                    trueVal(n,:) = [trueval1 trueval2];
                end;
                K.numReg(i,1)         = Opt.numReg(i);
                K.propCorrR2(i,1)     = mean(meancorrect1);
                K.propCorrCorr(i,1)   = mean(meancorrect2);
                K.meanR2(i,:)         = mean(trueVal(:,1));
                K.meanCorr(i,:)       = mean(trueVal(:,2));
            end;
            K.omega=ones(length(Opt.numReg),1)*om;
            L=addstruct(L,K);
            U=addstruct(U,T);
        end;
        T=L;
        subplot(3,1,1);
        lineplot(T.numReg,T.meanR2,'style_thickline');
        drawline(0,'dir','horz');
        subplot(3,1,2);
        lineplot(T.numReg,T.meanCorr,'style_thickline');
        drawline(0,'dir','horz');
        subplot(3,1,3);
        lineplot(T.numReg,[T.propCorrR2 T.propCorrCorr],'style_thickline','leg',{'R2','Corr'});
        drawline(0.5,'dir','horz');
        varargout={T,U};
        if (~isempty(Opt.outfile))
            save(fullfile(baseDir,Opt.outfile),'T','U');
        end;
    case 'encodeRegularise'
        % Test different dimensionalties of encoding models
        % rsa_testModelCompare('encodeRegularise','numSim',10000,'outfile','sim_encodeRegularise_Exp2.mat');
        Opt.model   = 'Model_chords.mat';
        Opt.outfile = [];
        Opt.numSim  = 30;
        Opt.numPart = 8;
        Opt.numVox  = 160;
        Opt.Omega   = 0.2;
        Opt.Omega_hat = [0.01 0.05:0.05:0.5];
        
        Opt=rsa.getUserOptions(varargin,Opt);
        
        L=[];   % Summary stats
        U=[];   % details on each simulation
        
        if (ischar(Opt.model))
            load(fullfile(baseDir,Opt.model))
        else
            M=Opt.model;
        end;
        
        numModels = length(M);
        
        % Experimental constants
        test=rsa.rdm.squareRDM(M(1).RDM);
        D.numCat = size(test,1);
        D.var_e   = 1;    % Noise variance
        D.numPart = Opt.numPart;
        D.numVox  = Opt.numVox;
        D.numSim  = Opt.numSim;
        
        % Prep some convenience matrices
        H = eye(D.numCat)-ones(D.numCat)/D.numCat;  % Centering matrix
        part  = kron([1:D.numPart]',ones(D.numCat,1));
        cond  = kron(ones(D.numPart,1),[1:D.numCat]');
        Xpart = indicatorMatrix('identity',part);
        Xcond = indicatorMatrix('identity',cond);
        
        % Prep the Variance components and regression matricies for each model
        for m=1:numModels
            % Prep variance components
            RDM{m}=rsa.rdm.squareRDM(M{m}.RDM);
            G{m}=-0.5*H*RDM{m}*H;
            % Standarsize ???
            
            % Get the design matrix for encoding models
            [V,Lam]=eig(G{m});
            [lambda{m},i]=sort(diag(Lam),'descend');
            X{m}=bsxfun(@times,V(:,i),sqrt(lambda{m}'));
            X{m}=X{m}(:,lambda{m}>eps);
        end;
        
        % Loop over different regularisation levels
        D.omega   = Opt.Omega;
        for om=Opt.Omega_hat;
            D.omega_hat=om;
            fprintf('%2.3f\n',om);
            T=[];
            % Generate data from the two models
            Y=[];
            for m=1:numModels
                [S,y,D]=rsa_testModelFit('simulate_data',M{m},D);
                S.truemodel=ones(D.numSim,1)*m;
                T=addstruct(T,S);
                Y=cat(3,Y,y);
            end;
            
            T.omega_hat = ones(D.numSim*numModels,1)*om;
            for m=1:numModels
                % Calculate the different forms of correlation,
                % Encoding model without regularisation
                for n=1:size(T.RDM,1)
                    if (mod(n,ceil(Opt.numSim/20))==0);
                        fprintf('.');
                    end;
                    [T.R2(n,m),~,~,T.corr(n,m)]= encode_crossval(Y(:,:,n),Xcond*X{m},part,'ridgeFixed','X',Xpart,'G',eye(size(X{m},2))*D.omega_hat,'sigma2',D.var_e);
                end;
                fprintf('model %d done\n',m);
            end;
            % Now calculate the proportion of cases the correct model wins
            N = size(T.truemodel,1);
            for n=1:N
                val1 = T.R2(n,:);
                trueval1 = val1(T.truemodel(n));
                val1(T.truemodel(n))=[];
                meancorrect1(n,:) = (sum(trueval1>val1) + ...
                    sum(trueval1==val1)*0.5)./(numModels-1);
                
                val2 = T.corr(n,:);
                trueval2 = val2(T.truemodel(n));
                val2(T.truemodel(n))=[];
                meancorrect2(n,:) = (sum(trueval2>val2) + ...
                    sum(trueval2==val2)*0.5)./(numModels-1);
                trueVal(n,:) = [trueval1 trueval2];
            end;
            K.omega_hat = om;
            K.propCorrR2     = mean(meancorrect1);
            K.propCorrCorr   = mean(meancorrect2);
            K.meanR2         = mean(trueVal(:,1));
            K.meanCorr       = mean(trueVal(:,2));
            K.omega=Opt.Omega;
            L=addstruct(L,K);
            U=addstruct(U,T);
        end;
        T=L;
        subplot(3,1,1);
        lineplot(T.omega_hat,T.meanR2,'style_thickline');
        drawline(0,'dir','horz');
        subplot(3,1,2);
        lineplot(T.omega_hat,T.meanCorr,'style_thickline');
        drawline(0,'dir','horz');
        subplot(3,1,3);
        lineplot(T.omega_hat,[T.propCorrR2 T.propCorrCorr],'style_thickline','leg',{'R2','Corr'});
        drawline(0.5,'dir','horz');
        varargout={T,U};
        if (~isempty(Opt.outfile))
            save(fullfile(baseDir,Opt.outfile),'T','U');
        end;
    case 'Figure_numReg' % This is the encoding model figure for different numbers of regressors
        cd(baseDir);
        methods={'encodeReg','encodePCM'};   %,'encodeRidge'
        indx = {1,[2 3],[4 5]};
        corrAx=[0.15 0.015 0.01];
        R2Ax = [0.015 0.00012 0.00008];
        for ex=1:3
            T=[];
            for m=1:length(methods)
                files=dir(sprintf('sim_%s_Exp%d*',methods{m},ex));
                for i=1:length(files)
                    R=load(files(i).name);
                    Num=size(R.U.sig_hat,1);
                    N = length(R.T.omega);
                    R.T.numSim  = ones(N,1)*Num;
                    R.T.method  = ones(N,1)*m;
                    T=addstruct(T,R.T);
                end;
            end;
            if (std(T.omega)>0.0001)
                warning('not the same simulation params');
            end;
            D=tapply(T,{'method','numReg'},{T.numSim,'sum','name','numSim'},...
                {T.propCorrR2.*T.numSim,'sum','name','numCorrR2'},...
                {T.propCorrCorr.*T.numSim,'sum','name','numCorrCorr'},...
                {T.meanR2.*T.numSim,'sum','name','sumR2'},...
                {T.meanCorr.*T.numSim,'sum','name','sumCorr'});
            D.propCorrR2   = D.numCorrR2 ./ D.numSim;
            D.propCorrCorr = D.numCorrCorr ./ D.numSim;
            D.meanR2       = D.sumR2 ./ D.numSim;
            D.meanCorr     = D.sumCorr ./ D.numSim;
            
            subplot(3,5,indx{ex});
            lineplot(D.numReg,D.meanR2,'split',D.method,'style_thickline');
            drawline(0,'dir','horz');
            set(gca,'YLim',[0 R2Ax(ex)]);
            
            
            subplot(3,5,indx{ex}+5);
            lineplot(D.numReg,D.meanCorr,'split',D.method,'style_thickline','linestyle',':');
            set(gca,'YLim',[0 corrAx(ex)]);
            
            subplot(3,5,indx{ex}+10);
            lineplot(D.numReg,[D.propCorrR2],'split',D.method,'style_thickline','leg','auto');
            hold on;
            lineplot(D.numReg,[D.propCorrCorr],'split',D.method,'style_thickline','linestyle',':');
            hold off;
            drawline(0.5,'dir','horz');
            set(gca,'YLim',[0.45 0.81]);
        end;
        set(gcf,'PaperPosition',[0 0 12 9]);
        wysiwyg;
    case 'Figure_rsaVarianceLDC'
        cd(baseDir);
        
        load(fullfile(baseDir,'Model_fiveFinger.mat'));
        subplot(3,2,1);
        rsa.fig.imageRDMs(M{1},'transformfcn','ssqrt','singleRDM',1,'colourScheme',gray);
        subplot(3,2,3);
        rsa.fig.imageRDMs(M{2},'transformfcn','ssqrt','singleRDM',1,'colourScheme',gray);
        
        load(fullfile(baseDir,'Model_chords.mat'))
        subplot(3,2,2);
        rsa.fig.imageRDMs(M{1},'transformfcn','ssqrt','singleRDM',1,'colourScheme',gray);
        subplot(3,2,4);
        rsa.fig.imageRDMs(M{2},'transformfcn','ssqrt','singleRDM',1,'colourScheme',gray);
        for ex=1:2
            
            T=[];
            S=[];
            
            % Load RSA files
            files=dir(sprintf('sim_rsa_Exp%d*',ex));
            for i=1:length(files)
                R=load(files(i).name);
                [Num,om]=pivottable(R.U.omega,[],R.U.sig_hat,'length');
                for o=1:length(om)
                    R.T.numSim(R.T.omega==om(o),1)=Num(o);
                end;
                T=addstruct(T,R.T);
            end;
            T.omega = round(T.omega,3);
            D=tapply(T,{'method','methodStr','omega'},...
                {T.propCorr.*T.numSim,'sum','name','numCorr'},...
                {T.numSim,'sum','name','numSim'},...
                'subset',T.method~=2);
            D.propCorr = D.numCorr ./ D.numSim;
            
            subplot(3,2,ex+4);
            color={[0.6 0.6 0.6],[0 0 0],[0 0 0]};
            lineplot(D.omega,D.propCorr,'split',D.methodStr,'linecolor',color,...
                'markercolor',color,'linestyle',{'-',':','-'},'linewidth',2,...
                'leg','auto','errorfcn',[]);
            set(gca,'YLim',[0.4 1]);
            drawline(0.5,'dir','horz');
        end;
        set(gcf,'PaperPosition',[0 0 6 10]);
        wysiwyg;
    case 'Figure_rsa_pcm'
        filesNames={'rsa','pcm'};
        cd(baseDir);
        for ex=1:3
            
            T=[];
            S=[];
            
            % Load RSA files
            for f=1:2
                files=dir(sprintf('sim_%s_Exp%d*',filesNames{f},ex));
                for i=1:length(files)
                    R=load(files(i).name);
                    [Num,om]=pivottable(R.U.omega,[],R.U.sig_hat,'length');
                    for o=1:length(om)
                        R.T.numSim(R.T.omega==om(o),1)=Num(o);
                    end;
                    T=addstruct(T,R.T);
                end;
                T.omega = round(T.omega,3);
                D=tapply(T,{'method','methodStr','omega'},{T.propCorr.*T.numSim,'sum','name','numCorr'},{T.numSim,'sum','name','numSim'});
                D.propCorr = D.numCorr ./ D.numSim;
                S=addstruct(S,D);
            end;
            
            subplot(1,3,ex);
            if (ex==1)
                lineplot(D.omega,D.propCorr,'split',D.methodStr,'style_thickline','leg','auto','subset',~strcmp(D.methodStr,'loglikPCM'),'errorfcn',[]);
            else
                lineplot(D.omega,D.propCorr,'split',D.methodStr,'style_thickline','subset',~strcmp(D.methodStr,'loglikPCM'),'errorfcn',[]);
            end;
            hold on;
            lineplot(D.omega,D.propCorr,'split',D.methodStr,'linestyle',':','linecolor','k','subset',strcmp(D.methodStr,'loglikPCM'),'errorfcn',[]);
            hold off;
            set(gca,'YLim',[0.4 1]);
            drawline(0.5,'dir','horz');
        end;
        set(gcf,'PaperPosition',[0 0 12 3]);
        wysiwyg;
    case 'Figure_rsa_new'
        filesNames={'rsan','pcm'};
        methodStr={'spearman','pearson','pearsonNc','pearsonSq','pearsonNcSq','cosine','loglikPCM'};
        cd(baseDir);
        for ex=1:3
            
            T=[];
            S=[];
            D=[];
            % Load RSA files
            for f=1:length(filesNames);
                files=dir(sprintf('sim_%s_Exp%d*',filesNames{f},ex));
                if (~isempty(files))
                    for i=1:length(files)
                        R=load(files(i).name);
                        [Num,om]=pivottable(R.U.omega,[],R.U.sig_hat,'length');
                        for o=1:length(om)
                            R.T.numSim(R.T.omega==om(o),1)=Num(o);
                        end;
                        T=addstruct(T,R.T);
                    end;
                    T.omega = round(T.omega,3);
                    S=tapply(T,{'method','methodStr','omega'},{T.propCorr.*T.numSim,'sum','name','numCorr'},{T.numSim,'sum','name','numSim'});
                    S.propCorr = S.numCorr ./ S.numSim;
                    D=addstruct(D,S);
                end;
            end;
            if (~isempty(D))
                for i=1:length(methodStr)
                    a=strcmp(methodStr{i},D.methodStr);
                    D.method(a)=i; 
                end;
                CAT.linecolor={'g','b','b','r','r','m','k'};
                CAT.markercolor={'g','b','b','r','r','m','k'};
                CAT.linewidth={2,2,2,2,2,2,1};
                CAT.linestyle={'-','-',':','-',':','-',':'};
                subplot(1,3,ex);
                if (ex==1)
                    lineplot(D.omega,D.propCorr,'split',D.method,'style_thickline',...
                        'leg',methodStr,'CAT',CAT,'errorfcn',[]);
                else
                    lineplot(D.omega,D.propCorr,'split',D.method,'style_thickline',...
                        'CAT',CAT,'errorfcn',[]);
                end;
                set(gca,'YLim',[0.4 1]);
                drawline(0.5,'dir','horz');
            end;
        end;
        set(gcf,'PaperPosition',[0 0 12 3]);
        wysiwyg;
     case 'Figure_rsa_prob'
        filesNames={'rsaProb','pcm'};
        methodStr={'pearson','cosine','cosineWNull','cosineWData','loglikIRLS','loglikPCM'};
        cd(baseDir);
        for ex=1:3
            
            T=[];
            S=[];
            D=[];
            % Load RSA files
            for f=1:length(filesNames);
                files=dir(sprintf('sim_%s_Exp%d*',filesNames{f},ex));
                if (~isempty(files))
                    for i=1:length(files)
                        R=load(files(i).name);
                        [Num,om]=pivottable(R.U.omega,[],R.U.sig_hat,'length');
                        for o=1:length(om)
                            R.T.numSim(R.T.omega==om(o),1)=Num(o);
                        end;
                        T=addstruct(T,R.T);
                    end;
                    T.omega = round(T.omega,3);
                    S=tapply(T,{'method','methodStr','omega'},{T.propCorr.*T.numSim,'sum','name','numCorr'},{T.numSim,'sum','name','numSim'});
                    S.propCorr = S.numCorr ./ S.numSim;
                    D=addstruct(D,S);
                end;
            end;
            if (~isempty(D))
                for i=1:length(methodStr)
                    a=strcmp(methodStr{i},D.methodStr);
                    D.method(a)=i; 
                end;
                CAT.linecolor={'b','m','r','r','r','k'};
                CAT.markercolor={'b','m','r','r','r','k'};
                CAT.linewidth={2,2,2,2,2,1};
                CAT.linestyle={'-','-',':','--','-',':'};
                subplot(1,3,ex);
                if (ex==1)
                    lineplot(D.omega,D.propCorr,'split',D.method,'style_thickline',...
                        'leg',methodStr,'CAT',CAT,'errorfcn',[]);
                else
                    lineplot(D.omega,D.propCorr,'split',D.method,'style_thickline',...
                        'CAT',CAT,'errorfcn',[]);
                end;
                set(gca,'YLim',[0.4 1]);
                drawline(0.5,'dir','horz');
            end;
        end;
        set(gcf,'PaperPosition',[0 0 12 3]);
        wysiwyg;
    case 'Figure_models'
        load(fullfile(baseDir,'Model_fiveFinger.mat'));
        subplot(2,3,1);
        rsa.fig.imageRDMs(M{1},'transformfcn','ssqrt','singleRDM',1);
        subplot(2,3,4);
        rsa.fig.imageRDMs(M{2},'transformfcn','ssqrt','singleRDM',1);
        
        load(fullfile(baseDir,'Model_chords.mat'))
        subplot(2,3,2);
        rsa.fig.imageRDMs(M{1},'transformfcn','ssqrt','singleRDM',1);
        subplot(2,3,5);
        rsa.fig.imageRDMs(M{2},'transformfcn','ssqrt','singleRDM',1);
        
        load(fullfile(baseDir,'START_compl.mat'));
        subplot(2,3,3);
        rsa.fig.imageRDMs(M{1},'transformfcn','ssqrt','singleRDM',1);
        subplot(2,3,6);
        rsa.fig.imageRDMs(M{7},'transformfcn','ssqrt','singleRDM',1);
    case 'Figure_regularisation'  % Figure 6 in paper
        T=[];
        files=dir('sim_encodeRegularise_Exp2*');
        for i=1:length(files)
            R=load(files(i).name);
            T=addstruct(T,R.T);
        end;
        D=tapply(T,{'omega_hat'},{'propCorrR2'},{'propCorrCorr'},{'meanR2'},{'meanCorr'});
        
        
        subplot(3,1,1);
        lineplot(D.omega_hat,D.meanR2,'style_thickline');
        set(gca,'YLim',[0 0.0006]);
        
        subplot(3,1,2);
        lineplot(D.omega_hat,D.meanCorr,'style_thickline','linestyle',':');
        set(gca,'YLim',[0 0.03],'YAxisLocation','right');
        subplot(3,1,3);
        lineplot(D.omega_hat,[D.propCorrR2 D.propCorrCorr],'style_thickline','leg',{'R2','Corr'});
        drawline(0.5,'dir','horz');
        set(gcf,'PaperPosition',[0 0 3 7])
        wysiwyg
    case 'Figure_pcm_rsa_encode'  % Figure 7 in paper
        filesNames={'encRidgeCorr','rsa','pcm'};
        cd(baseDir);
        for ex=1:3
            
            T=[];
            S=[];
            
            % Load RSA files
            for f=1:length(filesNames)
                files=dir(sprintf('sim_%s_Exp%d*',filesNames{f},ex));
                for i=1:length(files)
                    R=load(files(i).name);
                    if (isfield(R,'U'))  % Deal with missing U from last simulation
                        [Num,om]=pivottable(R.U.omega,[],R.U.sig_hat,'length');
                        for o=1:length(om)
                            R.T.numSim(R.T.omega==om(o),1)=Num(o);
                        end;
                    else
                        R.T.numSim=ones(length(R.T.omega),1)*1000;
                    end;
                    T=addstruct(T,R.T);
                end;
                T.omega = round(T.omega,3);
                D=tapply(T,{'method','methodStr','omega'},{T.propCorr.*T.numSim,'sum','name','numCorr'},{T.numSim,'sum','name','numSim'});
                D.propCorr = D.numCorr ./ D.numSim;
                S=addstruct(S,D);
            end;
            CAT.linecolor = {'b','g','r','g','g'};
            CAT.linestyle = {'-','-','-',':','--'};
            CAT.markercolor = {'b','g','r','g','g'};
            
            subplot(1,3,ex);
            if (ex==1)
                lineplot(D.omega,D.propCorr,'split',D.methodStr,'style_thickline',...
                    'CAT',CAT,'leg','auto','subset',~strcmp(D.methodStr,'fixed'),'errorfcn',[]);
            else
                lineplot(D.omega,D.propCorr,'split',D.methodStr,'style_thickline',...
                    'CAT',CAT,'subset',~strcmp(D.methodStr,'fixed') & D.omega<0.7,'errorfcn',[]);
            end;
            set(gca,'YLim',[0.4 1]);
            drawline(0.5,'dir','horz');
        end;
        set(gcf,'PaperPosition',[0 0 12 3]);
        wysiwyg;
    case 'Figure_bestMethods'     % Figure 8 in paper
        methods={'loglikIRLS','loglikPCM','encodePCM','encodePCMcorr'};
        methType = [1 2 3 3 4 4]';
        evalType = [1 1 1 2 1 2]';
        numMethods = length(methods);
        cd (baseDir);
        for ex=1:3
            
            T=[];
            S=[];
            D=[];
            % Load and concatinate all files
            files=dir(sprintf('sim_opt_Exp%d*',ex));
            for i=1:length(files)
                R=load(files(i).name);
                S=addstruct(S,R.T);
                T=addstruct(T,R.U);
            end;
            
            numModels = max(T.truemodel);
            % Now calculate the proportion of cases the correct model wins
            N = size(T.truemodel,1);
            for m=1:numMethods
                D.propCorr(m,1) = rsa_testModelCompare('calc_mean_correct',T.(methods{m}),T.truemodel);
                D.method(m,1)    = m;
                D.methodStr{m,1} = methods{m};
                indx = find(strcmp(S.methodStr,methods{m}));
                numSim = N/numModels/length(indx);
                D.avrgTime(m,1) = mean(S.avrgTime(indx))/numSim;
            end;
            
            % Test significance of percent correct decisions
            comp = {[1 2],[2 4],[1 4],[3 4]};
            for c=1:length(comp)
                N=length(T.truemodel);
                O(1,1) = D.propCorr(comp{c}(1))*N;
                O(2,1) = D.propCorr(comp{c}(2))*N;
                O(1,2) = (1-D.propCorr(comp{c}(1)))*N;
                O(2,2) = (1-D.propCorr(comp{c}(2)))*N;
                Prop   = (D.propCorr(comp{c}(1))+D.propCorr(comp{c}(2)))/2;
                E(1:2,1)  = Prop*N;
                E(1:2,2)  = (1-Prop)*N;
                G = 2*sum(O(:).*log(O(:)./E(:)));
                p=1-chi2cdf(G,1);
                fprintf('%s vs %s: diff=%2.2f G = %2.2f p = %2.3f\n',methods{comp{c}(1)},methods{comp{c}(2)},...
                    (D.propCorr(comp{c}(2))-D.propCorr(comp{c}(1)))*100,G,p);
            end;
            
            
            D.methType =methType(D.method);
            D.evalType =evalType(D.method);
            
            subplot(2,3,ex);
            barplot(D.methType,D.propCorr-0.6,'split',D.evalType,'subset',D.methType<4);
            set(gca,'XTickLabel',{methods{1:3},''},'YLim',[0 0.39],'YTick',[0 0.1 0.2 0.3],'YTickLabel',{'60','70','80','90'});
            
            subplot(2,3,ex+3);
            barplot(D.methType,D.avrgTime,'subset',D.evalType ==1);
            set(gca,'XTickLabel',{methods{1:3},''});
        end;
        set(gcf,'PaperPosition',[0 0 9 5]);
        wysiwyg;  
end;

function r=cosineW(A,B,Sig); % Weighted cosine similarity measure 
    % A: N x q vector 
    % B: M x q vector 
    % Sig: qxq variance matrix 
    % Output: 
    % N*M weighted inner product (cosineW) 
    [V,L]=eig(Sig); 
    if (sum(imag(diag(L)))>0.001)
        keyboard; % Should not happen if varD is correct
    end; 
    l=real(diag(L));
    sq = V*bsxfun(@rdivide,V',sqrt(l)); % Slightly faster than sq = V*diag(1./sqrt(l))*V';
    wA=A*sq; 
    wB=B*sq; 
    wA=bsxfun(@rdivide,wA,sqrt(sum(wA.^2,2)));
    wB=bsxfun(@rdivide,wB,sqrt(sum(wB.^2,2)));
    r=wA*wB';
