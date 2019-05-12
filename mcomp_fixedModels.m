function varargout=mcomp_fixedModels(what,varargin)
% RSA / PCM / Encoding analysis Modelcomparision examples 
% Diedrichsen & Kriegeskorte, 2017, PlosCompBio paper
% 
baseDir   = fileparts(which('mcomp_fixedModels.m')); 
resultDir = '/Users/jdiedrichsen/Dropbox (Diedrichsenlab)/Projects/modelCompare';

import rsa.*;
import rsa.util.*;
import rsa.stat.*;
import rsa.rdm.*;

%  Make the representational model matrices from features
switch (what)
    case 'simulate_data'      
        % Generate activity data from pcm model-structure 
        % 1. get G from model parameters or distances 
        % 2. get trueU (u=cholcov(G)*normrand(0,1))
        % 3. get Y (= Za*u + noise)
        % 4. get d_hat and sigma_hat by crossvalidation
        %
        % Note that G=UU'= P * sum(omega*Gc);
        % Noise*Noise' = P * var_e
        % This for a SNR of 1:1 omega needs to be same size as var_e
        Model  = varargin{1};  % Input Model (in ) 
                        
        % (Simulated) Experimental parameters: default
        D.numPart       = 9;                % Number of runs
        D.numVox        = 160;              % Number of voxels
        D.signal        = 0.1;   % hyperparamter on distance model
        D.noise         = 1;                % Noise variance: Like omega, this is cumulative over voxels
        D.numSim        = 30;             % Number of simulation
        D.theta         = 0; 
        %- allow to get user options
        D               = getUserOptions({varargin{2:end}},D);
        
        [Y,partvec,condvec]=pcm_generateData(Model,D.theta,'signal',D.signal,...
                                             'noise',D.noise,'numSim',D.numSim,...
                                             'numVox',D.numVox,'numPart',D.numPart); 
        for n=1:D.numSim 
            % Calc cross-validated distance and noise estimate
            % try
            [d_hat,Sig_hat] = rsa.distanceLDC(Y{n},partvec,condvec);
            % catch
            %     tmp1 = indicatorMatrix('identity',conditions);
            %     tmp2 = indicatorMatrix('allpairs',unique(conditions)');
            %     [d_hat,Sig_hat] = distance_ldc_sigma(Y(:,:,n),tmp1,tmp2,part);
            % end
            S.RDM(n,:) = d_hat';
            D.numCond = size(squareform(d_hat),1); 
            S.Sig_hat(n,:)= Sig_hat(tril(true(D.numCond),0))';  % Get vectorized version of the variance-covariance matrix
            S.sig_hat(n,1)= mean(diag(Sig_hat));               % Noise variance estimate
            Sigma(:,:,n) = Sig_hat;
        end;
        S.numVox = ones(size(S.RDM,1),1)*D.numVox;
        varargout = {Y, partvec,condvec,S,D};
    case 'modelCompare'
        % RSA simulations: 
        % methods={'spearman','pearson','fixed','loglikIRLS'};
        % mcomp_fixedModels('modelCompare','model','Model_fiveFinger.mat','numSim',1000,'outfile','sim_rsa_Exp1.mat','methods',methods,'Signal',[0:0.1:0.8]);
        % mcomp_fixedModels('modelCompare','model','Model_chords.mat','numSim',100,'outfile','sim_rsa_Exp2.mat','methods',methods,'Signal',[0:0.05:0.3]);
        % mcomp_fixedModels('modelCompare','model','START_compl.mat','numSim',1,'outfile','sim_rsa_Exp3.mat','methods',methods,'Signal',[0:0.1:0.8]);
        
        % PCM simulations
        %         mcomp_fixedModels('modelCompare','model','Model_fiveFinger.mat','numSim',500,'outfile','sim_pcm_Exp1.mat','methods',{'loglikPCM'},'Signal',[0:0.1:0.8]);
        %         mcomp_fixedModels('modelCompare','model','Model_chords.mat','numSim',100,'outfile','sim_pcm_Exp2.mat','methods',{'loglikPCM'},'Signal',[0:0.05:0.3]);
        %         mcomp_fixedModels('modelCompare','model','START_compl.mat','numSim',100,'outfile','sim_pcm_Exp3.mat','methods',{'loglikPCM'},'Signal',[0:0.1:0.8]);
        %
        % For best methoods in all categories: 
        % OPT_methods={'encodePCM','loglikIRLS','loglikPCM'};
        % mcomp_fixedModels('modelCompare','model','Model_fiveFinger.mat','numSim',1000,'outfile','sim_opt_Exp1a.mat','methods',OPT_methods,'Omega',0.3);
        % mcomp_fixedModels('modelCompare','model','Model_chords.mat','numSim',1000,'outfile','sim_opt_Exp2a.mat','methods',OPT_methods,'Omega',0.15);
        % mcomp_fixedModels('modelCompare','model','START_compl.mat','numSim',500,'outfile','sim_opt_Exp3a.mat','methods',OPT_methods,'Omega',0.5);
        Opt.methods ={'kendall','spearman','pearson','fixed','loglikIRLS','loglikPCM','encodeReg','encodePCM'};
        Opt.model   = 'Model_fiveFinger.mat';
        Opt.numSim  = 20;
        Opt.numPart = 8;
        Opt.numVox  = 60;
        Opt.outfile = [];
        Opt.theta   = 0; 
        Opt.signal  = [0 0.3 0.6];
        
        Opt=rsa.getUserOptions(varargin,Opt);
        
        L=[];   % Summary stats
        U=[];   % details on each simulation
        
        numMethods = length(Opt.methods);
        if (ischar(Opt.model))
            load(fullfile(baseDir,Opt.model))
        else
            M=Opt.model;
        end;
        
        numModels = length(M);
        
        % Experimental constants
        test=pcm_calculateG(M{1},[]); 
        D.numCond = size(test,1);
        D.noise   = 1;    % Noise variance
        D.numPart = Opt.numPart;
        D.numVox  = Opt.numVox;
        D.numSim  = Opt.numSim;
                
        % Prep the Variance components and regression matricies for each model
        for m=1:numModels
            M{m}.type='component'; 
            M{m}.numGparams=1; 
            
            % Prep variance components
            % Get the design matrix for encoding models
            [V,Lam]=eig(M{m}.Gc);
            [lambda{m},i]=sort(diag(Lam),'descend');
            X{m}=bsxfun(@times,V(:,i),sqrt(lambda{m}'));
            X{m}=X{m}(:,lambda{m}>eps);
            numReg(m) = size(X{m},2);
        end;
        
        % Loop over different signal levels
        for om=Opt.signal;
            D.signal = om;
            fprintf('%2.3f\n',om);
            T=[];
            % Generate data from the two models
            Y={};
            for m=1:numModels
                [y,partvec,condvec,S]=mcomp_fixedModels('simulate_data',M{m},D);
                S.truemodel=ones(D.numSim,1)*m;
                T=addstruct(T,S);
                Y=[Y(:);y(:)];
            end;
            T.signal = ones(D.numSim*numModels,1)*om;
            for m=1:numModels
                % Calculate the different forms of correlation,
                for meth = 1:numMethods
                    tic;
                    switch(Opt.methods{meth})
                        case 'kendall'
                            T.kendall(:,m)    = corr(T.RDM',M{m}.RDM','type','Kendall');
                        case 'spearman'
                            T.spearman(:,m)   = corr(T.RDM',M{m}.RDM','type','Spearman');
                        case 'pearson'
                            T.pearson(:,m)    = corr(T.RDM',M{m}.RDM');
                        case 'fixed'
                            [T.weight(:,m),T.fixed(:,m),T.loglikeFixed(:,m)] = ...
                                rsa.stat.fitModelOLS(M{m}.RDM',T.RDM);
                        case 'loglikIRLS'
                            % Likelihood under normal approximation with assumed Sigma
                            [T.weight(:,m),~,T.loglikIRLS(:,m)]=...
                                rsa_fitModelIRLS(M(m).RDM',T.RDM,T.sig_hat,8,D.numVox);
                        case 'loglikIRLSsig'
                            % Likelihood under normal approximation with
                            % inferred structure, but open sigma for
                            % scaling
                            [T.weight(:,m),~,T.loglikIRLSsig(:,m)]=...
                                rsa_fitModelIRLS(M(m).RDM',T.RDM,T.sig_hat,8,D.numVox,'assumeVoxelInd',0);
                        case 'loglikPCM'
                            % Log likelihood using PCM: 
                            [TT,theta]=pcm_fitModelIndivid(Y,M(m),partvec,condvec,'verbose',0,'fitScale',0);                                 % [~,theta1,~,la]=pcm_NR(Y(:,:,n),Xcond,'Gc',{G{m}},'X',Xpart,'hP',0);
                            T.loglikPCM(:,m) = TT.likelihood;
                            T.PCMe(:,m)=TT.noise;
                        case 'encodeReg'
                            % Encoding model without regularisation
                            for n=1:size(T.RDM,1)
                                [T.encodeReg(n,m),~,~,T.encodeRegCorr(n,m)] = encode_crossval(Y(:,:,n),Xcond*X{m}(:,1:2),part,'linregress','X',Xpart);
                            end;
                        case 'encodePCM'
                            % Encoding model with regularisation
                            for n=1:size(T.RDM,1)
                                [T.encodePCM(n,m),~,h,T.encodePCMcorr(n,m)] = encode_crossval(Y(:,:,n),Xcond*X{m},part,'pcm_NR_diag','X',Xpart,'Gd',ones(numReg(m),1));
                                T.encodePCMs(n,m)=h(1);
                                T.encodePCMe(n,m)=h(2);
                            end;
                        case 'encodeRidge'
                            % Encoding model with regularisation
                            for n=1:size(T.RDM,1)
                                [T.encodeRidge(n,m),~,~,T.encodeRidgeCorr(n,m)] = encode_crossval(Y(:,:,n),Xcond*X{m},part,'ridgeFixed','X',Xpart,'G',eye(numReg(m))*om,'sigma2',1);
                            end;
                    end;
                    exetime(m,meth) = toc;
                end;
                fprintf('model %d done\n',m);
            end;
            % Now calculate the proportion of cases the correct model wins
            for m=1:numMethods
                N = size(T.truemodel,1);                
                K.propCorr(m,1) = mcomp_fixedModels('calc_mean_correct',T.(Opt.methods{m}),T.truemodel);
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
end;

