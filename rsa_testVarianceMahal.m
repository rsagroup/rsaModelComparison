function varargout=rsa_testVarianceMahal(what,varargin)
% Tests the variance of the squared crossvalidated Mahalanobis distance
% by simulating data from the experiment in Homework 1.
%
tendigitDir='/Users/jdiedrichsen/Dropbox (Diedrichsenlab)/Teaching/AnalysisBrainImaging/homework04';


%  Make the representational model matrices from features
switch (what)
    
    case 'get_region_sphere'            % Get a 3d sphere of voxels and computes the
        % Voxel distance structure
        r=varargin{1};          % Radius
        [X,Y,Z]=ndgrid([-r:r],[-r:r],[-r:r]) ;
        D=[X(:) Y(:) Z(:)];
        d=sqrt(sum(D.*D,2));
        D=D(d<=r,:);
        C=surfing_eucldist(D',D');      % Computes Eucledian spatial distance
        varargout={C};
    case 'get_region_M1'                % get voxel distance structure from a real ROI from homework4 (XYZ)
        % Get the region
        R=load(fullfile(tendigitDir,'RegionOfInterest','reg_prewhiten.mat'));
        R=getrow(R,R.SN==2 & R.region==1);
        if (isempty(D.numVox))
            D.numVox = size(R.SN,1);
        else
            R=getrow(R,[1:D.numVox]');
        end;
        C=surfing_eucldist(R.xyz',R.xyz');
        vargout={C};
    case 'simulate_data_ts'             % Generate data with covariance between regressors in each run
        Xtask  = varargin{1};          % Design matrix
        C.dist = varargin{2};       % Voxel-to-voxel difference
        D.numCond = 10;             % Number of conditions per run (K)
        D.numPart = 8;              % Number of partitions (M)
        D.s_a     = 0.001;          % Width of signal kernel
        D.s_e     = 0.001;          % Width of noise kernel
        D.var_e   = 1;              % Variance of the noise
        D.goalG   = zeros(1,100);   % Goal for G
        D.numSim  = 120;            % Number of simulation
        
        D=rsa.getUserOptions({varargin{3:end}},D);
        
        D.numVox = size(C.dist,1);
        % Z-matrix
        part = kron([1:D.numPart]',ones(D.numCond,1));
        condition = kron(ones(D.numPart,1),[1:D.numCond]');
        Z=indicatorMatrix('identity',condition); % in Matlab/dataframe/util/indicatorMatrix.m
        
        goalG = reshape(D.goalG,D.numCond,D.numCond);
        con = indicatorMatrix('allpairs',[1:D.numCond]);
        
        [D.N,D.Q]= size(Xtask);
        
        C.SigA=exp(-C.dist.^2/D.s_a);
        C.SigE=exp(-C.dist.^2/D.s_e);
        C.cholSigA=cholcov(C.SigA);
        C.cholSigE=cholcov(C.SigE);
        C.cholG=cholcov(goalG);
        C.cholST=kron(C.cholG,C.cholSigA);
        
        trueU  = randn(1,size(C.cholST,1));
        trueU  = trueU*C.cholST;
        trueU  = reshape(trueU,D.numVox,D.numCond)';
        trueG  = trueU*trueU'./D.numVox;
        trueGw  = trueU/C.SigE*trueU'./D.numVox;
        D.trueD  = diag(con*trueG*con')';
        D.trueDw  = diag(con*trueGw*con')';
        D.trueG  = trueG';
        D.trueGw  = trueGw';
        for n=1:D.numSim
            err=randn(D.N,size(C.cholSigE,1))*C.cholSigE*sqrt(D.var_e);
            Y(:,:,n)=Xtask*Z*trueU+err;
        end;
        D.trueU = trueU;
        varargout = {Y,part,condition,D,C};
        
    case 'test_prewhiten'               % Test the LDC calculating with spatially correlated data
        load(fullfile(tendigitDir,'GLM_firstlevel_run/s02/SPM.mat'));
        Rdist = rsa_testVarianceLDC('get_region_sphere',5);
        
        SE=[1];           % Check different noise kernels
        T=[];
        VAR_A = [0:0.2:1];
        
        U = normrnd(0,1,10,20);
        G = U*U';
        G = G ./ mean(diag(G));
        
        
        Con = rsa.util.indicatorMatrix('allpairs',[1:10]);
        
        for var_a=VAR_A;
            GN = G * var_a;
            D.goalG   = GN(:)';
            for se = SE
                for n=1
                    [Y,part,condition,D,C]=rsa_testVarianceLDC('simulate_data_ts',SPM,Rdist,...
                        's_e',se,'s_a',se,'goalG',D.goalG(:)','numSim',100);                % Simualte data
                    for i=1:size(Y,3);
                        [S.dist(i,:),Sw,S.effVox(i,1),S.shrink(i,1),S.trSS(i,:)]=rsa_distanceLDCsep3(Y(:,:,i),SPM,condition,1);
                        S.Sig(i,:)=Sw(:)';          % Calcualte the distances
                    end;
                    D.meanDist  = mean(S.dist);
                    D.varDist   = var(S.dist);
                    Vest        = cov(S.dist);
                    D.covDist   = Vest(:)';
                    D.meanSig   = mean(S.Sig);
                    D.effVox    = mean(S.effVox);
                    D.shrink    = mean(S.shrink);
                    D.trSS      = mean(S.trSS);
                    
                    
                    Sig_hat = reshape(D.meanSig,D.numCond,D.numCond);
                    Vpred1 =   rsa_varianceLDC(zeros(45,1),Con,Sig_hat,D.numPart,D.numVox);
                    D.varPred1 = diag(Vpred1)';
                    Vpred2 =   rsa_varianceLDC(zeros(45,1),Con,Sig_hat,D.numPart,D.effVox);
                    D.varPred2 = diag(Vpred2)';
                    
                    % Get the mean distance
                    D.mmeanDist=mean(D.meanDist,2);
                    D.mvarDist=mean(D.varDist,2);
                    D.mvarPred1=mean(D.varPred1,2);
                    D.mvarPred2=mean(D.varPred2,2);
                    
                    T=addstruct(T,D);
                    fprintf('.\n');
                end;
            end;
        end;
        varargout={T};
        lineplot(T.s_e,[T.mvarDist T.mvarPred1 T.mvarPred2],'style_thickline',...
            'leg',{'data','naive','Sp-hat'},'leglocation','northwest');
        ylabel('Variance of distance');
        xlabel('Smoothness of noise (kernel width in voxels)');
        set(gcf,'PaperPosition',[2 2 5 5]);
        wysiwyg
    case 'test_prewhiten_bias'               % Test the LDC calculating with spatially correlated data
        
        numRuns = 8;
        numCond = 5;
        numVox  = 20;
        T=[];
        Rdist = rsa_testVarianceLDC('get_region_sphere',3);
        
        cond = kron(ones(numRuns,1),[[1:numCond]]');
        part  =kron([1:numRuns]',ones(numCond,1));
        
        condT = kron([[1:numCond] 0]',ones(7,1));
        B = indicatorMatrix('identity_p',condT);
        X = kron(eye(numRuns),B);
        partT = kron([1:numRuns]',ones(size(condT,1),1));
        
        se=[0.00001];           % Check different noise kernels
        var_a= [0.2];
        
        U = normrnd(0,1,numCond,numCond);
        G = U*U';
        G = G ./ mean(diag(G));
        
        GN = G * var_a;
        D.goalG   = GN(:)';
        [Y,part,condition,D,C]=rsa_testVarianceMahal('simulate_data_ts',X,Rdist,...
            's_e',se,'s_a',se,'goalG',D.goalG(:)','numSim',200,'numCond',numCond,'numVox',numVox);                % Simualte data
        for n=1:size(Y,3);
            
            beta_hat=pinv(X)*Y(:,:,n);                                       %%% ordinary least squares estimate of beta_hat = inv(X'*X)*X'*Y
            res=Y(:,:,n)-X*beta_hat;                               %%% residuals: res  = Y - X*beta
            
            u_hat   = zeros(size(beta_hat));
            shrink=zeros(numRuns,1);
            for i=1:numRuns
                idxT    = partT==i;             % Time points for this partition
                idxQ    = part==i;             % Regressors for this partition
                [Sw_reg(:,:,i),shrinkage(i),Sw_hat(:,:,i)]=rsa.stat.covdiag(res(idxT,:),sum(idxT)-sum(idxQ));                    %%% regularize Sw_hat through optimal shrinkage
                % Postmultiply by the inverse square root of the estimated matrix
                [V,L]=eig(Sw_reg(:,:,i));       % This is overall faster and numerical more stable than Sw_hat.^-1/2
                %  [V,L]=eig(C.SigE); 
                l=diag(L);
                sq(:,:,i) = V*bsxfun(@rdivide,V',sqrt(l)); % Slightly faster than sq = V*diag(1./sqrt(l))*V';
                u_hat(idxQ,:)=beta_hat(idxQ,:)*sq(:,:,i);
            end;
            shrinkage=mean(shrinkage);
            Sw_hat = mean(Sw_hat,3);
            Sw_reg = mean(Sw_reg,3);
            
            S.dist(n,:)=rsa.distanceLDC(u_hat,part,cond);
        end;
        con = indicatorMatrix('allpairs',[1:numCond]); 
        trueGr  = D.trueU/Sw_reg*D.trueU'./D.numVox;
        D.trueDr  = diag(con*trueGr*con')';

        
        
        T=addstruct(T,S);
        subplot(3,1,1); 
        scatterplot(D.trueD',mean(T.dist)','identity')
        subplot(3,1,2); 
        scatterplot(D.trueDw',mean(T.dist)','identity')
        subplot(3,1,3); 
        scatterplot(D.trueDr',mean(T.dist)','identity')

        fprintf('.\n');
        varargout={T}; 
end;




