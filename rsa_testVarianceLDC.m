function varargout=rsa_testVarianceLDC(what,varargin)
% Testing for contstructng and Fitting component models
%
%
% Define optional parameter
Opt.rootPath='/Users/joern/Desktop/rsaworkshop/rsatoolbox.v.2a/Demos/Demo_component';
Opt.rootPath='/Users/joern/Talks/2015/02_Demo_component';
% Opt.rootPath='/Users/jdiedrichsen/Talks/2015/02_Demo_component';
tendigitDir='/Users/joern/Projects/FingerPattern/tendigit1';
baseDir    = '/Users/joern/Dropbox (Diedrichsenlab)/Projects/distance/distance_variance';

% tendigitDir='/Users/jdiedrichsen/Data/tendigit1';

%  Make the representational model matrices from features
switch (what)
    case 'get_region_sphere'            % Get voxel distance structure for a sphere with radius R
        r=varargin{1};
        [X,Y,Z]=ndgrid([-r:r],[-r:r],[-r:r]) ;
        D=[X(:) Y(:) Z(:)];
        d=sqrt(sum(D.*D,2));
        D=D(d<=r,:);
        C=surfing_eucldist(D',D');
        varargout={C};
    case 'get_region'                % get voxel distance structure from a real ROI: M1
        % Get the region
        R=load(fullfile(tendigitDir,'RegionOfInterest','reg_prewhiten.mat'));
        R=getrow(R,R.SN==varargin{1} & R.region==varargin{2});
        C=surfing_eucldist(R.xyz',R.xyz');
        varargout={C};
    case 'change_SPM_autocorr'        % Changes the SPM autocorrelation structure of the noise
        load(fullfile(tendigitDir,'GLM_firstlevel_fast/s02/SPM.mat'));
        
        % Remake the temporal autocorrelation structure
        t=[0:122];
        r=0.5*exp(-t)+0.5*exp(-t/40);
        xVi.V=sparse((kron(eye(8),toeplitz(r))));
        xX=SPM.xX;
        
        W          = spm_sqrtm(spm_inv(xVi.V));
        W          = W.*(abs(W) > 1e-6);
        xX.W       = sparse(W);
        %-Design space and projector matrix [pseudoinverse] for WLS
        %--------------------------------------------------------------------------
        xX.xKXs        = spm_sp('Set',spm_filter(xX.K,W*xX.X));    % KWX
        xX.xKXs.X      = full(xX.xKXs.X);
        xX.pKX         = spm_sp('x-',xX.xKXs);                     % Projector
        erdf           = spm_SpUtil('trRV',xX.xKXs);               % error df
        
        %-Use non-sphericity xVi.V to compute [effective] degrees of freedom
        %--------------------------------------------------------------------------
        xX.V           = spm_filter(xX.K,spm_filter(xX.K,W*xVi.V*W')'); % KWVW'K'
        [trRV, trRVRV] = spm_SpUtil('trRV',xX.xKXs,xX.V);          % trRV (for X)
        xX.trRV        = trRV;                                     % <R'*y'*y*R>
        xX.trRVRV      = trRVRV;                                   %-Satterthwaite
        xX.erdf        = trRV^2/trRVRV;                            % approximation
        xX.Bcov        = xX.pKX*xX.V*xX.pKX';                      % Cov(beta)
        
        SPM.xX = xX;
        SPM.xVi= xVi;
        
        save(fullfile(baseDir,'tendigitSPM.mat'),'SPM');
    case 'simulate_data_ts'             % Generate data with covariance between regressors in each run
        SPM    = varargin{1};          % SPM Design matrix
        C.dist = varargin{2};       % Voxel-to-voxel difference
        D.numCond = 10;             % Number of conditions per run
        D.numPart = 8;              % Number of partitions
        D.s_a     = 0.001;          % Width of signal kernel
        D.s_e     = 0.001;          % Width of noise kernel
        D.var_e   = 1;              % Variance of the nosie
        D.goalG   = zeros(1,100);   % Goal for G
        D.numSim  = 120;            % Number of simulation
        D.Sw      = [];
        t=[0:122];
        D.tempCorr= 0.5*exp(-t)+0.5*exp(-t/40);
        
        D=rsa.getUserOptions({varargin{3:end}},D);
        
        D.numVox = size(C.dist,1);
        % Z-matrix
        part = kron([1:D.numPart]',ones(D.numCond,1));
        condition = kron(ones(D.numPart,1),[1:D.numCond]');
        Z=indicatorMatrix('identity',condition);
        
        goalG = reshape(D.goalG,D.numCond,D.numCond);
        con = indicatorMatrix('allpairs',[1:10]);
        
        [D.N,D.Q]= size(SPM.xX.X);
        D.Tm = SPM.nscan(1); % Assume scans of same length
        
        if (isempty(D.Sw))
            C.SigA=exp(-C.dist.^2/(D.s_a).^2);
            C.SigE=exp(-C.dist.^2/(D.s_e).^2);
            D.numVox=size(C.dist,1);
        else
            C.SigA=D.Sw;
            C.SigE=D.Sw;
            D.numVox=size(D.Sw,1);
            D.Sw=D.Sw(:)';
        end;
        
        errV   = toeplitz(D.tempCorr); % Temporal autocorrelation matrix
        C.cholSigA = cholcov(C.SigA);
        C.cholSigE = cholcov(C.SigE);
        C.cholG    = cholcov(goalG);
        C.cholE    = cholcov(errV);
        
        trueU  = mvnrnd_exact(goalG,D.numVox);
        trueG  = trueU*trueU'./D.numVox;
        trueD  = diag(con*trueG*con')';
        
        for n=1:D.numSim
            for r=1:D.numPart;
                err([1:D.Tm]+(r-1)*D.Tm,:)=C.cholE'*randn(D.Tm,D.numVox)*C.cholSigE*sqrt(D.var_e);
            end;
            Y(:,:,n)=SPM.xX.X*[Z*trueU;zeros(D.numPart,D.numVox)]+err;
        end;
        
        varargout = {Y,part,condition,D,C};
    case 'test_prewhiten'               % Test the LDC calculating with spatially correlated data
        load(fullfile(baseDir,'tendigitSPM.mat'));
        Rdist = rsa_testVarianceLDC('get_region_sphere',4);
        shrinkage = varargin{1};
        SE=[0.1 0.4 0.7 1 1.3 1.6];
        T=[];
        VAR_A = [0];
        
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
                        's_e',se,'s_a',se,'goalG',D.goalG(:)','numSim',1000);
                    for i=1:size(Y,3);
                        [S.dist(i,:),Sw,S.effVox(i,1),S.shrink(i,1),S.trSS(i,:)]=rsa_distanceLDCsep3(Y(:,:,i),SPM,condition,shrinkage);
                        S.Sig(i,:)=Sw(:)';
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
    case 'test_prewhiten_shrink'        % Test the LDC variance with spatially correlated data -  for different levels of prewhitening
        shrink = [0.2 0.4 0.6 1];
        for i=1:length(shrink)
            subplot(2,ceil(length(shrink)/2),i);
            rsa_testVarianceLDC('test_prewhiten',shrink(i));
            title(sprintf('Shrink: %2.2f',shrink(i)));
        end;
        set(gcf,'PaperPosition',[2 2 6 6]);
        wysiwyg;
    case 'test_sign_naturalSW'     % Test the LDC using a natural Sw estimate (from tendigit1), and natural noise.
        load(fullfile(baseDir,'tendigitSPM.mat'));
        load(fullfile(tendigitDir,'RegionOfInterest','reg_Sw_2_2.mat'));
        Rdist = rsa_testVarianceLDC('get_region',2,2);
        shrinkage = varargin{1};
        
        U = normrnd(0,1,10,20);
        G = U*U';
        G = G ./ mean(diag(G));
        
        % Contrast matrix
        Con = rsa.util.indicatorMatrix('allpairs',[1:10]);
        
        % Contrast to be tested
        c_all = [ones(45,1)];
        hand = [1 1 1 1 1 2 2 2 2 2];
        c1 = rsa.rdm.vectorizeRDM(bsxfun(@and,hand==1,hand'==1))'; 
        c2 = rsa.rdm.vectorizeRDM(bsxfun(@and,hand==2,hand'==2))'; 
        GN = G * 0;
        D.goalG   = GN(:)';
        T=[]; 
        for i=1:100
            fprintf('%d\n',i);
            [Y,part,condition,D,C]=rsa_testVarianceLDC('simulate_data_ts',SPM,Rdist,...
                'Sw',Sw,'goalG',D.goalG(:)','numSim',100);
            for i=1:size(Y,3);
                [S.dist(i,:),SSe,S.effVox(i,1),S.shrink(i,1),S.trSS(i,:)]=rsa_distanceLDCsep3(Y(:,:,i),SPM,condition,shrinkage);
                S.Sig(i,:)=SSe(:)';
                Vpred =   rsa_varianceLDC(zeros(45,1),Con,SSe,D.numPart,S.effVox(i));
                S.zDist(i,:) = S.dist(i,:)./sqrt(diag(Vpred))';
                S.zOverall(i,1)= S.dist(i,:)*c_all./sqrt(c_all'*Vpred*c_all);
                S.zHand(i,1)= S.dist(i,:)*c1./sqrt(c1'*Vpred*c1);
                S.zHand(i,2)= S.dist(i,:)*c2./sqrt(c2'*Vpred*c2);
                S.Vpred(i,:)=Vpred(:)';
            end;
            T=addstruct(T,S);
        end;
        varargout={T};
    case 'test_sign_diff'     % Test the LDC using a natural Sw estimate (from tendigit1), and natural noise.
        load(fullfile(baseDir,'tendigitSPM.mat'));
        load(fullfile(tendigitDir,'RegionOfInterest','reg_Sw_2_2.mat'));
        Rdist = rsa_testVarianceLDC('get_region',2,2);
        shrinkage = varargin{1};
      
        
        % Contrast matrix
        Con = rsa.util.indicatorMatrix('allpairs',[1:10]);
        
        % Contrast to be tested
        digit = [1 2 3 4 5 1 2 3 4 5]; 
        hand = [1 1 1 1 1 2 2 2 2 2];
        GN = eye(10);
        D.goalG   = GN(:)';
        T=[]; 
        for i=1:2
            fprintf('%d\n',i);
            [Y,part,condition,D,C]=rsa_testVarianceLDC('simulate_data_ts',SPM,Rdist,...
                'Sw',Sw,'goalG',D.goalG(:)','numSim',100);
            for i=1:size(Y,3);
                [S.dist(i,:),SSe,S.effVox(i,1),S.shrink(i,1),S.trSS(i,:)]=rsa_distanceLDCsep3(Y(:,:,i),SPM,condition,shrinkage);
                S.Sig(i,:)=SSe(:)';
                Vpred1 =   rsa_varianceLDC(zeros(45,1),Con,SSe,D.numPart,S.effVox(i));
                md=mean(S.dist(i,:)); % Compute mean distance  
                Vpred2 =   rsa_varianceLDC(ones(45,1)*md,Con,SSe,D.numPart,S.effVox(i));
                for j=2:45; 
                    c=zeros(45,1); 
                    c(1)=1; 
                    c(j)=-1; 
                    S.zDiff1(i,j-1) = S.dist(i,:)*c./sqrt(c'*Vpred1*c);
                    S.zDiff2(i,j-1) = S.dist(i,:)*c./sqrt(c'*Vpred2*c);
                end; 
                S.Vpred1(i,:) = Vpred1(:)'; 
                S.Vpred2(i,:) = Vpred2(:)'; 
            end;
            T=addstruct(T,S);
        end;
        varargout={T};
    case 'test_prewhiten_ROIsize'       % Test the LDC calculating with spatially correlated data
        load(fullfile(tendigitDir,'GLM_firstlevel_run/s02/SPM.mat'));
        
        SE     = [0.001 0.1 0.5 1 1.5 2 3];
        radius = [3 4 5];
        T=[];
        VAR_A = [0];
        
        U = normrnd(0,1,10,20);
        G = U*U';
        G = G ./ mean(diag(G));
        
        Con = rsa.util.indicatorMatrix('allpairs',[1:10]);
        
        for var_a=VAR_A;
            for r=radius
                Rdist = rsa_testVarianceLDC('get_region_sphere',r);
                
                GN = G * var_a;
                D.goalG   = GN(:)';
                for se = SE
                    for n=1
                        [Y,part,condition,D,C]=rsa_testVarianceLDC('simulate_data_ts',SPM,Rdist,...
                            's_e',se,'s_a',se,'goalG',D.goalG(:)','numSim',100);
                        for i=1:size(Y,3);
                            [S.dist(i,:),Sw,S.effVox(i,1),S.shrink(i,1),S.trSS(i,:)]=rsa_distanceLDCsep3(Y(:,:,i),SPM,condition);
                            S.Sig(i,:)=Sw(:)';
                        end;
                        D.meanDist  = mean(S.dist);
                        D.varDist   = var(S.dist);
                        Vest        = cov(S.dist);
                        D.covDist   = Vest(:)';
                        D.meanSig   = mean(S.Sig);
                        D.effVox    = mean(S.effVox);
                        D.shrink    = mean(S.shrink);
                        D.trSS      = mean(S.trSS);
                        D.numVox    = size(Rdist,1);
                        
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
        end;
        varargout={T};
    case 'test_prewhiten_ROIsize_plot'
        D=varargin{1};
        nV=unique(D.numVox);
        figure(1);
        for i=1:length(nV);
            subplot(length(nV),1,i);
            lineplot(D.trSS(:,1)./D.numVox,[D.mvarDist D.mvarPred1 D.mvarPred2],'style_thickline',...
                'leg',{'data','naive','corr Sw'},'subset',D.numVox==nV(i));
            set(gca,'XLim',[1 9]);
        end;
        figure(2);
        lineplot(D.s_e,D.trSS(:,1)./D.numVox,'split',D.numVox,'leg','auto','style_thickline')
    case 'varianceEstimator'         % Checks for measure that can serve to estimate the residual variance
        R=load(fullfile(tendigitDir,'RegionOfInterest','reg_prewhiten.mat'));
        R=getrow(R,R.SN==2 & R.region==2);
        D.numVox= size(R.SN,1);
        D.N = 100;
        D.s_e     = 50;
        
        
        D.var_e   = 1;
        D.numSim  = 60;
        
        D=rsa.getUserOptions(varargin,D);
        C.dist=surfing_eucldist(R.xyz',R.xyz');
        
        % Z-matrix
        C.SigE=exp(-C.dist.^2/D.s_e);
        C.cholSigE=cholcov(C.SigE);
        X=randn(D.N,D.numVox)*C.cholSigE*sqrt(D.var_e);
        X=bsxfun(@minus,X,mean(X));
        
        
        [u, S, v] = svd(X,0);
        s         = diag(S).^2;
        j         = find(s*length(s)/sum(s) > 1e-2);
        v         = v(:,j);
        u         = u(:,j);
        S         = S(j,j);
        
        keyboard;
    case 'LDCraw_LDCsvd'             % Compare distance prewhitened over covdiag and over svd
        
        load(fullfile(tendigitDir,'GLM_firstlevel_run/s02/SPM.mat'));
        R=load(fullfile(tendigitDir,'RegionOfInterest','reg_prewhiten.mat'));
        
        D.numCond = 10;
        D.numPart = 8;
        D.numVox  = [];
        D.s_a     = 10;
        D.s_e     = 20;
        D.var_e   = 3;
        X=normrnd(0,1,5,D.numCond);
        goalG     = X'*X/100;
        D.numSim  = 200;
        
        D=rsa.getUserOptions(varargin,D);
        
        
        % Get the region
        R=getrow(R,R.SN==2 & R.region==2);
        if (isempty(D.numVox))
            D.numVox = size(R.SN,1);
        else
            R=getrow(R,[1:D.numVox]');
        end;
        C.dist=surfing_eucldist(R.xyz',R.xyz');
        
        % Z-matrix
        part = kron([1:D.numPart]',ones(D.numCond,1));
        condition = kron(ones(D.numPart,1),[1:D.numCond]');
        Z=indicatorMatrix('identity',condition);
        
        [D.N,D.Q]= size(SPM.xX.X);
        
        C.SigA=exp(-C.dist.^2/D.s_a);
        C.SigE=exp(-C.dist.^2/D.s_e);
        C.cholSigA=cholcov(C.SigA);
        C.cholSigE=cholcov(C.SigE);
        C.cholG=cholcov(goalG);
        C.cholST=kron(C.cholG,C.cholSigA);
        
        % Determine trueG-matrix
        trueU=randn(1,size(C.cholST,1))*C.cholST;
        trueU=reshape(trueU,D.numVox,D.numCond)';
        trueG=trueU*trueU'./D.numVox;
        
        con = indicatorMatrix('allpairs',[1:10]);
        trueD = diag(con*trueG*con')';
        
        for n=1:D.numSim
            if (mod(n,10)==0)
                fprintf('.');
            end;
            err=randn(D.N,D.numVox)*C.cholSigE*sqrt(D.var_e);
            Y=SPM.xX.X*[Z*trueU;zeros(D.numPart,D.numVox)]+err;
            [wBeta,resMS,Svox,Beta]=rsa_noiseNormalizeBeta(Y,SPM);
            tBeta=bsxfun(@rdivide,Beta,sqrt(resMS));
            [S.dist1(n,:),Sw]=rsa_distanceLDC(tBeta,part,condition);
            [S.dist2(n,:),Sw]=rsa_distanceLDC(wBeta,part,condition);
            [S.dist3(n,:)]=rsa_distanceLDCraw(Y,SPM,condition);
            [S.dist4(n,:)]=rsa_distanceLDCraw1(Y,SPM,condition);
            [S.dist5(n,:),Sw,S.numVox(n,1)]=rsa_distanceLDCsvd(Y,SPM,condition);
        end;
        fprintf('\n');
        L.dist=[];
        L.method=[];
        for i=1:5
            f=sprintf('dist%d',i);
            CV(i)=mean(std(S.(f))./mean(S.(f)));
            R2(i)=mean(S.(f)*trueD'./sqrt(sum(S.(f).*S.(f),2).*(trueD*trueD')));
            L.dist=[L.dist;S.(f)];
            L.method=[L.method;ones(size(S.(f),1),1)*i];
        end;
        S.trueD = trueD;
        subplot(3,1,1);
        traceplot([1:size(L.dist,2)],L.dist,'split',L.method,'leg','auto');
        subplot(3,1,2);
        barplot([],CV);
        subplot(3,1,3);
        barplot([],1-R2);
        varargout={S,CV,R2};
    case 'LDCraw_shrinkage'             % Compare distances clalculated at different levels of shrinkage
        load(fullfile(tendigitDir,'GLM_firstlevel_run/s02/SPM.mat'));
        R=load(fullfile(tendigitDir,'RegionOfInterest','reg_prewhiten.mat'));
        
        D.numCond = 10;
        D.numPart = 8;
        D.numVox  = [];
        D.s_a     = 10;
        D.s_e     = 5;
        D.var_e   = 8;
        D.numSim  = 50;
        T=[];
        D=rsa.getUserOptions(varargin,D);
        Shrink = [ 0:0.05:1];
        
        % Get the region
        R=getrow(R,R.SN==2 & R.region==2);
        if (isempty(D.numVox))
            D.numVox = size(R.SN,1);
        else
            R=getrow(R,[1:D.numVox]');
        end;
        C.dist=surfing_eucldist(R.xyz',R.xyz');
        
        % Z-matrix
        part = kron([1:D.numPart]',ones(D.numCond,1));
        condition = kron(ones(D.numPart,1),[1:D.numCond]');
        Z=indicatorMatrix('identity',condition);
        
        [D.N,D.Q]= size(SPM.xX.X);
        
        C.SigA=exp(-C.dist.^2/D.s_a);
        C.SigE=exp(-C.dist.^2/D.s_e);
        C.cholSigA=cholcov(C.SigA);
        C.cholSigE=cholcov(C.SigE);
        
        
        con = indicatorMatrix('allpairs',[1:10]);
        
        for n=1:D.numSim
            if (mod(n,10)==0)
                fprintf('.');
            end;
            
            % Make the signal new every time
            X=normrnd(0,1,5,D.numCond);
            goalG     = X'*X/100;
            C.cholG=cholcov(goalG);
            C.cholST=kron(C.cholG,C.cholSigA);
            % Determine trueG-matrix
            trueU=randn(1,size(C.cholST,1))*C.cholST;
            trueU=reshape(trueU,D.numVox,D.numCond)';
            trueG=trueU*trueU'./D.numVox;
            trueD = diag(con*trueG*con')';
            
            err=randn(D.N,D.numVox)*C.cholSigE*sqrt(D.var_e);
            Y=SPM.xX.X*[Z*trueU;zeros(D.numPart,D.numVox)]+err;
            for s = 1:length(Shrink);
                S.cond   = s;
                S.shrink = Shrink(s);
                S.dist=rsa_distanceLDCraw1(Y,SPM,condition,Shrink(s));
                S.R2  =  (S.dist*trueD'./sqrt(sum(S.dist.*S.dist,2).*(trueD*trueD')));
                T = addstruct(T,S);
            end;
            S.cond = length(Shrink)+1;
            [S.dist,S.shrink]=rsa_distanceLDCraw1(Y,SPM,condition);
            S.R2  =  (S.dist*trueD'./sqrt(sum(S.dist.*S.dist,2).*(trueD*trueD')));
            T = addstruct(T,S);
        end;
        
        
        
        fprintf('\n');
        for c=1:max(T.cond)
            j = find(T.cond==c);
            d = T.dist(j,:);
            R.CV(c,1)=mean(std(d)./mean(d));
            R.Std(c,1)=mean(std(d));
            R.R2(c,1)=mean(T.R2(j,:));
            R.meanD(c,1) = mean(mean(d));
            R.sh(c,1) = mean(T.shrink(j));
        end;
        
        subplot(4,1,1);
        lineplot(R.sh(1:end-1),R.Std(1:end-1));
        hold on;
        plot(R.sh(end),R.Std(end),'*');
        hold off;
        ylabel('STD');
        
        subplot(4,1,2);
        lineplot(R.sh(1:end-1),R.meanD(1:end-1));
        hold on;
        plot(R.sh(end),R.meanD(end),'*');
        hold off;
        ylabel('mean Distance');
        
        subplot(4,1,3);
        lineplot(R.sh(1:end-1),R.CV(1:end-1));
        hold on;
        plot(R.sh(end),R.CV(end),'*');
        hold off;
        ylabel('coefficient of variation');
        
        subplot(4,1,4);
        lineplot(R.sh(1:end-1),1-R.R2(1:end-1));
        hold on;
        plot(R.sh(end),1-R.R2(end),'*');
        hold off;
        ylabel('1-R2');
        
        varargout={T,R};
    case 'test_ldc_svd'
        X=normrnd(0,1,10,3);
        X=bsxfun(@minus,X,mean(X));
        V = X'*X/size(X,1);
        Y=normrnd(0,1,10,3);
        
        Y=bsxfun(@minus,Y,mean(Y));
        Y1 = Y*(V^(-1/2));
        [u,s,u]=svd(V,0);
        Y2 = bsxfun(@rdivide,Y*u,sqrt(diag(s)'))*u';
        keyboard;
    case 'LDCvariance'             % Test if we can estimate the dependence across voxels
        load(fullfile(tendigitDir,'GLM_firstlevel_run/s02/SPM.mat'));
        
        
        if (trueD*trueD')==0
            d_hat   = zeros(1,45);
        else
            d_hat   = trueD*((trueD*trueD')\trueD*mean(S.dist)'); % Predicted distnance from model
        end;
        var_hat = mean(S.varB);
        Sig = reshape(var_hat,D.numCond,D.numCond);
        V1 = rsa_varianceLDC(d_hat,con,Sig,D.numPart,D.numVox);
        V2 = cov(S.dist);
        
        subplot(3,2,1);
        scatterplot(trueD',mean(S.dist)');
        
        m=max([V1(:);V2(:)]);
        subplot(3,2,3);
        imagesc_rectangle(V1,'YDir','reverse','scale',[0 m]);
        title('predicted distance');
        subplot(3,2,4);
        imagesc_rectangle(V2,'YDir','reverse','scale',[0 m]);
        title('realdistance');
        
        subplot(3,2,5);
        scatterplot(diag(V1),diag(V2),'identity','regression','linear','intercept',0,...
            'xaxisIncl',0,'yaxisIncl',0);
        title('variance');
        
        subplot(3,2,6);
        V1s=V1-diag(diag(V1));
        V2s=V2-diag(diag(V2));
        scatterplot(V1s(:),V2s(:),'identity','regression','linear','intercept',0);
        title('covariance');
        
        x=V1(:);
        y=V2(:);
        D.slope=(x'*x)\(x'*y);
        D.meanDist = mean(S.dist);
        D.predVar  = diag(V1)';
        D.simVar   = diag(V2)';
        
        varargout={D};
    case 'realVarianceBeta'             % Compare beta estimate based on X'X and on the betas from SH1/tendigit experiment
        
        % Switch between different experiments
        %         baseDir='/Users/joern/Projects/sh1';
        %         D=load(fullfile(baseDir,'RegionOfInterest/Allsub_prewhitened_data_all.mat');
        %         for s=1:12
        %             load(baseDir,'GLM_firstlevel_1',D.allDesign.name(s,:));
        %
        %         end;
        D=load(fullfile(tendigitDir,'RegionOfInterest','reg_prewhiten.mat'));
        
        part = kron([1:8]',ones(10,1));
        condition = kron(ones(8,1),[1:10]');
        
        for s=1:6
            load(fullfile(tendigitDir,'GLM_firstlevel_run',sprintf('s%2.2d',s),'SPM.mat'));
            
            reg=unique(D.region)';
            for r=reg
                i=find(D.region==r & D.SN==s);
                [d,Sw]=rsa_distanceLDC(D.wBeta(i,:)',part,condition);
                % Get average within-run covariance matrix
                for p=1:D.numPart
                    indx = part==p;
                    CV(:,:,p)=SPM.xX.Bcov(indx,indx);
                end;
                keyboard;
            end;
        end;
    case 'simulateVarianceBeta'         % Simulate variance of perwhitened Beta's
        % One reason that the predicted variance inv(X'X) does not fit the
        % observed variance across runs very well is that the noise model
        % is too optimistic - i.e. the temporal prewhitening may not fully
        % account for the real temporal co-dependency. Furthermore, true
        % inter-run variability of the activation is excluded.
        load(fullfile(tendigitDir,'GLM_firstlevel_run/s02/SPM.mat'));
        D.numVox = 100;
        D.numCond =10;
        D.numPart= 8;
        [N,Q]= size(SPM.xX.X);
        trueU = normrnd(0,1,10,D.numVox);
        part = kron([1:D.numPart]',ones(D.numCond,1));
        condition = kron(ones(D.numPart,1),[1:D.numCond]');
        Z=indicatorMatrix('identity',condition);
        Y=SPM.xX.X*[Z*trueU;zeros(D.numPart,D.numVox)]+normrnd(0,1,N,D.numVox);
        
        % Now do from scratch
        X=SPM.xX.xKXs.X;
        
        % The filtering of the data introduces some noise-dependency,
        % so that ultimately the variance estimation is off!
        KWY=spm_filter(SPM.xX.K,SPM.xX.W*Y);
        Beta1 = SPM.xX.pKX*KWY;
        res   = KWY-X*Beta1;
        
        resMS1 = sum(res.^2)./SPM.xX.erdf;
        Beta1 = Beta1(SPM.xX.iC,:);
        Bcov  = inv(X'*X);
        tBeta1=bsxfun(@rdivide,Beta1,sqrt(resMS1));
        
        
        [wBeta2,resMS2,Svox,Beta2]=rsa_noiseNormalizeBeta(Y,SPM);
        tBeta2=bsxfun(@rdivide,Beta2,sqrt(resMS2));
        [d,Sw]=rsa_distanceLDC(tBeta1,part,condition);
        % H=eye(D.numCond)-ones(D.numCond)/D.numCond;
        % Bcov=SPM.xX.Bcov;
        % Bcov=inv(SPM.xX.X'*SPM.xX.X);
        for p=1:D.numPart
            indx = part==p;
            CV(:,:,p)=Bcov(indx,indx); % *mean(resMS);
        end;
        CV = mean(CV,3);
        mean(diag(Sw./CV))
        varargout={Sw./CV};
    case 'test_LDCint_LDCsep'           % Compares the variance of the integrated and seperate LDC
        load(fullfile(tendigitDir,'GLM_firstlevel_fast/s02/SPM.mat'));
        Rdist = rsa_testVarianceLDC('get_region_sphere',5);
        D.numSim = 100;
        D.var_a  = 0.01;
        U = normrnd(0,1,10,20);
        G = U*U';
        G = G ./ mean(diag(G)) * D.var_a;
        
        D.goalG   = G(:)';
        D.distfcn = {'rsa_distanceLDCint','rsa_distanceLDCsep'};
        [S,D]     = rsa_testVarianceLDC('simulate_data_ts',SPM,Rdist,D);
        
        numDist=size(S.dist1,2);
        subplot(2,1,1);
        plot([1:numDist],mean(S.dist1),'k',[1:numDist],mean(S.dist2) ,'r'); % ,[1:numDist],mean(S.RDM3),'b');
        subplot(2,1,2);
        plot([1:numDist],std(S.dist1),'k',[1:numDist],std(S.dist2),'r'); % ,[1:numDist],std(S.RDM3),'b');
        
end;



