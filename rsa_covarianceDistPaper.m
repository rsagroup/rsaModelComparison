function varargout=rsa_covarianceDistPaper(what,varargin)
% Function that produces simulations and Figures for the paper
% "Comparing representational geometries: How to deal with bias and
% covariance of dissimilarity measures"
baseDir = '/Users/jdiedrichsen/Dropbox (Diedrichsenlab)/Projects/modelCompare';

% Use develop branch of the RSA toolbox
import rsa.*;
import rsa.util.*;
import rsa.stat.*;
import rsa.rdm.*;

% Determine color maps for RDM displays
cmapRDM = flipud(bone);
cmapCov = hot;

%  Make the representational model matrices from features
switch (what)
    % Prepare 3 different models and bring them into the correct format
    case 'make2CatModel'
        % Builds category model with n1 items in category 1 and n2 items in
        % category 2. The difference withing category is 1 and between categories 1 + betweenD
        N = [4 4];
        individVar = [1 1];
        commonVar = [1 1];
        vararginoptions(varargin,{'N','individVar','commonVar'});
        Z = [[ones(N(1),1) zeros(N(1),1)]*sqrt(commonVar(1));...
            [zeros(N(2),1) ones(N(2),1)]*sqrt(commonVar(2))];
        
        % Build the appropriate model
        M.type = 'fixed';
        M.numGparams= 0;
        M.Gc = Z*Z'+blockdiag(eye(N(1))*individVar(1),eye(N(2))*individVar(2));
        C=pcm_indicatorMatrix('allpairs',[1:sum(N)]);
        M.RDM = diag(C*M.Gc*C')';
        varargout={M};
    case 'numberOfK_1'
        % Simulation for systematically increasing the number of
        % conditions to the expense of the number of measurements
        RSA_methods={'spearman','pearson','pearsonNc','pearsonWNc','cosine','cosineWNull','loglikPCM'};
        load Model_fiveFinger.mat;
        nC = 5;
        TT=[];
        for i=1:2 % Five different levels of
            fac = 2^(i-1);
            numCond = nC*fac;
            numPart = 32/fac;
            
            % Build the appropriate model
            for m=1:2
                MM{m}=M{m};
                MM{m}.Gc = repmat(MM{m}.Gc,fac,fac);
                C=pcm_indicatorMatrix('allpairs',[1:numCond]);
                MM{m}.RDM = diag(C*MM{m}.Gc*C')';
            end;
            [~,T]=rsa_testModelCompare('modelCompare','model',MM,'numSim',1000,...
                'methods',RSA_methods,'Omega',0.05,'numPart',numPart);
            T.numPart = numPart*ones(size(T.omega,1),1);
            T.numCond = numCond*ones(size(T.omega,1),1);
            TT=addstruct(TT,T);
        end;
        T=TT;
        if (nargout==0)
            save(sprintf('sim_numberOfK_%s.mat',expStr,'T'));
        end;
        varargout={T};
    case 'numberOfK_2'
        % Simulation for systematically increasing the number of
        % conditions to the expense of the number of measurements
        RSA_methods={'spearman','pearson','pearsonNc','pearsonWNc','cosine','cosineWNull','loglikPCM'};
        TT=[];
        for i=1:4 % Five different levels of
            numItems = 2^i; % Number of items per category
            numCond = numItems*2;
            numPart = 128/numCond;
            [~,T]=rsa_testModelCompare('modelCompare','model',MM,'numSim',1000,...
                'methods',RSA_methods,'Omega',0.01,'numPart',numPart);
            T.numPart = numPart*ones(size(T.omega,1),1);
            T.numCond = numCond*ones(size(T.omega,1),1);
            TT=addstruct(TT,T);
        end;
        T=TT;
        if (nargout==0)
            save(sprintf('sim_numberOfK_%s.mat',expStr,'T'));
        end;
        varargout={T};
    case 'use_intercept'
        % M{1}.RDM = [1 2 2];
        % M{2}.RDM = [1 4 4];
        % M{1}.RDM = [1 2 1];
        % M{2}.RDM = [1 4 1];
        %M{1}.RDM = [1 2 2 2 2 1];
        % M{2}.RDM = [1 4 4 4 4 1];
        M{1}.RDM = [0 2 2 2 2 1];
        M{2}.RDM = [0 4 4 4 4 1];
        nC =4;
        H=eye(nC)-ones(nC)/nC;
        for i=1:2
            M{i}.type = 'fixed';
            M{i}.numGparams = 0;
            M{i}.Gc = -0.5*H*squareform(M{i}.RDM)*H;
        end;
    case 'Figure1'                     % Figure 1: RDMs
        G=diag([2 1 0 0 0]);
        nC =5;
        numVox = 20;
        H=eye(nC)-ones(nC)/nC;
        C=indicatorMatrix('allpairs',[1:nC]);
        G = H*G*H;
        d = (diag(C*G*C'));
        D = squareform(d);
        covD = rsa_varianceLDC(d,C,1,5,10);
        U=mvnrnd(zeros(nC,1),G,numVox)'+normrnd(0,0.5,5,numVox);
        Gemp = H*U*U'*H'/numVox;
        demp = (diag(C*Gemp*C'));
        Demp = squareform(demp);
        scm=max([D(:);Demp(:)]);
        subplot(1,2,1);
        imagesc_rectangle(D,'YDir','reverse','scale',[0 scm],'MAP',cmapRDM);
        axis equal;
        set(gca,'YTick',[],'XTick',[]);
        subplot(1,2,2);
        imagesc_rectangle(Demp,'YDir','reverse','scale',[0 scm],'MAP',cmapRDM);
        axis equal;
        set(gca,'YTick',[],'XTick',[]);
    case 'Figure_variancebias'         % Figure 3: Variance-bias plots
        % Get the simlation across different distance levels 
        % D=rsa_testVarianceBasic('dist_covariance_sim');
        D=varargin{1};  
        rsa_testVarianceBasic('Fig_variance',D)
        dist  = [1 2 5]; % Distances to consider 
        distV = [1 12 45];  % relevant entrees in the distance matrix 
        
        % Extract the dat a for the cross-validated distance
        x=D.pEc(:,dist);
        T.trueD = x(:); 
        T.pE=x(:); 
        x=D.Ec(:,dist);
        T.E=x(:); 
        x=D.pVc(:,distV);
        T.pV=x(:); 
        x=D.Vc(:,distV); 
        T.V=x(:); 
        T.crossval = ones(length(T.trueD),1); 
        % Now repeat this for non-cross validated distance 
        S=T; 
        x=D.pE(:,dist);
        S.pE=x(:); 
        x=D.E(:,dist);
        S.E=x(:); 
        x=D.pV(:,distV);
        S.pV=x(:); 
        x=D.V(:,distV); 
        S.V=x(:); 
        S.crossval = ones(length(T.trueD),1)*2; 
        T=addstruct(T,S); 
        
        % Do the two separate subplot 
        subplot(1,2,1); 
        lineplot(T.trueD,T.E,'split',T.crossval,'style_thickline','leg',{'unbiased','biased'});
        line([0;1.2],[0;1.2],'color','k'); 
        xlabel('True distance');
        ylabel('Expected value');       
        set(gca,'YLim',[0 max(T.pE)+0.05],'XLim',[0 1.3]);        

        subplot(1,2,2); 
        lineplot(T.trueD,T.V,'split',T.crossval,'style_thickline','leg',{'unbiased','biased'});
        xlabel('True distance');
        ylabel('Variance');       
        
        set(gca,'YLim',[0 max(T.pV)+0.04],'XLim',[0 1.3]);        
        set(gcf,'PaperPosition',[2 2 6.8 3]);
        wysiwyg;
        
        
    case 'Figure_covariances'          % Figure 5: Covariance matrices
        figure;
        sc = [0 0.14];
        subplot(2,3,1);
        G=diag([0 0 0 0 0]);
        nC =5;
        H=eye(nC)-ones(nC)/nC;
        C=indicatorMatrix('allpairs',[1:nC]);
        d = (diag(C*G*C'));
        D = squareform(d);
        covD = rsa_varianceLDC(d,C,1,8,10);
        imagesc_rectangle(sqrt(covD),'YDir','reverse','MAP',cmapCov,'scale',sc);
        axis equal;
        set(gca,'YTick',[],'XTick',[]);
        
        subplot(2,3,4);
        G=diag([1 0.5 0 0 0])/20;
        d = (diag(C*G*C'));
        D = squareform(d);
        covD = rsa_varianceLDC(d,C,1,8,10);
        imagesc_rectangle(sqrt(covD),'YDir','reverse','MAP',cmapCov,'scale',sc);
        axis equal;
        set(gca,'YTick',[],'XTick',[]);
        
        subplot(2,3,[2 3 5 6]); % 20 distance Figure
        nC =20;
        colormap(cmapCov);
        C=indicatorMatrix('allpairs',[1:nC]);
        covD = rsa_varianceLDC(zeros(size(C,1),1),C,1,8,10);
        imagesc(sqrt(covD),sc);
        axis equal;
        set(gca,'YTick',[],'XTick',[]);
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
    case 'Figure_rsa_crossval'
        filesNames={'rsan','pcm'};
        methodStr={'spearman','pearson','pearsonNc','pearsonSq','pearsonNcSq','cosine','loglikPCM'};
        CAT.linecolor={'g','b','b','r','r','m','k'};
        CAT.markercolor={'g','b','b','r','r','m','k'};
        CAT.linewidth={2,2,2,2,2,2,1};
        CAT.linestyle={'-','-',':','-',':','-',':'};
        cd(baseDir);
        for ex=1:3
            D=rsa_testModelCompare('Util_summarizeFiles',ex,filesNames);
            if (~isempty(D))
                D=rsa_testModelCompare('Util_relabelMethods',D,methodStr);
                subplot(1,3,ex);
                if (ex==1)
                    lineplot(D.omega,D.propCorr,'split',D.method,'style_thickline',...
                        'leg',methodStr,'CAT',CAT,'errorfcn',[],'subset',D.method>0);
                else
                    lineplot(D.omega,D.propCorr,'split',D.method,'style_thickline',...
                        'CAT',CAT,'errorfcn',[],'subset',D.method>0);
                end;
                set(gca,'YLim',[0.4 1]);
                drawline(0.5,'dir','horz');
            end;
        end;
        set(gcf,'PaperPosition',[0 0 12 3]);
        wysiwyg;
    case 'Figure_rsa_weight'           % Figure 6
        filesNames={'rsan','pcm'};
        methodStr={'pearsonNc','cosine','pearsonWNc','cosineWNull','loglikPCM'};
        CAT.linecolor={'b','r','b','r','k'};
        CAT.markercolor={'b','r','b','r','k'};
        CAT.linewidth={2,2,2,2,1};
        CAT.linestyle={'-','-','--','--',':'};
        cd(baseDir);
        for ex=1:3
            D=rsa_testModelCompare('Util_summarizeFiles',ex,filesNames,methodStr);
            
            if (~isempty(D))
                D=rsa_testModelCompare('Util_relabelMethods',D,methodStr);
                subplot(1,3,ex);
                if (ex==1)
                    lineplot(D.omega,D.propCorr,'split',D.method,'style_thickline',...
                        'leg',methodStr,'CAT',CAT,'errorfcn',[],'subset',D.method>0);
                else
                    lineplot(D.omega,D.propCorr,'split',D.method,'style_thickline',...
                        'CAT',CAT,'errorfcn',[],'subset',D.method>0);
                end;
                set(gca,'YLim',[0.4 1]);
                drawline(0.5,'dir','horz');
            end;
        end;
        set(gcf,'PaperPosition',[0 0 12 3]);
        wysiwyg;
    case 'Figure_numberOfK'            % Figure 7
        ex=1;
        methodStr={'spearman','pearson','pearsonNc','pearsonWNc','cosine','cosineWNull','loglikPCM'};
        CAT.linecolor={'m','b','b','b','r','r','k'};
        CAT.markercolor={'m','b','b','b','r','r','k'};
        CAT.linewidth={2,2,2,2,2,2,1};
        CAT.linestyle={'-',':','-','--','-','--',':'};
        
        T=[];
        files=dir(sprintf('sim_numberOfK_Exp%d*',ex));
        for i=1:length(files)
            R=load(files(i).name);
            T=addstruct(T,R.T);
        end;
        T=rsa_testModelCompare('Util_relabelMethods',T,methodStr);
        lineplot(T.numCond,T.propCorr,'split',T.method,'style_thickline',...
            'leg',methodStr,'CAT',CAT,'errorfcn',[],'subset',T.method>0);
        set(gca,'YLim',[0.55 0.8]);
        set(gcf,'PaperPosition',[0 0 4 4]);
        wysiwyg;
    case 'predict_bestMethods' % fast way to approximate which method will work better...
        N=10000;
        M=varargin{1};  % cell array of the two models to compare
        noise = 0.1;    % Noise variance on distances (under the null)
        wNoiseDist = [];      % Noise structure assumed for whitening
        numPart=5;      % Number of partitions
        P = 20;         % Number of voxels
        vararginoptions(varargin(2:end),{'noise','numPart','N','P','wNoiseDist'});
        K = size(squareform(M{1}.RDM));
        C=pcm_indicatorMatrix('allpairs',[1:K]);
        
        if (isscalar(noise))
            noise = noise * eye(K);
        end;
        
        % This is the covariance that is used to induce noise into the
        % distances
        XiM = C*noise*C';
        Var= 2/P * XiM.*XiM;
        A = cholcov(Var);
        
        % This is the covariance that is used to prewhiten the data
        if (isempty(wNoiseDist))
            wNoiseDist = (C*C').^2; % Use i.i.d assumption for whitening
        end
        [V,L]=eig(wNoiseDist);
        l=real(diag(L));
        sq = V*bsxfun(@rdivide,V',sqrt(l)); % Slightly faster than sq = V*diag(1./sqrt(l))*V';
        
        for i=1:3
            raw   = normrnd(0,1,N,size(M{1}.RDM,2));
            
            % Non-Crossvalidated
            epsilon = raw * A ;
            if (i<3)    % If simulation 1 or 2 - add true signal
                data = bsxfun(@plus,epsilon,M{i}.RDM); % Add the true signal 
            else 
                data = epsilon; 
            end;
            data = bsxfun(@plus,data,diag(XiM)'); % Add the bias
            model = [M{1}.RDM;M{2}.RDM];
            data0 = bsxfun(@minus,data,mean(data,2));
            model0 = bsxfun(@minus,model,mean(model,2));
            
            r{1}=corrN((data0)',model0');
            r{2}=corrN((data0*sq)',(model0*sq)');
            
            % Crossvalidated
            data = epsilon * sqrt(numPart/(numPart-1)); % Inflation of noise
            if (i<3)    % If simulation 1 or 2 - add true signal
                data = bsxfun(@plus,data,M{i}.RDM); % Add the true signal 
            end; 
            
            r{3}=corrN((data)',model');
            r{4}=corrN((data*sq)',(model*sq)');
            
            if (i<3) 
                tM = i;
                fM = 3-i;
            
                for j=1:4
                    r{j}=round(r{j},5);
                    correct(j,i)=sum(r{j}(:,tM)>r{j}(:,fM))+0.5*sum(r{j}(:,tM)==r{j}(:,fM));
                    model2(j,i)=sum(r{j}(:,2)>r{j}(:,1))+0.5*sum(r{j}(:,1)==r{j}(:,2));
                end;
            else % Null model 
                for j=1:4
                    bias(j) = sum(r{j}(:,2)>r{j}(:,1))+0.5*sum(r{j}(:,1)==r{j}(:,2)); 
                end; 
            end; 
        end;
        bias = bias/N; 
        correct  = correct/N;
        correctM = mean(correct,2);
        model2 = mean(model2/N,2);
        varargout={correctM,bias,model2,correct(:,1),correct(:,2)};
    case 'Figure_crossval_noncrossval' % Figure 4
        CAT.linecolor={'b','r','b','r'};
        CAT.markercolor={'b','r','b','r'};
        CAT.linewidth={2,2,2,2};
        CAT.linestyle={'-','-','--','--'};
        K=4;
        numPart=[2 3 4 5 6 8 10 12];
        c = 0.15;
        noise1 = [1 c 0 0;c 1 c 0;0 c 1 c;0 0 c 1]*2;
        noise = {noise1,0.8,1.5,2};
        ymin = [0.55 0.45 0.55 0.55];
        ymax = [0.75 0.65 0.75 0.75];
        H=eye(K)-ones(K)/K;
        
        % Make the wNoiseDist matrix used for whitening
        C=indicatorMatrix('allpairs',[1:K]);
        wNoiseDist=((C*C').^2)/2;
        [V,l]=eig(wNoiseDist);  % E are the eigenvalues of the variance-covariance matrix
        l = real(diag(l));
        
        for i=1:4
            switch(i)
                case {1,4}
                    M{1}.RDM=[0 -0.5 0 0 0 1]*V';
                    a = M{1}.RDM(1);
                    b = M{1}.RDM(2);
                    M{2}.RDM=[b a b b a b];
                case 2
                    M{1}.RDM=[0 -0.5 0 0 0 1]*V';  % Differs in the size of one component
                    M{2}.RDM=[0 -0.2 0 0 0 1]*V';
                case 3
                    M{1}.RDM=[0 -0.5 0 0 0 1]*V';
                    M{2}.RDM=[0.1 0 0 0 0 1]*V';
                case 5
                    M{1}.RDM=[0 -0.5 0 0 0 1]*V';  % Differs in sign on independent components
                    M{2}.RDM=[0  0  0.5 0 0 1]*V';
            end;
            for j=1:2
                M{j}.Gc = -0.5*H*squareform(M{j}.RDM)*H;
            end
            
            % Check the projection of the models 1
            % Note that when the projection of the model difference lands
            % on Eigenvectors with different values, prewhitening will make
            % a difference to the model decision.
            model0 = [M{1}.RDM-mean(M{1}.RDM);M{2}.RDM-mean(M{2}.RDM)];
            % model0 = bsxfun(@rdivide,model0,sqrt(sum(model0.^2,2)));
            p=model0*V;
            fprintf('\nSimulation %d\n',i);
            fprintf('Eigenvalues:% .2f  % .2f  % .2f  % .2f  % .2f  % .2f\n',l');
            fprintf('M1 biased:  % .2f  % .2f  % .2f  % .2f  % .2f  % .2f -> % .2f\n',p(1,:),model0(1,:)*model0(1,:)');
            fprintf('M2 biased:  % .2f  % .2f  % .2f  % .2f  % .2f  % .2f -> % .2f\n',p(2,:),model0(2,:)*model0(2,:)');
            
            % Check the projection of the models 1
            model = [M{1}.RDM;M{2}.RDM];
            % model = bsxfun(@rdivide,model,sqrt(sum(model.^2,2)));
            p=model*V;
            fprintf('----------------------------------------------------\n',p(1,:));
            fprintf('M1 unbiased:% .2f  % .2f  % .2f  % .2f  % .2f  % .2f -> % .2f\n',p(1,:),model(1,:)*model(1,:)');
            fprintf('M2 unbiased:% .2f  % .2f  % .2f  % .2f  % .2f  % .2f -> % .2f\n',p(2,:),model(2,:)*model(2,:)');
            
            
            for j=1:length(numPart)
                [Correct(j,:),Bias(j,:),Model2(j,:),C1(j,:),C2(j,:)]=rsa_covarianceDistPaper('predict_bestMethods',M,...
                    'numPart',numPart(j),'noise',noise{i},'N',500000,'wNoiseDist',wNoiseDist);
            end;
            figure(1);
            subplot(2,4,i);
            lineplot(numPart',Correct(:,[1 3]),'CAT',CAT,'errorfcn',[]);
            set(gca,'YLim',[ymin(i) ymax(i)]);
            drawline(0.5,'dir','horz');
            subplot(2,4,i+4);
            lineplot(numPart',Model2(:,[1 3]),'CAT',CAT,'errorfcn',[]);
            set(gca,'YLim',[0.2 0.7]);
            drawline(0.5,'dir','horz');
            
            fprintf('Biased p(M2): %2.3f  Unbiased p(M2): %2.3f\n',mean(Bias(:,1)), mean(Bias(:,3)));
            
            figure(2);
            for m=1:2
                subplot(2,4,(m-1)*4+i);
                D=squareform(M{m}.RDM);
                mD = max(D(:))*1.4;
                imagesc_rectangle(D,'YDir','reverse','MAP',cmapRDM,'scale',[0 mD]);
                set(gca,'YTick',[],'XTick',[]);
                axis equal;
            end;
        end;
        figure(1);
        set(gcf,'PaperPosition',[2 2 10 4]);
        wysiwyg;
        figure(2);
        set(gcf,'PaperPosition',[2 2 10 5]);
        wysiwyg;
    case 'Bias_Scenario3' % Why does a bias arise in Scenario 3 - despite the fact that noise is iid?
        K=4; 
        P=10; 
        N = 100000; 
        H=eye(K)-ones(K)/K;
        C=indicatorMatrix('allpairs',[1:K]);
        wNoiseDist=((C*C').^2)/2;
        [V,l]=eig(wNoiseDist);  % E are the eigenvalues of the variance-covariance matrix
        l = real(diag(l));
        A = cholcov(wNoiseDist);            % This matrix induces the covariance structure 
        sq = V*bsxfun(@rdivide,V',sqrt(l)); % Slightly faster than sq = V*diag(1./sqrt(l))*V';
        
        M{1}.RDM=[0 -0.5 0 0 0 1]*V';
        M{2}.RDM=[0.1 0 0 0 0 1]*V';
        for j=1:2
            M{j}.Gc = -0.5*H*squareform(M{j}.RDM)*H;
        end
        
        % Check the projection of the models 1
        % Note that when the projection of the model difference lands
        % on Eigenvectors with different values, prewhitening will make
        % a difference to the model decision.
        model0 = [M{1}.RDM-mean(M{1}.RDM);M{2}.RDM-mean(M{2}.RDM)];
        % model0 = bsxfun(@rdivide,model0,sqrt(sum(model0.^2,2)));
        p=model0*V;
        fprintf('Eigenvalues:% .2f  % .2f  % .2f  % .2f  % .2f  % .2f\n',l');
        fprintf('M1 biased:  % .2f  % .2f  % .2f  % .2f  % .2f  % .2f\n',p(1,:));
        fprintf('M2 biased:  % .2f  % .2f  % .2f  % .2f  % .2f  % .2f\n',p(2,:));
                
       % Non-Crossvalidated
       for i=1:2 
           data = normrnd(0,1,N,size(model0,2)) * A;
           data = bsxfun(@plus,data,mean(model0,1));
           % data = bsxfun(@plus,data,diag(XiM)'); % Add the bias
           data0 = bsxfun(@minus,data,mean(data,2)); % Subtract mean 

           r{1}=corrN((data0)',model0');
           r{2}=corrN((data0*sq)',(model0*sq)');
           m1n(i) = mean(r{1}(:,1)>r{1}(:,2)); 
           m1w(i) = mean(r{2}(:,1)>r{2}(:,2));            
       end; 
       fprintf('For Model 1 (N): %2.3f  %2.3f   %2.3f \n',m1n(1),m1n(2),mean(m1n)); 
       fprintf('For Model 1 (W): %2.3f  %2.3f   %2.3f \n',m1w(1),m1w(2),mean(m1w)); 
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
