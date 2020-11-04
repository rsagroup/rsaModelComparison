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
        numPart=5;      % Number of partitions
        P = 20;         % Number of voxels 
        vararginoptions(varargin(2:end),{'noise','numPart','N','P'}); 
        K = size(squareform(M{1}.RDM)); 
        C=pcm_indicatorMatrix('allpairs',[1:K]); 
        
        if (isscalar(noise))
            noise = noise * eye(K); 
        end; 

        XiM = C*noise*C';
        Var= 2/P * XiM.*XiM; 
        A = cholcov(Var); 
        [V,L]=eig(Var); 
                
        l=real(diag(L));
        sq = V*bsxfun(@rdivide,V',sqrt(l)); % Slightly faster than sq = V*diag(1./sqrt(l))*V';
            
        % Check the projection of the models 1 
        model = [M{1}.RDM-mean(M{1}.RDM);M{2}.RDM-mean(M{2}.RDM)]; 
        p=model*V;
        % [l';p]

        % Check the projection of the models 1 
        model = [M{1}.RDM/sqrt(sum(M{1}.RDM.^2));M{2}.RDM/sqrt(sum(M{2}.RDM.^2))]; 
        p=model*V;
        % [l';p]

        for i=1:2
            raw   = normrnd(0,1,N,size(M{1}.RDM,2));

            % Non-Crossvalidated 
            epsilon = raw * A ;
            data = bsxfun(@plus,epsilon,M{i}.RDM); 
            data = bsxfun(@plus,data,diag(XiM)'); % Add the bias  
            r{1}=corr((data)',model'); 
            r{2}=corr((data*sq)',(model*sq)'); 

            % Crossvalidated
            epsilon = raw * A * sqrt(numPart/(numPart-1)); % Inflation of noise 
            data = bsxfun(@plus,epsilon,M{i}.RDM); 
            r{3}=corrN((data)',model'); 
            r{4}=corrN((data*sq)',(model*sq)'); 
            
            tM = i; 
            fM = 3-i; 
            
            for j=1:4 
                r{j}=round(r{j},5); 
                correct(j,i)=sum(r{j}(:,tM)>r{j}(:,fM))+0.5*sum(r{j}(:,tM)==r{j}(:,fM)); 
                model2(j,i)=sum(r{j}(:,2)>r{j}(:,1))+0.5*sum(r{j}(:,1)==r{j}(:,2)); 
            end; 
        end; 
        correct = mean(correct/N,2); 
        model2 = mean(model2/N,2); 
        varargout={correct,model2}; 
    case 'Figure_crossval_noncrossval' % Figure 4 
        CAT.linecolor={'b','r','b','r'};
        CAT.markercolor={'b','r','b','r'};
        CAT.linewidth={2,2,2,2};
        CAT.linestyle={'-','-','--','--'};
        K=4; 
        numPart=[2 3 4 5 6 8 10 12]; 
        c = 0.15; 
        noise1 = [1 c 0 0;c 1 c 0;0 c 1 c;0 0 c 1]*1.3; 
        noise = {noise1,0.5,1,1.2}; 
        ymin = [0.5 0.4 0.5 0.5];         
        ymax = [0.8 0.7 0.8 0.8]; 
        H=eye(K)-ones(K)/K; 
        C=indicatorMatrix('allpairs',[1:K]); 
        V=((C*C').^2)/2; 
        [E,l]=eig(V);  % E are the eigenvalues of the variance-covariance matrix 
        for i=1:4 
            switch(i)
                case 1
                    M{1}.RDM=[0 -0.5 0 0 0 1]*E';
                    a = M{1}.RDM(1); 
                    b = M{1}.RDM(2); 
                    M{2}.RDM=[b a b b a b];
                case 2
                    M{1}.RDM=[0 -0.5 0 0 0 1]*E';  % Differs in the size of one component
                    M{2}.RDM=[0 -0.2 0 0 0 1]*E';
                case 3 
                    M{1}.RDM=[0 -0.5 0 0 0 1]*E';  
                    M{2}.RDM=[0.1 0 0 0 0 1]*E';
                case 4 
                    M{1}.RDM=[0 -0.5 0 0 0 1]*E';  % Differs in sign on independent components 
                    M{2}.RDM=[0  0  0.5 0 0 1]*E';
            end; 
            for j=1:2 
                M{j}.Gc = -0.5*H*squareform(M{j}.RDM)*H;
            end
            
            for j=1:length(numPart) 
                [Correct(j,:),Model2(j,:)]=rsa_covarianceDistPaper('predict_bestMethods',M,...
                    'numPart',numPart(j),'noise',noise{i},'N',1000); 
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
        set(gcf,'PaperPosition',[2 2 10 3]); 
        wysiwyg; 
        figure(2); 
        set(gcf,'PaperPosition',[2 2 8 5]); 
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
