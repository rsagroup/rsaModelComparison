function varargout=rsa_testVarianceBasic(what,varargin);
% Test for basic properties of the variance of inner products between
% random variables and of crossvalidated differences. The beta-weights are
% generally simulated directly. For a more realistic simulation of time
% series data with realistic Design matrices, etc, see rsa_testVarianceLDC.
switch(what)
    case '3d_sphere'            % Make a covariance matrix based on spatial prox
        r=varargin{1};
        [X,Y,Z]=ndgrid([-r:r],[-r:r],[-r:r]) ;
        D=[X(:) Y(:) Z(:)];
        d=sqrt(sum(D.*D,2));
        D=D(d<=r,:);
        C=surfing_eucldist(D',D');
        varargout={C};
    case 'design'               % Covariance for Kx2 design with matching similarity
        D=varargin{1};
        C=varargin{2};
        D.P=size(C.dist,1);
        
        % Make True G
        
        C.trueG=[1 0.5 0;0.5 1 0.2;0 0.2 1].*D.var_a;
        Z=kron(ones(D.b,1),eye(D.K));
        C.run=kron([1:D.b],ones(1,D.K));
        % make spatial cholinsky matrices
        C.SigA=exp(-C.dist.^2/D.s_a);
        C.SigE=exp(-C.dist.^2/D.s_e);
        C.cholSigA=cholcov(C.SigA);
        C.cholSigE=cholcov(C.SigE);
        C.cholG=cholcov(C.trueG);
        C.cholST=kron(C.cholG,C.cholSigA);
        
        % Make category membership
        varargout={D,C,Z};
    case 'design_matrix_even'   % Makes balanced design matrix
        D=varargin{1};
        X=[];
        for b=1:D.b
            tt=kron([1:D.K]',ones(D.numtrial,1));
            X=blockdiag(X,indicatorMatrix('identity',tt(randperm(size(tt,1)))));
        end;
        varargout={X};
    case 'design_matrix_uneven'
        D=varargin{1};
        X=[];
        for b=1:D.b
            tt=[];
            for k=1:D.K
                tt=[tt;ones(unidrnd(7),1)*k];
            end;
            X=blockdiag(X,indicatorMatrix('identity',tt(randperm(size(tt,1)))));
        end;
        varargout={X};
    case 'design_matrix_corr'   % Correlated design matrix
        D=varargin{1};
        X=[];
        for b=1:D.b
            tt=[];
            for k=1:D.K
                tt=[tt;ones(unidrnd(7),1)*k];
            end;
            x=indicatorMatrix('identity',tt(randperm(size(tt,1))));
            for i=1:size(x,2)
                x(:,i)=conv(x(:,i),[0.4 1],'same');
            end;
            X=blockdiag(X,x);
        end;
        varargout={X};
    case 'data'                 % Data for Kx2 design
        D=varargin{1};
        C=varargin{2};
        Z=varargin{3};
        D.N=size(Z,1);
        trueU=randn(1,size(C.cholST,1))*C.cholST;
        trueU=reshape(trueU,D.P,size(Z,2))';
        
        err=randn(D.N,D.P)*C.cholSigE*sqrt(D.var_e);
        
        y=Z*trueU+err;
        varargout={y,trueU};
    case 'vector_test'          % Test for the properties of the inner product between two vectors <x,y>
        % Case of indepdent vectors
        Mean=[0 0];
        Var_a=[0.2 0.2];
        Var_e=[1 1];
        
        P=100;                  % Number of independet voxels
        x1=normrnd(Mean(1),sqrt(Var_a(1)),P,1);    % true mean vector1
        x2=normrnd(Mean(2),sqrt(Var_a(2)),P,1);    % True mean vector2
        
        for i=1:100000
            X(:,1)=x1+normrnd(0,sqrt(Var_e(1)),P,1);
            X(:,2)=x2+normrnd(0,sqrt(Var_e(2)),P,1);
            T.c(i,1)=X(:,1)'*X(:,2);
            % You can also decompose into mean value and variance!
            % mX=mean(X);
            % X=bsxfun(@minus,X,mX);
            % T.c(i,1)=X(:,1)'*X(:,2)+mX(1)*mX(2)*P;
        end;
        mean_hat = x1'*x2;  % Predicted mean is unbiased
        var_hat  = Var_e(1)*Var_e(2)*P+(x1'*x1)*Var_e(2)+(x2'*x2)*Var_e(1);   % predicted variance
        fprintf('predicted m:%2.2f v:%2.2f\n',mean_hat,var_hat);
        fprintf('measured m:%2.2f v:%2.2f\n',mean(T.c),var(T.c));
    case 'vector_test_dep'      % Test for the properties of the inner product between two vectors dependent <x,y>
        % Case of indepdent vectors
        Mean=[0 0];
        Var_a=[0.2 0.2];
        Sig_e=[1 0.5;0.5 1]; % Variance-covariance across the two parts of the inner product
        A_e=cholcov(Sig_e);
        P=100;                  % Number of independet voxels
        xm(1,:)=normrnd(Mean(1),sqrt(Var_a(1)),1,P);    % true mean vector1
        xm(2,:)=normrnd(Mean(2),sqrt(Var_a(2)),1,P);    % True mean vector2
        N=200000;
        T.c = zeros(N,1);
        for i=1:N
            X = xm + A_e'*randn(2,P);
            T.c(i,1)=sum(X(1,:).*X(2,:));
            % You can also decompose into mean value and variance!
            % mX=mean(X);
            % X=bsxfun(@minus,X,mX);
            % T.c(i,1)=X(:,1)'*X(:,2)+mX(1)*mX(2)*P;
        end;
        mean_hat = xm(1,:)*xm(2,:)'+Sig_e(2,1)*P;  % Predicted mean is biased by covariance
        var_hat  = (xm(1,:)*xm(1,:)')*Sig_e(2,2)+(xm(2,:)*xm(2,:)')*Sig_e(1,1)+...
                    2*xm(1,:)*xm(2,:)'*Sig_e(2,1)+...
                    Sig_e(1,1)*Sig_e(2,2)*P+Sig_e(2,1)*Sig_e(2,1)*P;   % predicted variance
        fprintf('predicted m:%2.2f v:%2.2f\n',mean_hat,var_hat);
        fprintf('measured m:%2.2f v:%2.2f\n',mean(T.c),var(T.c));
    case 'matrix_test'          % Tests the properties of <X,Y> for a simple pair of independent crossvalidation folds
        D.P    = 50;            % Number of voxels
        D.K    = 3;             % Number of conditions
        D.eps  = 0.2;           % Variance for the noise
        Sig    = eye(D.K)*D.eps;
        
        D.G   = [1 0.2 0.1 0.2 0.5 0.05 0.1 0.05 0.1];         % Different variances of the constant patterns
        % D.G   = diag([1 0.5 0.2 0.1]);         % Different variances of the constant patterns
        
        % Make constant pattern that conform exactly to the desired G
        U = normrnd(0,1,D.K,D.P);
        E = (U*U')/D.P;
        G = reshape(D.G,D.K,D.K);
        A = cholcov(G)';
        U = A*E^(-0.5)*U;               % Now U*U'/P is exactly G
        
        sig= cholcov(Sig);
        for i=1:7000
            Y1=U+sig'*normrnd(0,1,D.K,D.P);         % Generate two indpendent folds
            Y2=U+sig'*normrnd(0,1,D.K,D.P);
            G_hat=Y1*Y2';                               % G-hat
            D.G_hat(i,:)=G_hat(:)';
        end;
        
        % Prediction of mean
        G  = U*U';       % This is the predicted mean
        g  = G(:);      % stretched out to a vector
        s  = Sig(:);    % Stretched out to a vector
        % Prediction of variances and covariance
        A=ones(D.K);
        V=(kron(A,G).*kron(Sig,A) + kron(G,A).*kron(A,Sig)  + kron(A,Sig).*kron(Sig,A).*D.P);
        subplot(2,1,1);
        imagesc(V);
        subplot(2,1,2);
        imagesc(cov(D.G_hat));
        keyboard;
    case 'matrix_test_center'   % Tests the properties of <X,Y> before and after mean subtraction
        P=50;               % Number of voxels
        K=4;                % Number of conditions
        V=[3 1 2 1]';       % Different variances
        VV=V*V';
        [I1,I2]=meshgrid([1:4],[1:4]);
        
        U=normrnd(0,1,K,P);
        U=bsxfun(@plus,U,normrnd(0,0,1,P)); %common activity
        
        Um=bsxfun(@minus,U,mean(U));       % Mean subtracted common activity patterns
        
        for i=1:7000
            for j=1:K
                Y1(j,:)=U(j,:)+normrnd(0,sqrt(V(j)),1,P);
                Y2(j,:)=U(j,:)+normrnd(0,sqrt(V(j)),1,P);
            end;
            G=Y1*Y2';                       % Without mean subtraction
            D.G(i,:)=G(:)';
            Y1m=bsxfun(@minus,Y1,mean(Y1));
            Y2m=bsxfun(@minus,Y2,mean(Y2));
            Gm=Y1m*Y2m';                    % With mean subtraction
            D.Gm(i,:)=Gm(:)';
        end;
        
        % Prediction of mean and variance without mean subtraction
        G  = U*U';       % This is the predicted mean
        g  = diag(G);    % This is <u_i,u_i>
        pV = VV(:)*P + V(I1(:)).*g(I2(:)) + V(I2(:)).*g(I1(:));
        
        % Prediction of mean and of variance with mean subtraction
        Gm  = Um*Um';       % This is the predicted mean
        gm  = diag(Gm);     % This is <u_i,u_i>
        % Reduction in variance for subtracting the mean this simplifies to
        % V*(1-1/K) for equal variances. The variance of the mean is
        % sum(V)/K^2 and the covariance V./K
        Vm  = V + sum(V)/K^2 - 2.*V./K;
        VVm = Vm*Vm';
        pVm = VVm(:)*P + Vm(I1(:)).*gm(I2(:)) + Vm(I2(:)).*gm(I1(:));
        
        d=[1:16];
        subplot(2,1,1);
        plot(d,mean(D.G) ,'k',d,G(:),'r:',...
            d,mean(D.Gm),'b',d,Gm(:),'g:','LineWidth',3);
        subplot(2,1,2);
        plot(d,var(D.G),'k',d,pV,'r:',...
            d,var(D.Gm),'b',d,pVm,'g:','LineWidth',3);
    case 'matrix_test_folds'    % Predicts the mean and variance of G-hat
        D.P = 50;
        D.K = 3;   % Number of categories
        D.M = 6;   % Number of runs
        D.G   = [1 0.2 0.1 0.2 0.5 0.05 0.1 0.05 0.1];         % Different variances of the constant patterns
        % D.G   = diag([1 0.5 0.2 0.1]);         % Different variances of the constant patterns
        % D.G  = [1 0.3;0.3 1];
        D.eps  = 5;                % Variance for the noise
        Sig  = eye(D.K)*D.eps;
        Sig(1,2)=0.5;
        Sig(2,1)=0.5;
        
        % Make constant pattern that conform exactly to the desired G
        U = normrnd(0,1,D.K,D.P);
        E = (U*U');
        G = reshape(D.G,D.K,D.K);
        A = cholcov(G)';
        U = A*E^(-0.5)*U;
        
        X=kron(ones(D.M,1),eye(D.K));
        part=kron([1:D.M]',ones(D.K,1));
        SSig = kron(eye(D.M),Sig); % And the error terms are assumeed to be independent
        ssig= cholcov(SSig);
        
        for i=1:7000
            Y = repmat(U,D.M,1);
            Y = Y + ssig'*normrnd(0,1,D.K*D.M,D.P);
            [G_hat,SIG] = crossval_estG(Y,X,part);
            GG_hat=Y*Y';
            D.sig(i,1)=mean(diag(SIG));
            D.G_hat(i,:) = G_hat(:)';
            D.GG_hat(i,:) = GG_hat(:)';
        end;
        
        % Prediction of mean
        G  = U*U';       % This is the predicted mean: note division by P
        g  = G(:)';      % stretched out to a vector
        M  = D.M;
        KK  = D.K*D.K;      % Number of possible pairs
        
        % The easiest thing to do is now to define variance-covariance
        % matrices over all conditions and folds and use Equation X to
        % Express this
        GG = kron(ones(D.M),G);    % The constsnt terms are assumed to be preserved across folds
        SSig = kron(eye(D.M),Sig); % And the error terms are assumeed to be independent
        
        % This is the full varaicne-covariance matrix over folds /
        % conditions
        T = size(GG,1);
        i=kron([1:T]',ones(T,1));
        j=kron(ones(T,1),[1:T]');
        
        VV=GG(i,i).*SSig(j,j) + GG(j,j).*SSig(i,i) +GG(i,j).*SSig(j,i)+ GG(j,i).*SSig(i,j) + (SSig(i,i).*SSig(j,j) + SSig(i,j).*SSig(j,i))*D.P;
        % VVb= cov(D.GG_hat);
        % Now use the pairs of cross-folds
        PF=1-eye(D.M);
        for i=1:KK            % Do all the possible sums
            for j=1:KK
                a= zeros(KK,1);
                a(i) = 1;
                Ci=kron(PF,reshape(a,D.K,D.K));
                Ci=Ci(:)./sum(Ci(:));
                a= zeros(KK,1);
                a(j) = 1;
                Cj=kron(PF,reshape(a,D.K,D.K));
                Cj=Cj(:)./sum(Cj(:));
                V(i,j) = Ci(:)'*VV*Cj(:);
            end;
        end;
        
        % You can also get V directly
        V1=crossval_varG(G,Sig,D.M,D.P);
        % Compare with real covariance matrix
        A=cov(D.G_hat);
        
        subplot(3,1,1);
        imagesc(V1);
        subplot(3,1,2);
        imagesc(A);
        subplot(3,1,3);
        scatterplot(V(:),A(:),'identity');
    case 'dist_test_folds'      % Predicts the mean and variance of Distances across folds
        D.P = 50;          % Number of voxels
        D.K = 8;           % Number of categories
        D.M = 6;           % Number of runs
        
        % Make random G
        U   = normrnd(0,1,D.K,D.K);
        G   = U*U'/D.K;
        % G   = diag([1 0.5 0.2 0.1]);         % Different variances of the constant patterns
        
        H   = eye(D.K) - ones(D.K)/D.K;
        G   = H*G*H;
        D.G    = G(:)';
        D.eps  = 1.1;                % Variance for the noise
        Sig  = eye(D.K)*D.eps;
        Sig  = H*Sig*H;
        
        % Make constant pattern that conform exactly to the desired G
        U = mvnrnd_exact(G,D.P);
        
        condvec = kron(ones(D.M,1),[1:D.K]');      % Make a simple design matrix
        part    = kron([1:D.M]',ones(D.K,1));
        SSig    = kron(eye(D.M),Sig); % And the error terms are assumeed to be independent
        ssig    = cholcov(SSig);
        
        C=indicatorMatrix('allpairs',[1:D.K]);
        
        for i=1:7000
            Y = repmat(U,D.M,1);
            Y = Y + ssig'*normrnd(0,1,size(ssig,1),D.P);
            [D_hat,Sig] = rsa_distanceLDC(Y,part,condvec);
            D.D_hat(i,:) = D_hat(:)';
        end;
        
        % Prediction of mean
        G  = U*U';       % This is the predicted mean
        d=diag(C*G*C');
        
        % Get the variances from the var(G(:))-matrix
        % By writing the distances as a linear combination of G:  d=CC'*G(:)
        for c=1:size(C,1);
            A=C(c,:)'*C(c,:);
            CC(:,c)=A(:);
        end;
        VV=crossval_varG(G,D.eps,D.M,D.P);
        V1 = CC'*VV*CC;
        
        % We can also get this directly
        V  = crossval_varDist(d,C,D.eps,D.M,D.P);
        
        % Compare with real covariance
        A=cov(D.D_hat);
        
        subplot(3,1,1);
        imagesc(V);
        subplot(3,1,2);
        imagesc(A);
        subplot(3,1,3);
        scatterplot(V(:),A(:),'identity');
        ylabel('predicted variance');
        xlabel('measured variance');
    case 'dist_test_oneout'     % Predicts the mean and variance of Distances using dependent leave-on-out crossvalidation
        D.P = 50;
        D.K = 5;   % Number of categories
        D.M = 6;   % Number of runs
        % D.G   = [1 0.2 0.1 0.2 0.5 0.05 0.1 0.05 0.1];         % Different variances of the constant patterns
        % Make random G
        % U   = normrnd(0,1,D.K,D.K);
        % G   = U*U'/D.K;
        G   = diag([1 0.5 0.2 0.1 0.3]);         % Different variances of the constant patterns
        G = zeros(5);
        
        H   = eye(D.K) - ones(D.K)/D.K;
        G   = H*G*H;
        % D.G  = [1 0.3;0.3 1];
        D.G    = G(:)';
        D.eps  = 1;                % Variance for the noise
        Sig  = eye(D.K)*D.eps;
        Sig  = H*Sig*H;
        
        % Make constant pattern that conform exactly to the desired G
        U = mvnrnd_exact(G,D.P);
        
        X=kron(ones(D.M,1),eye(D.K));
        part      = kron([1:D.M]',ones(D.K,1));
        condition = kron(ones(D.M,1),[1:D.K]');
        SSig = kron(eye(D.M),Sig); % And the error terms are assumeed to be independent
        ssig= cholcov(SSig);
        
        C=indicatorMatrix('allpairs',[1:D.K]);
        
        for i=1:1000
            Y = repmat(U,D.M,1);
            Y = Y + ssig'*normrnd(0,1,size(ssig,1),D.P);
            for m=1:D.M
                YA = Y(part~=m,:);
                XA = X(part~=m,:);
                YB = Y(part==m,:);
                XB = X(part==m,:);
                dA = C*pinv(XA)*YA;
                dB = C*pinv(XB)*YB;
                DD(m,:)  = sum(dA.*dB,2)';
            end;
            D.d_hat(i,:) = mean(DD)./D.P;           % Estimated distance of all partitions
            D.d11(i,:)   = DD(:,1)';                % Estimated distance only from partition 1 (not devided by P)
        end;
        
        % Prediction of mean
        G  = U*U'./D.P;       % This is the predicted mean
        d=diag(C*G*C');
        
        % Check the variance and covariance of 1 distance across folds:
        dd=d(1).*D.P;
        dSig = C*Sig*C';
        s = dSig(1,1);
        vardPred = dd*(s/(D.M-1) + s) + s*s/(D.M-1)*D.P;
        covdPred = dd*((D.M-2)*s/((D.M-1)^2) + 2*s/(D.M-1))+s/(D.M-1)*s/(D.M-1)*D.P;
        A=cov(D.d11);
        vard = mean(diag(A));
        covd = (sum(A(:))-trace(A))./(D.M*(D.M-1));
        fprintf('Variance of 1st fold:  %2.3f   pred: %2.3f\n',vard,vardPred);
        fprintf('covaiance between folds:  %2.3f   pred: %2.3f\n',covd,covdPred);
        
        % now do this more systematic for the fist difference
        [V,N] = rsa_sigmaLDCraw(eye(D.K*D.M),part,condition);
        V1    = rsa_varianceLDCraw(d,C,V,N,D.M,D.P);
        
        % We can also get this directly
        V2  = rsa.stat.varianceLDC(d,C,D.eps,D.M,D.P);
        
        % Compare with real covariance
        A=cov(D.d_hat);
        
        subplot(3,1,1);
        imagesc(V1);
        subplot(3,1,2);
        imagesc(A);
        subplot(3,1,3);
        scatterplot(V1(:),A(:),'identity');
        keyboard;
    case 'dist_distribution'    % Distribution of distance in a 3x3
        D.P = 50;  % Number of voxels
        D.K = 3;   % Number of categories
        D.M = 3;   % Number of runs
        % D.G   = [1 0.2 0.1 0.2 0.5 0.05 0.1 0.05 0.1];         % Different variances of the constant patterns
        % Make G
        G   = zeros(3);         % Different variances of the constant patterns
        G   = [1 0 0;0 1 0;0 0 1];         % Different variances of the constant patterns
        G   = [3 0 0;0 1 0;0 0 1];         % Different variances of the constant patterns
        G   = [1 0 0;0 1 0;0 0 1];         % Different variances of the constant patterns
        G   = [1 1 0;1 1 0;0 0 1];         % Different variances of the constant patterns
        G   = zeros(3);         % Different variances of the constant patterns
        G   = [1.3 0.7 0;0.7 1 0;0 0 1];         % Different variances of the constant patterns
        G   = zeros(3);         % Different variances of the constant patterns
        G   = [1.4 0 0.5;0 1.3 0;0.5 0 1];         % Different variances of the constant patterns
        H   = eye(D.K) - ones(D.K)/D.K;
        G   = H*G*H;
        
        D.G    = G(:)';
        D.eps  = 1;                % Variance for the noise
        Sig  = eye(D.K)*D.eps;
        Sig  = H*Sig*H;
        
        % Make constant pattern that conform exactly to the desired G
        U = mvnrnd_exact(G,D.P);
        
        conditionVec=kron(ones(D.M,1),[1:D.K]');
        partitionVec=kron([1:D.M]',ones(D.K,1));
        SSig = kron(eye(D.M),Sig); % And the error terms are assumeed to be independent
        ssig= cholcov(SSig);
        
        C=indicatorMatrix('allpairs',[1:D.K]);
        
        for i=1:10000
            Y = repmat(U,D.M,1);
            Y = Y + ssig'*normrnd(0,1,size(ssig,1),D.P);
            D.D_hat(i,:) = rsa.distanceLDC(Y,partitionVec,conditionVec);
        end;
        
        % Prediction of mean
        G  = U*U'/D.P;       % This is the predicted mean: note division by P
        d=diag(C*G*C');
        
        % We can also get this directly
        V  = rsa_varianceLDC(d',C,D.eps,D.M,D.P);
        
        % Plot figure
        xlim=[-1 4];
        ylim=[0 1600];
        catX=[-1:0.05:4];
        c=[0.5 0.5 0.5];
        
        subplot(3,3,[2 3 5 6]);
        scatterplot(D.D_hat(:,1),D.D_hat(:,2),'markertype','.');
        set(gca,'XLim',xlim,'YLim',xlim);
        drawline(0);
        drawline(0,'dir','horz');
        drawline(d(1),'linestyle',':');
        drawline(d(2),'linestyle',':','dir','horz');
        
        subplot(3,3,[1 4]);
        histplot(D.D_hat(:,2),'catX',catX,'linecolor',c,'facecolor',c);
        view(-90,90);
        set(gca,'XLim',xlim,'YLim',ylim);
        drawline(0);drawline(d(2),'linestyle',':');
        
        subplot(3,3,[8 9]);
        histplot(D.D_hat(:,1),'catX',catX,'linecolor',c,'facecolor',c);
        set(gca,'XLim',xlim,'YLim',ylim);
        drawline(0);drawline(d(1),'linestyle',':');
        
        set(gcf,'PaperPosition',[0 0 6 5]);
        wysiwyg;
        
        corr(D.D_hat(:,1),D.D_hat(:,2))
        d
    case 'dist_covariance'      % Variance-covariance matrix for spatially independent data - normal and crossvalidated
        D.P = 30;  % Number of voxels
        D.K = 5;   % Number of categories
        D.M = 5;   % Number of runs
        D.var_a=0; % Amount of signal
        D.eps  = 1;                % Variance for the noise
        D.numSim = 1000;
        G    = zeros(5);
        G(1,1)=1;
        G(2,2)=0.5;
        D.G = G(:)';
        D=rsa.getUserOptions(varargin,D);
        
        G = reshape(D.G,D.K,D.K);
        H   = eye(D.K) - ones(D.K)/D.K;
        G   = H*G*H;
        
        Sig  = eye(D.K)*D.eps;
        Sig  = H*Sig*H;
        
        % Make constant pattern that conform exactly to the desired G
        U = mvnrnd_exact(G*D.var_a,D.P);
        
        conditionVec=kron(ones(D.M,1),[1:D.K]');
        partitionVec=kron([1:D.M]',ones(D.K,1));
        Z=indicatorMatrix('identity',conditionVec); 
        SSig = kron(eye(D.M),Sig); % And the error terms are assumeed to be independent
        ssig= cholcov(SSig);
        
        % make data and compute crossvalidated and non-crossvalidated
        % distances 
        C=indicatorMatrix('allpairs',[1:D.K]);
        for i=1:D.numSim
            Y = repmat(U,D.M,1);
            Y = Y + ssig'*normrnd(0,1,size(ssig,1),D.P);
            T.D_hatc(i,:) = rsa.distanceLDC(Y,partitionVec,conditionVec);
            mY=pinv(Z)*Y;
            T.D_hat(i,:)  = diag(C*mY*mY'*C')./D.P; 
        end;
        
        % Prediction of mean
        G  = U*U'/D.P;       % This is the predicted mean: note division by P

        % Get the predicted mean and variance of the crossvalidated distance
        Xi = C*Sig*C'; 
        De = C*U*U'*C'/D.P; 

        D.E = mean(T.D_hat); 
        D.Ec = mean(T.D_hatc); 
        V = cov(T.D_hat); 
        D.V = V(:)'; 
        Vc = cov(T.D_hatc);
        D.Vc = Vc(:)'; 
        V = cov(ssqrt(T.D_hat)); 
        D.Vs = V(:)'; 
        Vc = cov(ssqrt(T.D_hatc));
        D.Vcs = Vc(:)'; 
        
        D.pEc  = diag(C*G*C')';
        % pVc  = rsa_varianceLDC(G,C,D.eps,D.M,D.P)};
        pVc  = 1./(D.P.^2)*(4*D.P/D.M * De + 2*D.P/(D.M*(D.M-1))*Xi).*Xi;
        D.pVc = pVc(:)'; 
        D.pE  = diag(C*G*C')' + diag(Xi)'./D.M;
        pV  = 1./(D.P.^2)*(4*D.P/D.M * De + 2*D.P/(D.M^2)*Xi).*Xi; 
        D.pV = pV(:)'; 
        varargout={D,T};
    case 'dist_covariance_sim'; 
        Var_a=[0 0.05 0.1 0.15 0.2]*4;
        D=[];
        for i=1:length(Var_a);
            fprintf('%d\n',i);
            [T]=rsa_testVarianceBasic('dist_covariance',...
                'var_a',Var_a(i),'numSim',10000,'P',30,'M',5,'eps',2);
            D=addstruct(D,T);
        end;
        varargout={D}; 
    case 'Fig_histogram'; 
        var_a=0;
        eps = 2; 
        catX = [-1:0.1:1.5];
        vararginoptions(varargin,{'var_a','eps','catX'}); 
        [D,T]=rsa_testVarianceBasic('dist_covariance',...
                'var_a',var_a,'numSim',10000,'P',30,'M',5,'eps',eps);
        subplot(4,1,1); 
        histplot(T.D_hatc(:,1),'catX',catX);
        set(gca,'XLim',[catX(1)-0.2 catX(end)+0.2]); 
        
        subplot(4,1,2); 
        histplot(T.D_hat(:,1),'catX',catX);
        set(gca,'XLim',[catX(1)-0.2 catX(end)+0.2]); 
        
        subplot(4,1,3); 
        histplot(ssqrt(T.D_hatc(:,1)),'catX',catX);
        set(gca,'XLim',[catX(1)-0.2 catX(end)+0.2]); 
        
        subplot(4,1,4); 
        histplot(ssqrt(T.D_hat(:,1)),'catX',catX);
        set(gca,'XLim',[catX(1)-0.2 catX(end)+0.2]); 
        
        set(gcf,'PaperPosition',[2 2 5 7]);
        wysiwyg;

    case 'Fig_variance'
        D=varargin{1};   
        dist  = [1 2 5]; 
        distV = [1 12 45]; 
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
        % Now repeat this for non-cros
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
        xyplot(T.E,T.V,T.trueD,'split',T.crossval,'style_thickline','leg',{'crossvalidated','standard'});
        hold on; 
        xyplot(T.pE,T.pV,T.trueD,'split',T.crossval,'style_thickline','linestyle',':');
        hold off; 
        xlabel('Expected value');
        ylabel('Variance');       
        set(gcf,'PaperPosition',[2 2 4 3]);
        wysiwyg;
        set(gca,'YLim',[0 max(T.pV)+0.01]);
        
    case 'Fig_sqrt'
        D=varargin{1};  
        numDist = size(D.E,2); 
        for i=1:length(D.var_a)
            V=reshape(D.pV(i,:),numDist,numDist); 
            pV=1./(4*sqrt(D.pE(i,:)'*D.pE(i,:))).*V; 
            D.pVs(i,:)=pV(:)';
            V=reshape(D.pVc(i,:),numDist,numDist); 
            EE=D.pEc(i,:)'*D.pEc(i,:); 
            EE(EE<eps)=NaN; 
            pV=V./(4*sqrt(EE)); 
            D.pVcs(i,:)=pV(:)';            
        end; 
        % Order the data into a new data frame 
        dist  = [1 2 5]; 
        distV = [1 12 45]; 
        x=sqrt(D.pEc(:,dist));
        T.trueD = x(:); 
        T.pE=x(:); 
        x=ssqrt(D.Ec(:,dist));
        T.E=x(:); 
        x=D.pVcs(:,distV);
        T.pV=x(:); 
        x=D.Vcs(:,distV); 
        T.V=x(:); 
        T.crossval = ones(length(T.trueD),1); 
        % Now repeat this for non-cros
        S=T; 
        x=sqrt(D.pE(:,dist));
        S.pE=x(:); 
        x=sqrt(D.E(:,dist));
        S.E=x(:); 
        x=D.pVs(:,distV);
        S.pV=x(:); 
        x=D.Vs(:,distV); 
        S.V=x(:); 
        S.crossval = ones(length(T.trueD),1)*2; 
        T=addstruct(T,S); 
        xyplot(T.E,T.V,T.trueD,'split',T.crossval,'style_thickline','leg',{'crossvalidated','standard'});
        hold on; 
        xyplot(T.pE,T.pV,T.trueD,'split',T.crossval,'style_thickline','linestyle',':');
        hold off; 
        xlabel('Expected value');
        ylabel('Variance');       
        set(gcf,'PaperPosition',[2 2 4 3]);
        wysiwyg;
        set(gca,'YLim',[0 max(T.V)+0.02]);
        
    case 'Figure1'
        label={'12','13','14','15','23','24','25','34','35','45'};
        numDist = size(D.E,2); 
        lineplot(D.var_a,[D.V(:,1) D.Vs(:,1)],'leg','auto','style_thickline'); 
        subplot
        
        set(gcf,'PaperPosition',[2 2 5.5 12]);
        wysiwyg;
        
        colormap(hot);
        clim=[0 0.03];
        subplot(4,2,1);
        imagesc(reshape(T.V(1,:),numDist,numDist),clim);
        set(gca,'YTick',[1:10],'YTickLabel',label,'Box','off');
        subplot(4,2,2);
        imagesc(reshape(T.V(end,:),numDist,numDist),clim);
        set(gca,'YTick',[1:10],'YTickLabel',label,'Box','off');
        
        subplot(4,2,3);
        imagesc(reshape(T.pV(1,:),numDist,numDist),clim);
        set(gca,'YTick',[1:10],'YTickLabel',label,'Box','off');
        subplot(4,2,4);
        imagesc(reshape(T.pV(end,:),numDist,numDist),clim);
        set(gca,'YTick',[1:10],'YTickLabel',label,'Box','off');
        
        subplot(4,2,5);
        p=[0.01:0.02:0.99];
        h=quantile(D_hat{1}(:,1),p);
        q=norminv(p,D.d(1,1),sqrt(D.pVar12(1)));
        scatterplot(q',h','identity');
        set(gca,'XLim',[-0.31 0.31],'YLim',[-0.31 0.31],'XTick',[-0.3:0.1:0.3],'YTick',[-0.3:0.1:0.3]);
        
        subplot(4,2,6);
        h=quantile(D_hat{end}(:,1),p);
        q=norminv(p,D.d(end,1),sqrt(D.pVar12(end)));
        scatterplot(q',h','identity');
        set(gca,'XLim',[-0.15 0.8],'YLim',[-0.15 0.8],'XTick',[0:0.2:0.8],'YTick',[0:0.2:0.8]);
        
        subplot(4,2,[7 8]);
        lineplot(D.var_a,[D.var12 D.var23 D.var45 D.cov12_23 D.cov12_45],'style_thickline','leg',{'v1','v5','v10','c1-5','c1-10'});
        hold on;
        lineplot(D.var_a,[D.pVar12 D.pVar23 D.pVar45 D.pCov12_23 D.pCov12_45],'style_thickline','leg',{'v1','v5','v10','c1-5','c1-10'},'linestyle',':');
        hold off;
        
        varargout={D};       
    case 'voxel_dependence'
        D.K = 2;           % Number of categories
        D.M = 2;           % Number of runs
        D.s_a     = 0.001; % Width of spatial noise kernel
        D.s_e     = 0.001; % Width of spatial signal kernel
        D.var_a   = 0;     % Signal variance
        D.var_e   = 1;     % Noise variance
        D.numSim  = 5000;  % Number of simulations
        D=rsa.getUserOptions(varargin,D);
        D.N = D.K*D.M;     % Number of trials
        
        C.dist=rsa_testVarianceBasic('3d_sphere',3);
        D.P=size(C.dist,1);
        Z=kron(ones(D.M,1),eye(D.K));
        
        C.trueG=eye(D.K).*D.var_a;
        
        C.run  = kron([1:D.M],ones(1,D.K));
        C.cond = kron(ones(1,D.M),[1:D.K]);
        
        % make spatial cholinsky matrices to generate spatially dependent
        % random variables
        C.SigA = exp(-C.dist.^2/D.s_a);
        C.SigE = exp(-C.dist.^2/D.s_e);
        C.cholSigA = cholcov(C.SigA);
        C.cholSigE = cholcov(C.SigE);
        C.cholG  = cholcov(C.trueG);
        C.cholST = kron(C.cholG,C.cholSigA);
        
        % Contrast matrix
        Con=indicatorMatrix('allpairs',[1:D.K]);
        
        % Generate the true patterns
        trueU=randn(1,size(C.cholST,1))*C.cholST;
        trueU=reshape(trueU,D.P,size(Z,2))';
        
        % Generate data with spatially correlated noise
        for i=1:D.numSim
            err=randn(D.N,D.P)*C.cholSigE*sqrt(D.var_e);
            y=Z*trueU+err;
            V.dist(i,:)=rsa.distanceLDC(y,C.run',C.cond');
        end;
        
        % Prediction of mean
        G  = trueU*trueU'./D.P;
        D.dist_true = diag(Con*G*Con');
        D.dist      = mean(V.dist);
        D.var       = var(V.dist);
        D.var_pred  = rsa_varianceLDC(D.dist_true,Con,D.var_e,D.K,D.P);
        G  = squareform(-0.5 * D.dist_true);
        dG = Con*G*Con';                         % Predicted variance-covariance of the signals
        Sig = eye(D.K).*D.var_e;
        dSig   = Con*Sig*Con';
        nFolds = D.M*(D.M-1);
        D.var_pred2=(4*(dG.*dSig)/D.M+2*(dSig.*dSig)/nFolds)*(sum(sum(C.SigE.*C.SigE))/D.P^2);
        varargout={D};
    case 'voxel_dependence_run'
        SigE=[0.1 0.5 1 2 3];
        varA=[0]; %  0.05 0.1 0.15];
        T=[];
        for j=1:length(varA);
            for i=1:length(SigE);
                D=rsa_testVarianceBasic('voxel_dependence','s_e',SigE(i),'s_a',1,'var_a',varA(j));
                T=addstruct(T,D);
            end;
        end;
        
        subplot(2,1,1);
        lineplot([T.s_e],[T.var],'split',T.var_a,'style_thickline');
        hold on;
        lineplot([T.s_e],[T.var_pred],'split',T.var_a,'style_thickline','linestyle',':');
        hold off;
        
        subplot(2,1,2);
        lineplot([T.s_e],[T.var],'split',T.var_a,'style_thickline');
        hold on;
        lineplot([T.s_e],[T.var_pred2],'split',T.var_a,'style_thickline','linestyle',':');
        hold off;
        varargout={T};
    case 'voxel_dependence_det' % Measure of multivariate dependence and it's distribution (see Reddon et al. 1985)
        N=200;
        T=[];
        for numVox = [5 10 50 100 150]
            rho = -N+(2*numVox+5)/6;
            for n=1:100
                X=normrnd(0,1,N,numVox);
                R=corr(X);
                D.logdetR(n,1) = logdet(R);
                D.rho(n,1) = rho;
                D.numVox(n,1) = numVox;
                D.m(n,1)= numVox*(numVox-1)/2;
                D.dep(n,1) = 1-exp(logdet(R)/numVox); % Dependence measure of Pena et al.
            end;
            D.depPred = 1-exp(D.m./D.rho./D.numVox);
            T=addstruct(T,D);
        end;
        varargout={T};
    case 'wishart1'             % Checks the mean and variance of the central and non-central wishart distr. (df=1)
        V=[1 0.3;0.3 1];         % Variance-covariance matrix of the variables
        mu=[0.3 0.2]';
        A=cholcov(V);
        X=randn(100000,2)*A;     % Generate normal random variables
        X=bsxfun(@plus,X,mu');
        Y=[X.*X X(:,1).*X(:,2)]; % Take the products
        % Predicted mean and variance for the central wishart
        % v = diag(V);
        % pM = V;
        % pV = V.*V + v*v';
        % MeanCentral1 = [pM(1,1) pM(2,2) pM(1,2)];
        % VarCentral1  = [pV(1,1) pV(2,2) pV(1,2)];
        % Predicted mean and variance for the non-central wishart
        pM = mu*mu'+V;
        pV = 2*(mu*mu').*V+(mu.*mu)*v'+v*(mu.*mu)'+V.*V+v*v';
        MeanNonC1 = [pM(1,1) pM(2,2) pM(1,2)];
        VarNonC1  = [pV(1,1) pV(2,2) pV(1,2)];
        fprintf('Measured mean : %2.3f %2.3f %2.3f\n',mean(Y));
        fprintf('Predicted mean: %2.3f %2.3f %2.3f\n',MeanNonC1);
        fprintf('Measured var  : %2.3f %2.3f %2.3f\n',var(Y));
        fprintf('Predicted var : %2.3f %2.3f %2.3f\n',VarNonC1);
    case 'wishartN'             % Checks the mean and variance of the central and non-central wishart distr. (df=N)
        V=[1.4 0.3;0.3 1];         % Variance-covariance matrix of the variables
        N=10;
        samples = 10000;
        mu=[0.3 0.2]';
        M=kron(ones(1,N),mu);  % this has the Means the same
        M=normrnd(0,2,2,N);    % This has the means different across voxels
        A=cholcov(V);
        for i=1:samples
            X=A'*randn(2,N)+M;     % Generate normal random variables
            s=X*X';
            S(i,1)=s(1,1);
            S(i,2)=s(2,2);
            S(i,3)=s(1,2);
        end;
        % Predicted mean and variance for the non-central wishart
        v = diag(V);
        pM = M*M'+N*V;
        sM=sum(M.*M,2);
        pV = 2*(M*M').*V+sM*v'+v*sM'+N*(V.*V+v*v');
        MeanNonC1 = [pM(1,1) pM(2,2) pM(1,2)];
        VarNonC1  = [pV(1,1) pV(2,2) pV(1,2)];
        fprintf('Measured mean : %2.3f %2.3f %2.3f\n',mean(S));
        fprintf('Predicted mean: %2.3f %2.3f %2.3f\n',MeanNonC1);
        fprintf('Measured var  : %2.3f %2.3f %2.3f\n',var(S));
        fprintf('Predicted var : %2.3f %2.3f %2.3f\n',VarNonC1);
    case 'wishartS'             % Checks the mean and variance of the non-central wishart distr (df=N), if the observations are not independent
        P   = 2;  % Number of variates (rows)
        N   = 10; % Number of observations (columns)
        V   = [1.4 0.3;0.3 1];  % Variance-covariance matrix of the variables
        width = 4;              % Generate correltation matrix of
        Sig = bsxfun(@minus,[1:N],[1:N]');
        Sig = exp(-Sig.^2/width)*2;
        samples = 50000;
        mu=[1 2]';
        M=kron(ones(1,N),mu);  % this has the Means the same
        M=normrnd(0,2,P,N);    % This has the means different across voxels
        % Perpare the generation of the matrix multivariate normal N(M,V,S);
        VS = kron(Sig,V);
        cVS=cholcov(VS);
        for i=1:samples
            X = cVS'*randn(N*P,1);
            X = reshape(X,P,N);
            X = X+M;     % Generate normal random variables
            s = X*X';    % Claculate outer product and store in vectorised form
            S(i,1)=s(1,1);
            S(i,2)=s(2,2);
            S(i,3)=s(1,2);
        end;
        % Predicted mean and variance for the non-central wishart
        v = diag(V);
        % 
        KK=M*Sig*M'; 
        sK=diag(KK);
        pM = M*M'+trace(Sig)*V;
        pV = 2*(KK).*V+sK*v'+v*sK'+trace(Sig*Sig)*(V.*V+v*v');
        % sM=sum(M.*M,2);
        % pV = 2*(M*M').*V+sM*v'+v*sM'+N*V.*V+N*v*v';
        MeanNonC1 = [pM(1,1) pM(2,2) pM(1,2)];
        VarNonC1  = [pV(1,1) pV(2,2) pV(1,2)];
        fprintf('Measured mean : %2.3f %2.3f %2.3f\n',mean(S));
        fprintf('Predicted mean: %2.3f %2.3f %2.3f\n',MeanNonC1);
        fprintf('Measured var  : %2.3f %2.3f %2.3f\n',var(S));
        fprintf('Predicted var : %2.3f %2.3f %2.3f\n',VarNonC1);
    case 'covarianceD'     % Checks the predicted (co-)variance matrix of the diagonal of XX' for non-zero mean and dependent voxels 
        P   = 3;  % Number of variates (rows)
        N   = 10; % Number of observations (columns)
        V   = [1.4 0.3 0;0.3 1 0.5;0 0.5 1];  % Variance-covariance matrix of the variables
        width = 4;              % Generate correltation matrix of
        Sig = bsxfun(@minus,[1:N],[1:N]');
        Sig = exp(-Sig.^2/width)*2;
        samples = 50000;
        mu=[1 2]';
        M=kron(ones(1,N),mu);  % this has the Means the same
        M=normrnd(0,2,P,N);    % This has the means different across voxels
        % Perpare the generation of the matrix multivariate normal N(M,V,S);
        VS = kron(Sig,V);
        cVS=cholcov(VS);
        for i=1:samples
            X = cVS'*randn(N*P,1);
            X = reshape(X,P,N);
            X = X+M;     % Generate normal random variables
            D(i,:)=diag(X*X');  % Here we care only about the diagonal of XX'
        end;
        % Predicted mean and co-variance for diagonal elements 
        v = diag(V);
        m=diag(M*Sig*M'); 
        pM = diag(M*M')+trace(Sig)*diag(V);
        pC = 4*M*Sig*M'.*V+2*trace(Sig*Sig)*(V.*V);
        C = cov(D);
        subplot(1,3,1); 
        plot([pM mean(D)']); 
        subplot(1,3,[2:3]); 
        imagesc([pC C]); 
end;
