function varargout = sh1_rsa(what,varargin)
%% Analysis script for the hierarchical sequence learning project
%  mostly based on: sl4_imana.m  by ?
%                   cpd1_imana.m by Naveed Ejaz
%					slb1_imana.m by Atsushi Yokoi
% July 2014 Atsushi Yokoi


%% Directory definitions
volumenames = {'/Volumes/HD-PEBU2/Windows/myBackup/ICN/',...
    '/Volumes/Macintosh HD/Users/Atsushi',...
    '/Volumes/DATA/',...
    '/Users/Atsushi/',...
    'F:/Windows/myBackup/Dropbox/work/UK/Research',...
    'H:/Dropbox/work/UK/Research/',...
    '/Users/joern/Projects'};%

baseDir         = '/Users/joern/Projects/sh1';   % For Joern
% baseDir         = '/Volumes/DATA/MotorControl/data/SequenceLearning/sh1';   % For AY

fieldmapsDir    = fullfile(baseDir, 'fieldmaps');
behaviourDir    = fullfile(baseDir, 'data');
analyzeDir 		= fullfile(baseDir, 'analyze');
anatomicalDir   = fullfile(baseDir, 'anatomicals');
imagingRawDir   = fullfile(baseDir, 'imaging_data_raw');
imagingDir      = fullfile(baseDir, 'imaging_data');
freesurferDir   = fullfile(baseDir, 'surfaceFreesurfer');
caretDir        = fullfile(baseDir, 'surfaceCaret');
regDir          = fullfile(baseDir, 'RegionOfInterest');
analysisDir     = fullfile(baseDir, 'analysis');
figDir          = fullfile(baseDir, 'figures');
cerebDir        = fullfile(anatomicalDir,'SUIT');
BGDir           = fullfile(anatomicalDir,'basal_ganglia');

%% Subject specific information (file names etc.)
subj_name       = {'p01','t01','t02','s03','s04','s07',... % s05
    's09','s10','s14','s15', 's17','s19','s22','s23','s25'};
subj_name_behav = {'sh1_p01_scan','sh1_t01_scan','sh1_t02_scan',...
    'sh1_s03_scan','sh1_s04_scan','sh1_s07_scan',...
    'sh1_s09_scan','sh1_s10_scan','sh1_s14_scan','sh1_s15_scan',...
    'sh1_s17_scan','sh1_s19_scan','sh1_s22_scan','sh1_s23_scan','sh1_s25_scan'};

%% Experiment specific parameters (TR info etc.)
use3D           = 1;            % in certain envinronment, using 3D image makes analysis faster
startTR         = 6;            % ceil(startSlice(trial=1)/32) + 1
endTR           = 134;%124      % ceil(startSlice(trial=end)/32 + complete(trial=end)/2720)
additionalTR    = 1;            % arbitrary: to cover the slow HRF responce at last trial
nTR             = endTR+additionalTR-startTR;
hrf_nTR         = 0;
nSequence       = 8;           % number of sequence

%% Subject specific run/session information (for in case we need to pull subject out during experiment)
nRun = 9;% 10
run             = {{'1','2','3','4','5','6','7','8','9'}};
sess            = {'sess1'};

subj_sess       = {sess, sess, sess, sess, sess, sess,... % set A
    {'sess1','sess2'},sess, sess, sess, sess,{'sess1','sess2'},... % set B
    sess,sess,sess}; % set A (s22),setB (s23,s25)
subj_runs       = {{{'1','2','3','4','5','6','7','8','9','10'}},... % p01
    run, run, run, run, run,...
    {{'1','2','3','4','5','6'},{'7','8','9'}},...% s09
    run, run, run, run,...
    {{'1','2','3','4','5','6','7'},{'8','9'}},... % s19
    run,run,... % s22,23
    {{'1','3','4','5','6','7','8','9','11'}}}; % s25
chunk_sets      = {'A','A','A','A','A','A',...
    'B','B','B','B','B','B',...
    'A','B','B'};
chunk_set       = [1 1 1 1 1 1 2 2 2 2 2 2 1 2 2];

%% 1st-level GLM, Freesurfer & ROI parameters
atlasA      = {'x'};
atlasname   = {'fsaverage_sym'};
glmName     = {'GLM_firstlevel_1','GLM_firstlevel_2','GLM_firstlevel_3','GLM_firstlevel_4'};
radius      = 12;
hem         = {'lh','rh'};
hemName     = {'LeftHem','RightHem'};
hemName_s   = {'L','R'};
numregions_surf = 8;
numregions_BG   = 4;
numregions      = numregions_surf+numregions_BG;
regSide     = [1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 2 2 2 2];
regType     = [1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 9 10 11 12 9 10 11 12];
regname         = {'S1','M1','PMd','PMv','SMA','V12','SPLa','SPLp','CaudateN' 'Pallidum', 'Putamen' 'Thalamus'};
regname_surf    = {'S1','M1','PMd','PMv','SMA','V12','SPLa','SPLp'};
regname_BG      = {'CaudateN' 'Pallidum', 'Putamen' 'Thalamus'};

%% Analysis parameters (chunk and sequence definitions)
% Chunk & Seq definition
% chunk
Chunk        = {[1 3];     % 1
    [5 2 4];
    [2 3 2];
    [5 1 4]
    [3 5];     % 5
    [4 2 1];
    [2 5 2];
    [1 4 3]
    };

Chunk_B      = {[1 4];
    [3 2 5];
    [2 4 2];
    [1 3 5];
    [5 1];
    [4 2 3];
    [2 5 2];
    [4 1 3];
    };

% chunk length
Clength     = [2 3 3 3 2 3 3 3];

% sequence def for real exp
Sequence         = [1 2 3 4;...
    4 3 2 1;...
    5 6 7 8;...
    8 7 6 5;...
    1 2 7 8;...
    8 7 2 1;...
    5 6 3 4;...
    4 3 6 5];

Sequence_B =   [4 3 2 1;
    5 6 7 8;
    8 7 6 5;
    1 2 3 4;
    4 3 6 5;
    1 2 7 8;
    8 7 2 1;
    5 6 3 4];
Contrasts      = {'seq1','seq2','seq3','seq4','Seq5','seq6','seq7','seq8',...
    'allMovement', 'F-contrast', '2333','3332'};

varargout = {};

%% Graph parameters
CAT.facecolor   = {[0 0 0.7],[0.3 0.3 1],[0.7 0 0],[1 0.6 0.6]}; %[1 0.1 0.1],
CAT.fillcolor   = CAT.facecolor;
CAT.linecolor   = CAT.facecolor;
CAT.markertype  = {'o'};
CAT.markersize  = {5};
CAT.markercolor = CAT.fillcolor;
CAT.markerfill  = {[1 1 1]};
CAT.gapwidth    = {[0.2 0.2 0.2 0.2]};
CAT.linewidth   = {1.5};
CAT.errorwidth  = {1.5};
CAT.errorcolor  = {[0 0 0]};

omega_names = {'One-digit','Two-digit','Chunk','Sequence'};

%% Main operation
switch(what)
    case 'ROI_distraw'  % Extracts ROI data and calculates distances
        T=[];
        sn = varargin{1};
        chunkset = 'A';
        glm = 3;
        ROI = 'all';
        fname 	= 'movement';       % using metric file which compares any movement against rest
        regions = 1:8;
        hemi = 1:2;
        vararginoptions({varargin{2:end}},{'chunkset','glm','ROI','fname','regions','hemi'});
        
        Data=[];
        for s=sn
            fprintf('subj = %d\n',s)
            
            % load SPM.mat
            glmDirSubj=fullfile(baseDir,glmName{glm},subj_name{s});
            load(fullfile(glmDirSubj,'SPM.mat'));
            SPM = spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data',subj_name{s},'sess1'));
            
            % choose ROI
            switch (ROI)
                case 'func'
                    load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_reg800_%s_1.mat',fname)]));
                case 'all'
                    load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_regAll_%s_1.mat',fname)]));
            end
            
            % Get design
            E = load(fullfile(glmDirSubj,'SPM_info.mat'));
            E = getrow(E,E.run<11);
            
            % Now get prewhitened betas from these regions
            NRegions = numel(regions)*numel(hemi);
            nRegions = numel(R);
            for reg = regions
                for hem = hemi
                    if reg<9
                        roi = numregions_surf*(hem-1)+reg;
                    elseif reg>=9
                        roi = 2*numregions_surf + numregions_BG*(hem-1)+(reg-8);
                    end
                    if (nRegions<NRegions)
                        warning('number of ROIs doesn''t much!');
                        break;
                    end
                    if (~isempty(R{roi}))
                        fprintf('extracting: %s ...',R{roi}.name);
                        
                        % get data
                        data = region_getdata(SPM.xY.VY,R{roi});
                        numVox = size(R{roi}.data,1);
                        % cut out NaN or all-zero voxels
                        idx = sum(data.*data,1)==0;
                        if sum(idx)==0
                            fprintf('... done (%d voxels)\n',numVox);
                        else
                            fprintf('... done (%d voxels, %d were discarded since containing no data)\n',P,sum(idx));
                        end
                        data = data(:,~idx);
                        numVox = size(data,2);
                        
                        % get the distances
                        [S.RDM,Sw,S.effVox,S.trSS]=rsa_distanceLDCsepPerm(data,SPM,E.seqType);
                        S.Sigma  = Sw(:)';
                        S.region = reg;
                        S.subj  = s;
                        S.hemis = hem;
                        S.numVox = numVox;
                        
                        T=addstruct(T,S);
                    end;
                end;
            end;
        end;
        varargout={T};
        save(fullfile(regDir,'distances_sepPerm.mat'),'-struct','T');
    case 'searchlight_distraw'         % **new** Calc G and pattern distance (prewhiten data using covariance matrix estimated from raw time series data)
        sn              = varargin{1};
        glm=3; 
        numSigma 		= nSequence^2; % length of G matrix
        numDist= nSequence*(nSequence-1)/2; % length of distance
        suffix=[]; 
        
        home = cd;
        for s = sn
            glmDirSubj=fullfile(baseDir,glmName{glm},subj_name{s});
            cd(glmDirSubj);
            % load SPM.mat
            
            nii_out = {};
            % Define output files
            % dist -> 4D
            for elem = 1:numDist
                nii_out{end+1} = fullfile(glmDirSubj,[subj_name{s},'_dist_raw', suffix, '.nii,', num2str(elem)]);
            end
            % sigma -> 4D
            for elem = 1:numSigma
                nii_out{end+1} = fullfile(glmDirSubj,[subj_name{s},'_sigma_raw', suffix, '.nii,', num2str(elem)]);
            end
            
            % effVox
            nii_out{end+1} = fullfile(glmDirSubj,[subj_name{s},'effVox', suffix, '.nii']);
            
            % load search light definition, etc.
            load(fullfile(glmDirSubj,'SPM.mat'));
            SPM = spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data',subj_name{s},'sess1'));
            T   = load(fullfile(glmDirSubj,'SPM_info.mat'));
            
            %SPM = spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data',subj_name{s},'sess1'));
            SL  = load(fullfile(glmDirSubj,['vol_roi_160vox.mat']));
            idx = find(T.regType==1);
            
            params = {SPM,T.seqNum(idx,1)};
            
            % run searchlight analysis
            lmva_spm(SL,SPM.xY.P,nii_out,@sh1_calcDistancePerm,'params',params,'isNP',1);
        end % s
        
    case 'ROI_reliability_btw' % Calculate between subject reliability of RDMs
        
        C = [];
        for glm = [1 3]
            % load data set
            T = load(fullfile(regDir,sprintf('distances_sepPerm.glm%d.mat',glm)));
            
            chunksets = {'A','B'};
            for reg = [1:8]
                for h = [1 2]
                    for chunkset = [1 2]
                        % calc within chunkset value
                        D = getrow(T,T.region==reg&T.hemis==h&T.chunkset==chunksets{chunkset});
                        D_ = getrow(T,T.region==reg&T.hemis==h&T.chunkset~=chunksets{chunkset});
                        
                        rdm_all = nanmean(D.RDM,1); % all rdm for the same chunkset
                        rdm_all_ = nanmean(D_.RDM,1); % all rdm for other chunkset
                        subj = unique(D.subj);
                        for s=1:numel(subj)
                            idx = D.subj==subj(s);
                            rdm = D.RDM(idx,:);
                            rdm_loo = nanmean(D.RDM(~idx,:),1); % leave-one-subject-out rdm
                            
                            R.within_u  = corr(rdm',rdm_all','type','Pearson'); % upper bound
                            R.within_l  = corr(rdm',rdm_loo','type','Pearson'); % lower bound
                            
                            % we must choose only the distance between physically identical sequences here
                            tmprdm      = getCommon(rdm,[3 5 7]);
                            tmprdm_     = getCommon(rdm_all_,[3 5 7]);
                            R.across    = corr(tmprdm',tmprdm_','type','Pearson'); % across chunk set corr
                            
                            tmprdm      = getCommon(rdm,[3 5 7]);
                            tmprdm_     = getCommon(rdm_loo,[3 5 7]);
                            R.within_l_reduced  = corr(tmprdm',tmprdm_','type','Pearson'); % lower bound (calculated with reduced RDM)
                            
                            R.subj      = subj(s);
                            R.chunkset  = chunkset;
                            R.region    = reg;
                            R.hemis     = h;
                            R.glm       = glm;
                            
                            C = addstruct(C,R);
                        end
                    end
                end
            end
        end
        
        for glm=[1 3]
            figure('name',sprintf('reliability (glm=%d)',glm));
            D = getrow(C,C.glm==glm);
            sh1_rsa('plot_reliability_btw',D);
        end
        
        varargout = {C};
        
    case 'fit_model_EB_lin'
        T=load(fullfile(regDir,'distances_sepPerm.mat'));
        regions=unique(T.region);
        for r=regions'
            for h=[1 2]
                indx = find(T.region==r & T.hemis==h);
                D=getrow(T,indx);
                for cs = 1:2
                    Xm(1,:,cs)=sh1_getRDMtempw(cs,1,-2,'sqEuclidean'); % One digit
                    Xm(2,:,cs)=sh1_getRDMtempw(cs,2,1,'sqEuclidean'); % two digit
                    Xm(3,:,cs)=sh1_getRDMtempw(cs,4,0,'sqEuclidean'); % chunk
                    Xm(4,:,cs)=sh1_getRDMtempw(cs,6,1,'sqEuclidean'); % sequence
                    Xm(:,:,cs)=bsxfun(@rdivide,Xm(:,:,cs),sqrt(sum(Xm(:,:,cs).^2,2)));
                end;
                for s=1:length(D.subj)
                    Sigma(:,:,s)=reshape(D.Sigma(s,:),8,8);
                end;
                Model.X=Xm(:,:,chunk_set(D.subj));
                [omega(indx,:),logEvidence(indx,:),lt]=rsa_fitModelHierarchEB(Model,D.RDM,Sigma,9,D.effVox);
                logtheta(indx,:)=repmat(lt,length(indx),1);
            end;
        end;
        T.omega=omega;
        T.logtheta=logtheta;
        varargout={T};
    case 'fit_model_EB_nonlin'
        glm = 1;
        useolddata = 1;
        vararginoptions(varargin(:),{'glm','useolddata'})
        
        % load data set
        if useolddata
            %T       = load(fullfile(regDir,'distances_glm1.mat'));
            T       = load(fullfile(regDir,'allRDMs_ROI_all_glm3.mat'));
            %T       = T.RDM1;
            T.RDM   = bsxfun(@rdivide,T.RDM,T.numVox);
        else
            T       = load(fullfile(regDir,sprintf('distances_sepPerm.glm%d.mat',glm)));
        end
        
        %regions=unique(T.region);
        regions=[1:8]';
        mRDM = [];
        for r=regions'
            for h=[1 2]
                indx = find(T.region==r & T.hemis==h);
                D=getrow(T,indx);
                for s=1:length(D.subj)
                    Sigma(:,:,s)=reshape(D.Sigma(s,:),8,8);
                end;
                
                % organize model RDMs
                Model.fcn = @sh1_getRDMmodelTau1;
                Model.constantParams = {chunk_set(D.subj),...       % Model terms
                    [1 2 4 6],...               % Chunk set
                    'sqEuclidean'...            % Distance term
                    };
                Model.numComp   = numel(Model.constantParams{2}); % Number of linear components
                Model.numPrior  = numel(Model.constantParams{2}); % Number of prior variances on parameters
                Model.numNonlin = numel(Model.constantParams{2}); % Number of nonlinear parameters
                Model.nonlinP0  = [-10 1 10 0];%zeros(size(Model.constantParams{2})); % Starting value of nonlinear(mixing) parameters
                
                % estimate
                [w,logEvidence(indx,:),lt] = rsa_fitModelHierarchEB(Model,D.RDM,Sigma,9,D.effVox);
                
                %
                logtheta(indx,:)=repmat(lt,length(indx),1);
                
                % normalize omega with mean of each RDM
                M = feval(Model.fcn,lt(Model.numPrior+1:Model.numPrior+Model.numNonlin),...
                    Model.constantParams{:});
                omega(indx,:) = w.*permute(mean(M.RDM,2),[3 1 2]);
                
                % save resultant RDMs with estimated tau
                A = feval(Model.fcn,lt(Model.numPrior+1:Model.numPrior+Model.numNonlin),...
                    1,Model.constantParams{2:end}); % chunk set A
                B = feval(Model.fcn,lt(Model.numPrior+1:Model.numPrior+Model.numNonlin),...
                    2,Model.constantParams{2:end}); % chunk set B
                mrdm.RDM    = [A.RDM;B.RDM];
                mrdm.set    = kron([1 2]',ones(Model.numComp,1));
                mrdm.comp   = repmat(Model.constantParams{2}',2,1);
                mrdm.name   = {A.name{:};B.name{:}};
                mrdm.region = repmat(r,size(mrdm.set));
                mrdm.hemis  = repmat(h,size(mrdm.set));
                mRDM        = addstruct(mRDM,mrdm);
                %
            end;
        end;
        T.omega=omega;
        T.logtheta=logtheta;
        T.logEvidence=logEvidence;
        %T.logEvidenceSplit = logEvidenceSplit;
        
        sh1_rsa('plot_omega_prior',T);
        sh1_rsa('plot_prior',T);
        sh1_rsa('plot_tau',T);
        varargout={T,mRDM};
    case 'fit_model_EB_nonlin_poolROI'
        glm = 3;
        useolddata = 0;
        modelTerms = [1 2 4 6];
        Niter       = 50; % iteration from different initial parameter values
        Ntake       = 50; % how many best models you choose
        region      = [3 4];
        vararginoptions(varargin(:),{'glm','useolddata','modelTerms','Niter','Ntake','region'})
        
        TT=[];
        % load data set
        if useolddata
            %T       = load(fullfile(regDir,'distances_glm1.mat'));
            T       = load(fullfile(regDir,'allRDMs_ROI_all_glm3.mat'));
            %T       = T.RDM1;
            T.RDM   = bsxfun(@rdivide,T.RDM,T.numVox);
        else
            T       = load(fullfile(regDir,sprintf('distances_sepPerm.glm%d.mat',glm)));
        end
        
        % pool all ROI data to estimate single tau parameter for whole brain
        mRDM = [];
        indx = find(ismember(T.region,region)&T.hemis==1);%true(size(T.subj));%find(T.region==1&T.hemis==1);
        D=getrow(T,indx);
        for s=1:length(D.subj)
            Sigma(:,:,s)=reshape(D.Sigma(s,:),8,8);
        end;
        
        for iter = 1:Niter
            iter
            initialguess = rand(size(modelTerms))*3-1.5;
            initialguess(end) = 0
            
            % organize model RDMs
            Model.fcn = @sh1_getRDMmodelTau1;
            Model.constantParams = {chunk_set(D.subj),...       % Model terms
                modelTerms,...               % Chunk set
                'sqEuclidean'...            % Distance term
                };
            Model.numComp   = numel(Model.constantParams{2}); % Number of linear components
            Model.numPrior  = numel(Model.constantParams{2}); % Number of prior variances on parameters
            Model.numNonlin = numel(Model.constantParams{2}); % Number of nonlinear parameters
            Model.nonlinP0  = initialguess; % Starting value of nonlinear(mixing) parameters
            
            % estimate
            [w(:,:,iter),levid(:,iter),lt(iter,:)] = rsa_fitModelHierarchEB(Model,D.RDM,Sigma,9,D.effVox);
            
            mevid = mean(levid,1);
        end
        [~,indices]= sort(mevid,2,'descend');
        maxindices = indices(1:Ntake);%find(mevid==max(mevid),Ntake,'first');
        
        for m=1:length(maxindices)
            
            maxidx=maxindices(m);
            % choose
            logtheta(:,:)    = repmat(lt(maxidx,:),length(indx),1);
            ww                  = squeeze(w(:,:,maxidx));
            logEvidence(:,:) = levid(:,maxidx);
            
            % normalize omega with mean of each RDM
            M = feval(Model.fcn,lt(Model.numPrior+1:Model.numPrior+Model.numNonlin),...
                Model.constantParams{:});
            normX = sqrt(mean(M.RDM.^2,2));
            omega(:,:)   = ww;
            omegaN(:,:)  = ww.*permute(normX,[3 1 2]);%bsxfun(@times,ww,normX');
            
            % save resultant RDMs with estimated tau
            %         A = feval(Model.fcn,lt(Model.numPrior+1:Model.numPrior+Model.numNonlin),...
            %             1,Model.constantParams{2:end}); % chunk set A
            %         B = feval(Model.fcn,lt(Model.numPrior+1:Model.numPrior+Model.numNonlin),...
            %             2,Model.constantParams{2:end}); % chunk set B
            %         mrdm.RDM    = [A.RDM;B.RDM];
            %         mrdm.set    = kron([1 2]',ones(Model.numComp,1));
            %         mrdm.comp   = repmat(Model.constantParams{2}',2,1);
            %         mrdm.name   = {A.name{:};B.name{:}};
            %         mrdm.region = repmat(r,size(mrdm.set));
            %         mrdm.hemis  = repmat(h,size(mrdm.set));
            %         mRDM        = addstruct(mRDM,mrdm);
            %         %
            
            D.maxidx        = repmat(m,length(indx),1);
            D.omega         = omega;
            D.omegaN        = omegaN;
            D.logtheta      = logtheta;
            D.logEvidence   = logEvidence;
            %D.logEvidenceSplit = logEvidenceSplit;
            
            TT = addstruct(TT,D);
            
            %sh1_rsa('plot_omega_prior',D);
            %sh1_rsa('plot_prior',D);
            %sh1_rsa('plot_tau',D);
        end
        figure;
        subplot(3,2,[1 2])
        lineplot(TT.maxidx,TT.logEvidence,'style_shade');ylabel('logEvidence')
        
        subplot(3,2,[3 4]);
        barplot(TT.maxidx,TT.logtheta(:,5:end),...
            'subset',ismember(TT.maxidx,[1:10:100]));ylabel('logTau')
        
        subplot(3,2,[5 6]);
        barplot(TT.maxidx,TT.omegaN,...
            'subset',ismember(TT.maxidx,[1:10:100]));ylabel('Omega')
        
        
        varargout={TT,mRDM};
    case 'fit_model_EB_nonlin_singlemodel' % fit single model to compare log-evidence
        glm = 3;
        useolddata = 0;
        modelTerms = [1 2 4 6];
        nonlinP0 = [0 0 0 0];
        vararginoptions(varargin(:),{'glm','useolddata','modelTerms','nonlinP0'})
        
        % load data set
        if useolddata
            %T       = load(fullfile(regDir,'distances_glm1.mat'));
            T       = load(fullfile(regDir,'allRDMs_ROI_all_glm3.mat'));
            %T       = T.RDM1;
            T.RDM   = bsxfun(@rdivide,T.RDM,T.numVox);
        else
            T       = load(fullfile(regDir,sprintf('distances_sepPerm.glm%d.mat',glm)));
        end
        
        %regions=unique(T.region);
        regions=[2 3 4 5 7]';
        for r=regions'
            fprintf('Estimating ROI distance %s...\n',regname{r});
            for h=[1 2]
                for m=1:numel(modelTerms)
                    indx = find(T.region==r & T.hemis==h);
                    D=getrow(T,indx);
                    for s=1:length(D.subj)
                        Sigma(:,:,s)=reshape(D.Sigma(s,:),8,8);
                    end;
                    
                    % organize model RDMs
                    Model.fcn = @sh1_getRDMmodelTau1;
                    Model.constantParams = {chunk_set(D.subj),...       % Model terms
                        modelTerms(m),...               % Chunk set
                        'sqEuclidean'...            % Distance term
                        };
                    Model.numComp   = numel(Model.constantParams{2}); % Number of linear components
                    Model.numPrior  = numel(Model.constantParams{2}); % Number of prior variances on parameters
                    Model.numNonlin = numel(Model.constantParams{2}); % Number of nonlinear parameters
                    Model.nonlinP0  = nonlinP0(m);%zeros(size(Model.constantParams{2})); % Starting value of nonlinear(mixing) parameters
                    
                    % estimate
                    [w,logEvidence(indx,m),lt] = rsa_fitModelHierarchEB(Model,D.RDM,Sigma,9,D.effVox);
                    
                    logprior(indx,m)  = repmat(lt(1),length(indx),1);
                    logtau(indx,m)    = repmat(lt(2),length(indx),1);
                    term(indx,m)      = repmat(modelTerms(m),length(indx),1);
                    
                    % normalize omega with mean of each RDM
                    M = feval(Model.fcn,lt(Model.numPrior+1:Model.numPrior+Model.numNonlin),...
                        Model.constantParams{:});
                    omega(indx,m) = w.*permute(sqrt(mean(M.RDM.^2,2)),[3 1 2]);
                    
                end % model term
            end;
        end;
        T.omega = omega;
        T.logtheta = [logprior, logtau];
        T.logEvidence = logEvidence;
        T.term = term;
        %T.logEvidenceSplit = logEvidenceSplit;
        
        % plot result
        sh1_rsa('plot_logEvidence',T);
        sh1_rsa('plot_omega',T);
        sh1_rsa('plot_prior',T);
        sh1_rsa('plot_tau',T);
        varargout={T};
    case 'fit_model_EB_nonlin_multiguess'
        glm         = 1;
        useolddata  = 1;
        modelTerms  = [1 2 4 6];
        Niter       = 100;
        vararginoptions(varargin(:),{'glm','useolddata','modelTerms'})
        
        % load data set
        if useolddata
            %T       = load(fullfile(regDir,'distances_glm1.mat'));
            T       = load(fullfile(regDir,'allRDMs_ROI_all_glm3.mat'));
            %T       = T.RDM1;
            T.RDM   = bsxfun(@rdivide,T.RDM,T.numVox);
        else
            T       = load(fullfile(regDir,sprintf('distances_sepPerm.glm%d.mat',glm)));
        end
        
        %regions=unique(T.region);
        regions=[1:2]';
        mRDM = [];
        for r=regions'
            for h=[1 2]
                indx = find(T.region==r & T.hemis==h);
                D=getrow(T,indx);
                for s=1:length(D.subj)
                    Sigma(:,:,s)=reshape(D.Sigma(s,:),8,8);
                end;
                
                for iter = 1:Niter
                    iter
                    initialguess = rand(size(modelTerms))*100-50;
                    initialguess(end) = 0;
                    
                    % organize model RDMs
                    Model.fcn = @sh1_getRDMmodelTau1;
                    Model.constantParams = {chunk_set(D.subj),...       % Model terms
                        modelTerms,...               % Chunk set
                        'sqEuclidean'...            % Distance term
                        };
                    Model.numComp   = numel(Model.constantParams{2}); % Number of linear components
                    Model.numPrior  = numel(Model.constantParams{2}); % Number of prior variances on parameters
                    Model.numNonlin = numel(Model.constantParams{2}); % Number of nonlinear parameters
                    Model.nonlinP0  = initialguess; % Starting value of nonlinear(mixing) parameters
                    
                    % estimate
                    [w(:,:,iter),levid(:,iter),lt(iter,:)] = rsa_fitModelHierarchEB(Model,D.RDM,Sigma,9,D.effVox);
                    
                    mevid = mean(levid,1);
                end
                maxidx = find(mevid==max(mevid),1,'first');
                % choose one
                logtheta(indx,:)    = repmat(lt(maxidx,:),length(indx),1);
                w                   = squeeze(w(:,:,maxidx));
                logEvidence(indx,:) = levid(:,maxidx);
                
                % normalize omega with mean of each RDM
                M = feval(Model.fcn,lt(Model.numPrior+1:Model.numPrior+Model.numNonlin),...
                    Model.constantParams{:});
                normX = sqrt(mean(M.RDM.^2,2));
                omega(indx,:) = w;
                omegaN(indx,:) = bsxfun(@times,w,normX');
                
                % save resultant RDMs with estimated tau
                A = feval(Model.fcn,lt(Model.numPrior+1:Model.numPrior+Model.numNonlin),...
                    1,Model.constantParams{2:end}); % chunk set A
                B = feval(Model.fcn,lt(Model.numPrior+1:Model.numPrior+Model.numNonlin),...
                    2,Model.constantParams{2:end}); % chunk set B
                mrdm.RDM    = [A.RDM;B.RDM];
                mrdm.set    = kron([1 2]',ones(Model.numComp,1));
                mrdm.comp   = repmat(Model.constantParams{2}',2,1);
                mrdm.name   = {A.name{:};B.name{:}};
                mrdm.region = repmat(r,size(mrdm.set));
                mrdm.hemis  = repmat(h,size(mrdm.set));
                mRDM        = addstruct(mRDM,mrdm);
                %
            end;
        end;
        T.omega=omega;
        T.logtheta=logtheta;
        T.logEvidence=logEvidence;
        %T.logEvidenceSplit = logEvidenceSplit;
        
        sh1_rsa('plot_omega_prior',T);
        sh1_rsa('plot_prior',T);
        sh1_rsa('plot_tau',T);
        varargout={T,mRDM};
        
    case 'plot_omega'
        T=varargin{1};
        
        figure('name','omega');
        for reg=1:8;
            subplot(2,4,reg);
            xpos=myboxplot(T.hemis,T.omega,'subset',T.region==reg,...
                'CAT',CAT);
            drawline(0,'dir','horz');
            title(regname{reg},'fontsize',12);
            %set(gca,'XTick',[xpos(1:4)],'XTicklabel',omega_names,'fontsize',12);
            %rotateXLabels(gca,45);
        end
    case 'plot_prior'
        T=varargin{1};
        
        figure('name','prior');
        for reg=1:8;
            subplot(2,4,reg);
            xpos=barplot(T.hemis,T.logtheta(:,1:4),'subset',T.region==reg,...
                'CAT',CAT);
            drawline(0,'dir','horz');
            title(regname{reg},'fontsize',12);
            %set(gca,'XTick',[xpos(1:4)],'XTicklabel',omega_names,'fontsize',12);
            %rotateXLabels(gca,45);
        end
    case 'plot_tau'
        T=varargin{1};
        
        figure('name','tau');
        for reg=1:8;
            subplot(2,4,reg);
            xpos=barplot(T.hemis,T.logtheta(:,5:end),'subset',T.region==reg,...
                'CAT',CAT);
            drawline(0,'dir','horz');
            title(regname{reg},'fontsize',12);
            %set(gca,'XTick',[xpos(1:4)],'XTicklabel',omega_names,'fontsize',12);
            %rotateXLabels(gca,45);
        end
    case 'plot_logEvidence'
        T=varargin{1};
        
        figure('name','logEvidence');
        for reg=1:8;
            subplot(2,4,reg);
            xpos=myboxplot(T.hemis,T.logEvidence,'subset',T.region==reg,...
                'CAT',CAT);
            drawline(0,'dir','horz');
            title(regname{reg},'fontsize',12);
            %set(gca,'XTick',[xpos(1:4)],'XTicklabel',omega_names,'fontsize',12);
            %rotateXLabels(gca,45);
        end
    case 'plot_omega_prior'
        T=varargin{1};
        
        figure('name','omega');
        for reg=1:8;
            subplot(2,8,reg);
            xpos=myboxplot(T.hemis,T.omega,'subset',T.region==reg,...
                'CAT',CAT);
            drawline(0,'dir','horz');
            title(regname{reg},'fontsize',12);
            %set(gca,'XTick',[xpos(1:4)],'XTicklabel',omega_names,'fontsize',12);
            %rotateXLabels(gca,45);
            
            % plot log prior
            subplot(2,8,reg+8);
            data  = getrow(T,T.region==reg);
            scale = max(max(abs(data.logtheta(:,1:4))));
            Ldata = getrow(T,T.region==reg&T.hemis==1);
            Rdata = getrow(T,T.region==reg&T.hemis==2);
            imagesc_rectangle([mean(Ldata.logtheta(:,1:4),1),mean(Rdata.logtheta(:,1:4),1)],...
                'scale',[-scale, scale]);
        end
    case 'plot_reliability_btw' % plot between subject reliabitliy of RDM
        C = varargin{1};
        region = [1:8];
        barwidth = 0.8;
        
        Nregion = numel(region);
        figure('name',sprintf('Between subject reliability of RDM'),'color','w');
        for reg = region;
            fprintf('%s',regname{reg})
            count=0;
            for h=[2 1]
                % right hemi
                subplot(2,Nregion,reg+count*Nregion);
                subset = C.region==reg&C.hemis==h;
                upperbound  = nanmean(C.within_u(subset,:));
                se          = nanmean(C.within_u(subset,:))/sqrt(nancount(C.within_u(subset,:)));
                
                xpos = barplot([],[C.within_l,C.across,C.within_l_reduced],'subset',subset,'barwidth',barwidth);
                [tval, p]=ttest_mc(C.within_l,C.across,1,'paired','subset',subset)
                
                drawline(0,'dir','horz','lim',[xpos(1)-barwidth,xpos(2)+barwidth]); hold on
                drawline(upperbound,'dir','horz','color',[0.5 0.5 0.5],'lim',...
                    [xpos(1)-barwidth,xpos(2)+barwidth],...
                    'linewidth',2,'error',se,'errorcolor',[0.8 0.8 0.8])
                xpos = barplot([],[C.within_l,C.across,C.within_l_reduced],'subset',subset,'barwidth',barwidth); hold on
                
                title([regname{reg},'-',hemName{h}],'fontsize',12);
                ylabel('Pearson r','fontsize',12);xlabel('');
                set(gca,'ylim',[-0.2 0.8],'xticklabel',{'W','B'},'fontsize',12)
                if p<0.05
                    drawline(0.78,'dir','horz','lim',xpos([1,2]),'color',[1 0 0]);
                elseif p<0.01
                    drawline(0.78,'dir','horz','lim',xpos([1,2]),'color',[0 1 0]);
                elseif p<0.001
                    drawline(0.78,'dir','horz','lim',xpos([1,2]),'color',[0 0 1]);
                end
                count=count+1;
            end
        end
        
        
        
        
    case 'showRDMs'
        T=varargin{1};
        glm=varargin{2};
        for reg=1:8
            for h=1:2
                figName = sprintf('RDMs.%s.%s.glm%d)',regname{reg},hemName{h},glm);
                figure('name',figName);
                tmp = getrow(T,T.region==reg&T.hemis==h);
                rsa.fig.imageRDMs(tmp);
                
                %saveas(h,fullfile(figDir,figName),'fig');
            end
        end
    otherwise
        disp('there is no such case.')
end;

function RDMout = getCommon(RDMin, seqnum);
% RDMin : vectorised RDM
% RDMout: vectorised RDM
% seqnum: seq of no interest

I = true(nSequence);
I(diag(diag(I))) = false;
I(seqnum,:) = false;
I(:,seqnum) = false;

RDMout = RDMin(:,squareform(I));


function out = sh1_calcDistancePerm(data,SPM,seqType);
[RDM,Sw,effVox]=rsa_distanceLDCsepPerm(data,SPM,seqType);
out=[RDM(:);Sw(:);effVox];
