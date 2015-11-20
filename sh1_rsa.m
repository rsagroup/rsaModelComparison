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

%baseDir         = '/Users/joern/Projects/sh1';   % For Joern
baseDir         = '/Volumes/DATA/MotorControl/data/SequenceLearning/sh1';   % For AY
%baseDir         = '/Users/Atsushi/Dropbox/work/UK/Research/MotorControl/data/SequenceLearning/sh1';   % For AY 2

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
Sequence         = [1 2 3 4;... % 2 3 3 3
    4 3 2 1;...
    5 6 7 8;...
    8 7 6 5;...
    1 2 7 8;...
    8 7 2 1;... % 3 3 3 2
    5 6 3 4;...
    4 3 6 5]; % 3 3 3 2

Sequence_B =   [4 3 2 1; % 3 3 3 2
    5 6 7 8;
    8 7 6 5;
    1 2 3 4;
    4 3 6 5;
    1 2 7 8; % 2 3 3 3
    8 7 2 1;
    5 6 3 4];% 2 3 3 3
Contrasts      = {'seq1','seq2','seq3','seq4','Seq5','seq6','seq7','seq8',...
    'allMovement', 'F-contrast', '2333','3332'};

varargout = {};

%% Graph parameters
CAT.facecolor   = {[0 0 0.7],[0.3 0.3 1],[0.7 0 0],[1 0.6 0.6]}; %[1 0.1 0.1],
CAT.fillcolor   = CAT.facecolor;
CAT.linecolor   = CAT.facecolor;
CAT.edgecolor   = CAT.facecolor;
CAT.capwidth    = 0.01;
CAT.markertype  = {'o'};
CAT.markersize  = {8};
CAT.markercolor = CAT.fillcolor;
CAT.markerfill  = {[1 1 1]};
CAT.gapwidth    = {[0.2 0.2 0.2 0.2]};
CAT.linewidth   = {1};
CAT.errorwidth  = {1};
CAT.errorcolor  = {[0 0 0]};

omega_names = {'One-digit','Two-digit','Chunk','Sequence'};

%% Main operation
switch(what)
    case 'ROI_distraw'                                                      % Extracts ROI data and calculates distances
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
    case 'ROI_distraw_effVox_subspace'                 % Extracts ROI data and calculates distances
        T=[];
        sn = varargin{1};
        numVox_subspace = 160;
        Niter = 100;
        glm = 3;
        ROI = 'all';
        fname 	= 'movement';       % using metric file which compares any movement against rest
        regions = 1:12;
        hemi = 1:2;
        vararginoptions({varargin{2:end}},{'chunkset','glm','ROI','fname','regions','hemi'});
        
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
                        
                        % apply random subspace method
                        if numVox<numVox_subspace
                            index = randperm(numVox,numVox);
                        else                            
                            index = randperm(numVox,numVox_subspace);
                        end
                        for i = Niter
                            % get the distances
                            [S.RDM,Sw,S.effVox,S.trSS]=rsa_distanceLDCsepPerm(data(:,index),SPM,E.seqType);
                            S.Sigma  = Sw(:)';
                            S.region = reg;
                            S.subj  = s;
                            S.hemis = hem;
                            S.numVox = numVox_subspace;
                            S.iter = i;
                            
                            T=addstruct(T,S);
                            
                            if numVox<numVox_subspace
                                S.numVox = numVox;
                                break;
                            end
                        end
                    end;
                end;
            end;
        end;
        varargout={T};
        save(fullfile(regDir,'distances_sepPerm_subspace.mat'),'-struct','T');
    case 'searchlight_distraw'                                              % Calc distance, sigma, effective no. of voxesl from time series data
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
    case 'searchlight_fit_EB_lin'
        sn = varargin{1};
        method = 'ZellnerIndivid';
        hemi = [1 2];
        %logtau = [-6.18 -0.243 -0.898 1];
        %logtau = [-10 2.0217 2.0217 0]; % estimated from behaviour (error rate)
        reguse = [1 2 4 6];
        glm =3;
        vararginoptions(varargin(2:end),{'method','hemi','logtau','reguse'});                
        
        Model = sh1_getRDMmodelTau1(logtau,chunk_set(sn),reguse,'normSqEuclidean');
        Model.numNonlin = 0;
        Model.prior = method;
        for h=hemi
            for s=1:length(sn)
                distFiles{s}=fullfile(baseDir,glmName{glm},subj_name{sn(s)},sprintf('%s_dist_raw.nii',subj_name{sn(s)}));
                sigmFiles{s}=fullfile(baseDir,glmName{glm},subj_name{sn(s)},sprintf('%s_sigma_raw.nii',subj_name{sn(s)}));
                nvoxFiles{s}=fullfile(baseDir,glmName{glm},subj_name{sn(s)},sprintf('%seffVox.nii',subj_name{sn(s)}));
                whiteFiles{s}=fullfile(caretDir,['x' subj_name{sn(s)}],hemName{h},[hem{h} '.WHITE.coord']);
                pialFiles{s}=fullfile(caretDir,['x' subj_name{sn(s)}], hemName{h},[hem{h} '.PIAL.coord']);
            end;
            for i=1:4
                outFiles{i}=fullfile(caretDir,'fsaverage_sym', hemName{h},sprintf('%s.%sOmega%d.metric',hem{h},method,reguse(i)));
            end;
            outFiles{5}= fullfile(caretDir,'fsaverage_sym', hemName{h},sprintf('%s.%sNumSubj.metric',hem{h},method));
            rsa_fitModelSurfEB(Model,distFiles,sigmFiles,9,nvoxFiles,whiteFiles,pialFiles,outFiles);
        end;
            
    case 'ROI_reliability_btw'                                              % Calculate between subject reliability of RDMs
        
        useallRDM = 0; % remove non-shared sequences or not
        glm = [3];
        vararginoptions(varargin,{'useallRDM','glm'});
        
        C = [];
        for g = glm
            % load data set
            T = load(fullfile(regDir,sprintf('distances_sepPerm.glm%d.mat',g)));
            
            chunksets = {'A','B'};
            for reg = [1:8]
                for h = [1 2]
                    for chunkset = [1 2]
                        % calc within chunkset value
                        D = getrow(T,T.region==reg&T.hemis==h&T.chunkset==chunksets{chunkset});
                        D_ = getrow(T,T.region==reg&T.hemis==h&T.chunkset~=chunksets{chunkset});
                        
                        rdm_all = nanmean(D.RDM,1); % all rdm for the same chunkset
                        rdm_all_ = nanmean(D_.RDM,1); % all rdm for other chunkset
                        
                        if useallRDM==0
                            % we must choose only the distance between physically identical sequences here
                            rdm_all     = getCommon(rdm_all,[3 5 7]);
                            rdm_all_    = getCommon(rdm_all_,[3 5 7]);
                        end
                                                
                        subj = unique(D.subj);
                        for s=1:numel(subj)
                            idx = D.subj==subj(s);
                            
                            rdm = D.RDM(idx,:);
                            rdm_loo = nanmean(D.RDM(~idx,:),1); % leave-one-subject-out rdm
                                                        
                            if useallRDM==0
                                % we must choose only the distance between physically identical sequences here                                
                                rdm      = getCommon(rdm,[3 5 7]);
                                rdm_loo  = getCommon(rdm_loo,[3 5 7]);                                
                            end
                            
                            % visualise RDMs
                            Tmp.RDM = [rdm_all;rdm_loo;rdm_all_;rdm];
                            Tmp.name = {sprintf('avrgRDM (%s)',chunksets{chunkset});...
                                        sprintf('avrgRDM w.o. subj:%s (%s)',subj_name{subj(s)},chunksets{chunkset});...
                                        sprintf('avrgRDM (%s)',chunksets{3-chunkset});...
                                        sprintf('RDM of subj:%s (%s)',subj_name{subj(s)},chunksets{chunkset})};
                            %rsa.fig.imageRDMs(Tmp);pause(0.005);title(sprintf('ROI=%s',regname{reg}))
                                    
                            R.within_u          = corr(rdm',rdm_all','type','Pearson'); % upper bound
                            R.within_l          = corr(rdm',rdm_loo','type','Pearson'); % lower bound                            
                            R.across            = corr(rdm',rdm_all_','type','Pearson'); % across chunk set corr                                                                                    
                                                        
                            R.subj      = subj(s);
                            R.chunkset  = chunkset;
                            R.region    = reg;
                            R.hemis     = h;
                            R.glm       = g;
                            
                            C = addstruct(C,R);
                        end
                    end
                end
            end
        end
        
        % get Fisher-z-transformed version of C
        Z = C;
        Z.within_u = fisherz(Z.within_u);
        Z.within_l = fisherz(Z.within_l);        
        Z.across   = fisherz(Z.across);
        
        for g=glm
            %figure('name',sprintf('reliability (glm=%d)',glm));
%             D = getrow(Z,Z.glm==g&Z.chunkset==1);
%             sh1_rsa('plot_reliability_btw',D);
%             
%             D = getrow(Z,Z.glm==g&Z.chunkset==2);
%             sh1_rsa('plot_reliability_btw',D);
            
            D = getrow(Z,Z.glm==g);
            sh1_rsa('plot_reliability_btw',D);
        end
        
        varargout = {C};
    case 'ROI_btwchunkset_corr'                                             % Calculate between subject reliability of RDMs
        
        useallRDM = 0; % remove non-shared sequences or not
        glm = [3];
        removeSeq = [3 5 7];
        corrfun = 'Pearson';
        ztransform = 1;
        vararginoptions(varargin,{'useallRDM','glm','corrfun','removeSeq','ztransform'});
        
        % get model rdm for chunk
        MA = sh1_getRDMmodelTau1(10,1,4,'sqEuclidean'); MA.name{1} = 'chunk-setA';
        MB = sh1_getRDMmodelTau1(10,2,4,'sqEuclidean'); MB.name{1} = 'chunk-setB';
        if useallRDM==0
            MA.RDM = getCommon(MA.RDM,removeSeq);
            MB.RDM = getCommon(MB.RDM,removeSeq);
        end
        switch corrfun
            case 'corrN'
                Corr_model = real(corrN(MA.RDM',MB.RDM'));
            otherwise
                Corr_model = corr(MA.RDM',MB.RDM','type',corrfun);
        end
        figure(1); rsa.fig.imageRDMs(addstruct(MA,MB));        
        fprintf('Correlation between chunk RDM for set A and B (%s) : %1.2f\n',corrfun,Corr_model);
            
        G = []; All=[];
        for g = glm
            % load data set
            T = load(fullfile(regDir,sprintf('distances_sepPerm.glm%d.mat',g)));
                        
            for reg = [1:8]
                for h = [1 2]
                    D = getrow(T,T.region==reg&T.hemis==h);
                    
                    if useallRDM==0
                        D.RDM = getCommon(D.RDM,removeSeq);
                    end
                    
                    % get logical mask
                    idxA        = double(D.chunkset=='A');
                    idxB        = double(D.chunkset=='B');
                    idxWithin   = logical(idxA*idxA'+idxB*idxB');
                    idxBetween  = ~idxWithin;
                    idxAll      = true(size(idxWithin));
                    diagonal    = logical(eye(size(idxWithin)));
                    idxWithin(diagonal) = false;
                    idxBetween(diagonal) = false;
                    idxAll(diagonal) = false;
                    
                    % get corelation matrix
                    switch corrfun
                        case 'corrN'
                            Ctmp = real(corrN(D.RDM',D.RDM'));
                        otherwise
                            Ctmp = corr(D.RDM','type',corrfun);
                    end
                    Ctmp(diagonal) = NaN;
                    if ztransform==1
                        Ctmp = fisherz(Ctmp);
                    end

                    % apply mask (df is the number of subjects)                    
                    Wtmp = Ctmp;
                    Wtmp(idxBetween) = NaN;
                    Btmp = Ctmp;
                    Btmp(idxWithin) = NaN;
                    
                    
                    
%                     % vectorise and give indicator
%                     R.r = Ctmp(tril(idxAll));
%                     R.z = fisherz(R.r);
%                     R.type(idxWithin(tril(idxAll)),1) = 1;
%                     R.type(idxBetween(tril(idxAll)),1) = 2;

                    R.corr = [nanmean(Wtmp,2);nanmean(Btmp,2)];
                    R.type = [ones(size(Wtmp,1),1);2*ones(size(Wtmp,1),1)];
                    R.subj = [D.subj;D.subj];
                    R.region = repmat(reg,size(R.corr));
                    R.hemis = repmat(h,size(R.corr));
                    R.glm = repmat(g,size(R.corr));
                    
                    All = addstruct(All,R);

                end
            end
        end
        % ttest
        S = [];
        for r=unique(All.region)'
            for h=unique(All.hemis)'
                fprintf('paired t-test (%s-%s)\n',regname{r},hemName{h})
                subsetW = All.region==r&All.hemis==h&All.type==1;
                subsetB = All.region==r&All.hemis==h&All.type==2;
                [t,p] = ttest_mc(All.corr(subsetW),All.corr(subsetB),1,'paired');
                fprintf('t(%d) = %2.3f (p = %2.4f)\n---\n',numel(unique(All.subj))-1,t,p);
                
                sta.t=t;
                sta.p=p;
                sta.region = r;
                sta.hemis = h;
                
                S = addstruct(S,sta);
            end
        end
        
        
        figure;
        subplot(2,1,1);barplot([All.hemis All.region],All.corr,'split',All.type,...
            'leg',{'within','between'},'leglocation','northeast','gapwidth',[2 0.5 0.2]);
        ylabel(sprintf('z-transformed correlation (%s)',corrfun),'fontsize',13);
        xlabel('ROI','fontsize',13);
        set(gca,'ylim',[-0.1 0.5])
        
        subplot(2,1,2);barplot([S.hemis,S.region],log10(S.p),'gapwidth',[2 0.5 0.2]);
        ylabel('log(p) (one-sided paired-ttest)','fontsize',13);
        xlabel('ROI','fontsize',13);
        drawline([log10(0.001) log10(0.01) log10(0.05) log10(0.1)],'dir','horiz','color',[1 0 0])
        set(gca,'ylim',[-2 0],'ytick',[log10(0.001) log10(0.01) log10(0.05) log10(0.1)])
        
        
        varargout = {All};
    case 'ROI_fit_model_noise_ceiling'                                      % Calculate noise ceiling with linear model
        method = 'GLS';
        glm = 3;
        vararginoptions(varargin,{'method','ceiling','glm'}); 
        
        % load data set
        T = load(fullfile(regDir,sprintf('distances_sepPerm.glm%d.mat',glm)));
        
        T.chunk_set = chunk_set(T.subj)'; 
        regions=unique(T.region);
        T.ceiling = NaN(size(T.subj)); % col-1: upper bound, col-2: lower bound        
        
        for r=regions'
            for h=[1 2]
                % Get the data to be fitted 
                indx = find(T.region==r & T.hemis==h);

                % Generate the model for all subject you are fitting 
                for c = [1 2]
                    % for upper bound
                    dist    = getrow(T,T.region==r&T.hemis==h&T.chunk_set==c);
                    Xu      = nanmean(dist.RDM,1);
                    Xu      = Xu / sqrt(sum(Xu.*Xu));
                    Upper(c).RDM = Xu;
                    % for lower bound
                    for s = 1:length(dist.subj)
                        idx    = true(size(dist.subj));
                        idx(s) = false;
                        Xl     = nanmean(dist.RDM(idx,:),1);
                        Xl     = Xl / sqrt(sum(Xl.*Xl));
                        Lower(c,s).RDM = Xl;
                    end                    
                end
                
                % Assemble the design matrices 
                for s=1:length(indx)
                    Sigma(:,:,s)=reshape(T.Sigma(indx(s),:),8,8);
                end;

                % Do the fit
                switch(method)
                    case 'OLS'
                        for i=1:2 
                            dist = getrow(T,T.region==r&T.hemis==h&T.chunk_set==i);                            
                            for s = 1:length(dist.subj)                                
                                idx = find(T.chunk_set==i&T.region==r&T.hemis==h&T.subj==dist.subj(s));
                                T.ceiling(idx,1) = rsa_fitModelOLS(Upper(i),T.RDM(idx,:)); 
                                T.ceiling(idx,2) = rsa_fitModelOLS(Lower(i,s),T.RDM(idx,:)); 
                            end
                        end; 
                    case 'GLS' 
                        for i=1:2 
                            dist = getrow(T,T.region==r&T.hemis==h&T.chunk_set==i);
                            for s = 1:length(dist.subj)                                
                                idx = find(T.chunk_set==i&T.region==r&T.hemis==h&T.subj==dist.subj(s));
                                T.ceiling(idx,1) = rsa_fitModelGLS(Upper(i),T.RDM(idx,:),T.Sigma(idx,:),9,T.effVox(idx,1));
                                T.ceiling(idx,2) = rsa_fitModelGLS(Lower(i,s),T.RDM(idx,:),T.Sigma(idx,:),9,T.effVox(idx,1));
                            end                                                        
                        end; 
                    case 'IRLS' 
                        for i=1:2 
                            dist = getrow(T,T.region==r&T.hemis==h&T.chunk_set==i);
                            for s = 1:length(dist.subj)                                
                                idx = find(T.chunk_set==i&T.region==r&T.hemis==h&T.subj==dist.subj(s));
                                T.ceiling(idx,1) = rsa_fitModelIRLS(Upper(i),T.RDM(idx,:),T.Sigma(idx,:),9,T.effVox(idx,1));
                                T.ceiling(idx,2) = rsa_fitModelIRLS(Lower(i,s),T.RDM(idx,:),T.Sigma(idx,:),9,T.effVox(idx,1));
                            end                            
                        end;
                end; 
            end;
        end;
        subplot(2,1,1); 
        barplot(T.region,[T.ceiling],'subset',T.hemis==1); 
        subplot(2,1,2); 
        barplot(T.region,[T.ceiling],'subset',T.hemis==2); 
        set(gcf,'Name',method); 
        T.method = repmat({method},size(T.subj));
        varargout={T};    
        
    case 'fit_model_EB_lin'
        method = 'ZellnerIndivid'; 
        %logtau = [-10 2.0217 2.0217 0]; % estimated from behaviour (error rate)
        logtau = [-6.265 -0.243 -0.898 1]; % estimated by EB nonlin
        ceiling = 1;
        glm = 3;
        sn = [2:6,8,10:15];
        regions = [1:8];
        vararginoptions(varargin,{'method','logtau','ceiling','glm','sn','regions'}); 
        
        % load data set
        T = load(fullfile(regDir,sprintf('distances_sepPerm.glm%d.mat',glm)));
        T.chunk_set = chunk_set(T.subj)'; 
        %regions=unique(T.region);
        
        T = getrow(T,ismember(T.subj,sn));
        
        % Seperate models for the two chunk sets 
        Mod(1) = sh1_getRDMmodelTau1(logtau,1,[1 2 4 6],'normSqEuclidean'); 
        Mod(2) = sh1_getRDMmodelTau1(logtau,2,[1 2 4 6],'normSqEuclidean'); 
        Mod(1).numNonlin = 0;   
        Mod(2).numNonlin = 0;   

        for r=regions
            for h=[1 2]
                % Get the data to be fitted 
                indx = find(T.region==r & T.hemis==h);

                % Generate the model for all subject you are fitting 
                Model = sh1_getRDMmodelTau1(logtau,chunk_set(T.subj(indx)),[1 2 4 6],'normSqEuclidean'); 
                Model.numNonlin = 0;   

                % Assemble the design matrices 
                for s=1:length(indx)
                    Sigma(:,:,s)=reshape(T.Sigma(indx(s),:),8,8);
                end;

                % Do the fit
                switch(method)
                    case 'OLS'
                        for i=1:2 
                            indx = find(T.chunk_set==i & T.region==r & T.hemis==h);
                            T.omega(indx,:)=rsa_fitModelOLS(Mod(i),T.RDM(indx,:)); 
                        end; 
                    case 'GLS' 
                        for i=1:2 
                            indx = find(T.chunk_set==i & T.region==r & T.hemis==h);
                            T.omega(indx,:)=rsa_fitModelGLS(Mod(i),T.RDM(indx,:),T.Sigma(indx,:),9,T.effVox(indx,1)); 
                        end; 
                    case 'IRLS' 
                        for i=1:2 
                            indx = find(T.chunk_set==i & T.region==r & T.hemis==h);
                            T.omega(indx,:)=rsa_fitModelIRLS(Mod(i),T.RDM(indx,:),T.Sigma(indx,:),9,T.effVox(indx,1)); 
                        end;
                        
                   case {'Ridge','Zellner','RidgeIndivid','ZellnerIndivid'}
                        Model.prior=method; 
                        [T.omega(indx,:),T.logEvidence(indx,:),lt]=rsa_fitModelHierarchEB(Model,T.RDM(indx,:),Sigma,9,T.effVox(indx,1));
                        T.theta(indx,:)=repmat(lt,length(indx),1);
                end; 
            end;
        end;
        
        C = [];
        if ceiling==1
            C = sh1_rsa('ROI_fit_model_noise_ceiling','glm',glm,'method','GLS');
        end
        T.method = repmat({method},size(T.subj));
                
        sh1_rsa('plot_omega_new',T,'ceiling',C);
        
        % summary logEvidence
        pivottable(T.region,T.hemis,T.logEvidence,'nanmean');
        fprintf('total logEvidence: %2.3f\n',sum(T.logEvidence));
        
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
    case 'fit_model_EB_nonlin_singlemodel'                                  % fit single model to compare log-evidence from multiple starting points
        glm = 3;
        useolddata = 0;
        modelTerms = [1 2 4 6];
        nonlinP0 = [0 0 0 0];
        Niter = 4;
        regions = [1 2 3 4 5 7 8];
        hemi = [1 ];
        sn = [2:15];
        
        vararginoptions(varargin(:),{'glm','useolddata','modelTerms','nonlinP0','Niter','regions','hemi','sn'})
        
        % load data set
        if useolddata
            %T       = load(fullfile(regDir,'distances_glm1.mat'));
            T       = load(fullfile(regDir,'allRDMs_ROI_all_glm3.mat'));
            %T       = T.RDM1;
            T.RDM   = bsxfun(@rdivide,T.RDM,T.numVox);
        else
            T       = load(fullfile(regDir,sprintf('distances_sepPerm.glm%d.mat',glm)));
        end
        
        T = getrow(T,ismember(T.subj,sn));
        
        B=[];
        
        for r=regions
            fprintf('Estimating ROI distance %s...\n',regname{r});
            for h=hemi
                % Get the data and sigma for each person
                indx = find(T.region==r & T.hemis==h);
                D = getrow(T,indx);
                for s=1:length(D.subj)
                    Sigma(:,:,s)=reshape(D.Sigma(s,:),8,8);
                end;
                
                % Fit the Null model
                logEvidence0 = rsa_fitModelNull(D.RDM,D.Sigma,9,D.effVox);
                
                for m = 1:numel(modelTerms)
                    m
                    for i=1:Niter % try multiple starting points
                        % organize model RDMs
                        Model.fcn = @sh1_getRDMmodelTau1;
                        Model.constantParams = {chunk_set(D.subj),...       % Model terms
                            modelTerms(m),...               % Chunk set
                            'sqEuclidean'...            % Distance term
                            };
                        Model.numComp   = numel(Model.constantParams{2}); % Number of linear components
                        Model.numPrior  = numel(Model.constantParams{2}); % Number of prior variances on parameters
                        Model.numNonlin = numel(Model.constantParams{2}); % Number of nonlinear parameters
                        %Model.nonlinP0  = nonlinP0(m);%zeros(size(Model.constantParams{2})); % Starting value of nonlinear(mixing) parameters
                        Model.nonlinP0  = 3*rand(1)-1.5;
                        
                        % estimate
                        [w,levid,lt] = ...
                            rsa_fitModelHierarchEB(Model,D.RDM,Sigma,9,D.effVox,'minimizeLength',1000);
                        
                        A.logprior(:,i)  = repmat(lt(1),length(indx),1);
                        A.logtau(:,i)    = repmat(lt(2),length(indx),1);
                        A.logEvidence(:,i) = levid-logEvidence0; % increase of logEvidence from null model
                        A.term(:,1)      = repmat(modelTerms(m),length(indx),1);
                        
                        % normalize omega with norm of each RDM
                        M = feval(Model.fcn,lt(Model.numPrior+1:Model.numPrior+Model.numNonlin),...
                            Model.constantParams{:});
                        A.omega(:,i) = w.*permute(sqrt(mean(M.RDM.^2,2)),[3 1 2]);
                        A.omegaN(:,i) = w;
                        
                        A.iter = repmat(i,length(indx),1);
                        A.region = repmat(r,length(indx),1);
                        A.hemis = repmat(h,length(indx),1);
                        A.subj = D.subj;
                    end
                    % Before adding, pick the best one
                    [~,i]=max(sum(A.logEvidence));
                    A.logprior = A.logprior(:,i);
                    A.logtau  = A.logtau(:,i);
                    A.logEvidence = A.logEvidence(:,i); % increase of logEvidence from null model
                    A.term(:,1)      = repmat(modelTerms(m),length(indx),1);
                    
                    % normalize omega with norm of each RDM
                    M = feval(Model.fcn,lt(Model.numPrior+1:Model.numPrior+Model.numNonlin),...
                        Model.constantParams{:});
                    A.omega = A.omega(:,i);
                    A.omegaN = A.omegaN(:,i);
                    
                    B = addstruct(B,A);
                end % model term
            end;
        end;
        pivottable([B.hemis B.region],B.term,B.logEvidence,'mean');
        pivottable([B.hemis B.region],B.term,B.logtau,'log(mean(exp(x)))');
        L  = pivottable([B.hemis B.region],B.term,B.logEvidence,'mean');
        [~,m]=max(L'); 
        lT = pivottable([B.hemis B.region],B.term,B.logtau,'log(mean(exp(x)))');
        for i=1:3 
            fprintf('logtau (%d):%2.3f\n',i,log(mean(exp(lT(m==i,i))))); 
        end; 
        
        varargout={B};
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
    case 'fit_model_family_singlemodel'                                     % fit single model
        glm         = 3;
        logtaus     = [-1.5:0.3:1.5]; % parameter which determines the number of model families
        regions     = [2];
        hemi        = [1];
        modelTerms  = [2 4];
        vararginoptions(varargin(:),{'glm', 'logtaus','regions','hemi','modelTerms'});
        
        % Get model family
        [ModelFamily, idx]  = sh1_getRDMTau('getfamily','logtaus',logtaus,'modelTerms',modelTerms,'chunkset',[1 2]);
        ModelFamily.prior   = 'Ridge';%'Zellner';
        ModelFamily.X       = ModelFamily.RDM;        
        ModelFamily.Nmodel  = size(ModelFamily.X,1);
        ModelFamily.numComp = size(ModelFamily.X,1);        
        ModelFamily.idx     = idx;
        ModelFamily.numPrior = size(ModelFamily.RDM,1);
        
        % load data set
        T = load(fullfile(regDir,sprintf('distances_sepPerm.glm%d.mat',glm)));
        
        % initialise
        T.omega_hat = zeros(length(T.region),ModelFamily.numComp);
        T.omega_hatn = T.omega_hat;
        T.logEvidence = T.omega_hat;
        T.logtheta = T.omega_hat;
        
        %regions=unique(T.region);        
        for r = regions
            fprintf('Estimating ROI distance %s...\n',regname{r});
            for h = hemi                                                
                indx = find(T.region==r & T.hemis==h);
                D = getrow(T,indx);
                numSubj = length(unique(D.subj));
                for s=1:numSubj
                    Sigma(:,:,s)=reshape(D.Sigma(s,:),8,8);
                    sig(s) = mean(diag(reshape(D.Sigma(s,:),8,8)));
                    X(:,:,s) = ModelFamily.X(:,:,chunk_set(D.subj(s)));
                end;
                F = ModelFamily;
                F.X = X;
                
                % get logEvidence for null model                 
                logEvidence0 = rsa.stat.fitModelNull(D.RDM,sig,9,D.effVox);
                                
                % fit single model from model families
                for m = 1:ModelFamily.Nmodel
                    fprintf('   - Fitting model %s...',F.name{m})
                    F.X = X(m,:,:);
                    %F.X = F.X/sqrt(mean(F.X.^2));
                    
                    [T.omega_hat(indx,m),T.logEvidence(indx,m),theta]=...,S.logEvidenceSplit(indx,:)] = ...
                        rsa_fitModelHierarchEB(F,D.RDM,Sigma,9,D.effVox,'minimiseLength',10);
                    
                    % Summary result
                    T.logtheta(indx,m) = kron(ones(numSubj,1),theta(1));
                    
                    % post-process of omega to adjust by mean of regressor
                    normX = squeeze(sqrt(mean(F.X.^2,2)));
                    T.omega_hatn(indx,m) = bsxfun(@times,T.omega_hat(indx,m),normX);
                    fprintf('... done\n')
                end
                
                % adjust by logEvidence0
                T.logEvidence(indx,:) = T.logEvidence(indx,:)-repmat(logEvidence0,1,F.numComp); % difference from null model
                
            end            
        end
                
        for r=[7 8]
            for s=[2:15]
            subset = T.region==r&T.subj==s;
            
                figure('name',[what,'-',regname{r}]);
                subplot(1,2,1)
                barplot([],T.logEvidence(:,F.idx'),'barwidth',0.8,'subset',subset&T.hemis==1);
                set(gca,'xticklabel',ModelFamily.name(ModelFamily.idx),'view',[270 90],'fontsize',12);
                ylabel('logEvidence','fontsize',12);
                title(sprintf('%s-%s',regname{r},hemName{1}),'fontsize',13)
                
                subplot(1,2,2)
                barplot([],T.logEvidence(:,F.idx'),'barwidth',0.8,'subset',subset&T.hemis==2);
                set(gca,'xticklabel',ModelFamily.name(ModelFamily.idx),'view',[90 270],'fontsize',12);
                ylabel('logEvidence','fontsize',12);
                title(sprintf('%s-%s',regname{r},hemName{2}),'fontsize',13)
            end
        end
                        
        varargout = {T};
    
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
    case 'plot_reliability_btw'                                             % plot between subject reliabitliy of RDM
        C = varargin{1};
        region = [1:8];
        barwidth = 0.8;
                        
        Nsubj = length(unique(C.subj));
        col   = jet(Nsubj);
        for i=1:Nsubj
            Col{i} = col(i,:);
        end
        
        T = [];
        Nregion = numel(region);
        figure('name',sprintf('Between subject reliability of RDM'),'color','w');
        for reg = region;
            fprintf('%s\n',regname{reg})
            count=0;
            for h=[2 1]
                % right hemi
                subplot(2,Nregion,reg+count*Nregion);
                subset = C.region==reg&C.hemis==h;
                upperbound  = nanmean(C.within_u(subset,:));
                se          = nanmean(C.within_u(subset,:))/sqrt(nancount(C.within_u(subset,:)));
                
                % bar to get xpos
                xpos = barplot([],[C.within_l,C.across],'subset',subset,'barwidth',barwidth); hold off
                                
                % ttest
                [t.tval, t.p] = ttest_mc(C.within_l,C.across,1,'paired','subset',subset);
                
                % ranksum test
                [t.prank] = ranksum(C.within_l(subset,:),C.across(subset,:),'tail','both','method','exact');
               
                % drawlines
                drawline(0,'dir','horz','lim',[xpos(1)-barwidth,xpos(2)+barwidth]); hold on
                drawline(upperbound,'dir','horz','color',[1 0.5 0.5],'lim',...
                    [xpos(1)-barwidth,xpos(2)+barwidth],'linewidth',2,'linestyle','--')
                
                % bar
                xpos = barplot([],[C.within_l,C.across],'subset',subset,'barwidth',barwidth); hold on
                
                % line plot
                plotdata.x = [xpos(1)*ones(size(C.within_l));xpos(2)*ones(size(C.across))];
                plotdata.y = [C.within_l;C.across];
                plotdata.s = [C.subj;C.subj];
                lineplot(plotdata.x,plotdata.y,'subset',[subset;subset],'split',plotdata.s,...
                    'linecolor',Col,'markerfill',Col,'markercolor',Col,'linewidth',1.5,'markersize',5);
                
                
                title([regname{reg},'-',hemName{h}],'fontsize',12);
                ylabel('Pearson r','fontsize',12);xlabel('');
                set(gca,'xticklabel',{'W','B'},'fontsize',12)
                %set(gca,'ylim',[-1 1],)
                if t.p<0.05
                    drawline(0.78,'dir','horz','lim',xpos([1,2]),'color',[1 0 0]);
                elseif t.p<0.01
                    drawline(0.78,'dir','horz','lim',xpos([1,2]),'color',[0 1 0]);
                elseif t.p<0.001
                    drawline(0.78,'dir','horz','lim',xpos([1,2]),'color',[0 0 1]);
                end
                count=count+1;
                
                t.hemi = h;
                t.region = reg;
                T = addstruct(T,t);
            end
        end
        
        % t and p values
        figure;
        subplot(3,1,1);
        barplot([T.hemi T.region],T.tval);
        title('t-value')
        
        subplot(3,1,2);
        barplot([T.hemi T.region],log10(T.p));
        drawline(log10(0.05),'dir','horz','color',[1 0 0])
        title('p-value by t-test')
        
        subplot(3,1,3);
        barplot([T.hemi T.region],log10(T.prank));
        drawline(log10(0.05),'dir','horz','color',[1 0 0])
        title('p-value by ranksum test')
        
        varargout = {T};
    case 'plot_omega_new'
        
        T = varargin{1};
        ceiling = [];
        doSQRT  = 1; 
        regions = [1 2 3 4 5 6 7 8 9 10 11 12];
        reguse = [1 2 4 6];
          
        vararginoptions(varargin(2:end),{'ceiling','doSQRT','regions','reguse'});
        
        if ~isempty(ceiling)
            if doSQRT==1
                C = ssqrt(ceiling.ceiling);
            else
                C = ceiling.ceiling;
            end
        end
        
        % ylim setting
        data = T.omega;
        if (doSQRT) 
            D = ssqrt(data); 
            ylim1 = ssqrt([-1.5 30]*0.001); 
            ylim2 = ssqrt([-1.5 30]*0.001); 
            ylim3 = ssqrt([-1.5 30]*0.001); 
        else 
            D = data;
            ylim1 = [-1.5 30]*0.001; 
            ylim2 = [-1.5 30]*0.001;
            ylim3 = [-1.5 30]*0.001; 
        end;
        
        Stats = [];
        % plot
        barwidth = 0.75;
        figure('name',sprintf('EB lin regress (%s)',T.method{1}),'position',[50 50 900/sqrt(2)*0.8 700],'color','w')
        plots = 1;
        for r = 1:length(regions) 
            if ((plots>8))
                figure('name',sprintf('EB lin regress (%s)',T.method{1}),'position',[50 50 900/sqrt(2)*0.8 700],'color','w')
                plots = 1;
            end
            subplot(4,2,plots); 
            plots=plots+1;
                        
            CAT.edgecolor = [1 1 1];
            [x,y,err] = barplot(T.hemis,D,'subset',T.region==regions(r),...
                'CAT',CAT,'barwidth',barwidth); hold on;
            
            midpoint = x(numel(reguse)) + (x(numel(reguse)+1)-x(numel(reguse)))/2;
            
            % draw noise ceiling if passed
            if ~isempty(ceiling)
                tmpL = nanmean(C(T.region==regions(r)&T.hemis==1,:),1); % LH
                h=fill([x(1)-barwidth/2 x(4)+barwidth/2 x(4)+barwidth/2 x(1)-barwidth/2],...
                    [tmpL(2) tmpL(2) tmpL(1) tmpL(1)],[0.8 0.8 0.8]); hold on                
                %set(h,'edgecolor','non','facealpha',0.5);
                tmpR = nanmean(C(T.region==regions(r)&T.hemis==2,:),1); % RH
                h=fill([x(5)-barwidth/2 x(8)+barwidth/2 x(8)+barwidth/2 x(5)-barwidth/2],...
                    [tmpR(2) tmpR(2) tmpR(1) tmpR(1)],[0.8 0.8 0.8]); hold on
                %set(h,'edgecolor','non','facealpha',0.5);
            end
            ylim1(2) = 0.28;%max([tmpL,tmpR]);
            ylim2(2) = 0.28;%max([tmpL,tmpR]);
            ylim3(2) = 0.1;%max([tmpL,tmpR]);
            
            
            title(sprintf('LH<-%s->RH',regname{regions(r)}),'fontsize',15);
            ylabel('Regression weight','fontsize',15)
            set(gca,'xlim',[x(1)-1, x(end)+1],'XTicklabel',{},'fontsize',15);
            set(gca,'fontsize',15,'tickdir','out','ticklength',[0.025 0.025])
                
            % choosing ylim
            if (r<3) 
                set(gca,'YLim',ylim1); 
                drawline(midpoint, 'dir','vert','linestyle','--','lim',ylim1);
            elseif r>8
                set(gca,'YLim',ylim3); 
                drawline(midpoint, 'dir','vert','linestyle','--','lim',ylim3);
            else 
                set(gca,'YLim',ylim2); 
                drawline(midpoint, 'dir','vert','linestyle','--','lim',ylim2);
            end;
            set(gca,'ytick',[0 0.1 0.2 0.25 0.3]);
            
            % one-sample t-test
            for s=1:2 
                for co = 1:size(T.omega,2) 
                    Y = pivottable(T.subj,[],D(:,co),'mean','subset',T.hemis==s & T.region==regions(r));
                    [t,p] = ttest_mc(Y,[],1,'onesample'); 
                    indx = co+(s-1)*length(reguse); 
                    if (p<0.001)
                        text(x(indx),y(indx)+1.5*err(indx),'***',...
                            'HorizontalAlignment','center'); 
                    elseif (p<0.05) 
                        text(x(indx),y(indx)+1.5*err(indx),'*',...
                            'HorizontalAlignment','center'); 
                    end;
                    
                    stat.t = t;
                    stat.df= numel(Y);
                    stat.p = p;
                    stat.region = r;
                    stat.hemis = s;
                    stat.omega = co;
                    
                    Stats = addstruct(Stats,stat);                    
                end; 
            end; 
        end;     
        
        pivottable([Stats.region Stats.hemis],Stats.omega,Stats.p,'(x)');
        
        varargtout = {Stats};
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
nSequence = 8;

I = true(nSequence);
I(diag(diag(I))) = false;     
I(seqnum,:) = false;
I(:,seqnum) = false;

RDMout = RDMin(:,squareform(I));


function out = sh1_calcDistancePerm(data,SPM,seqType);
[RDM,Sw,effVox]=rsa_distanceLDCsepPerm(data,SPM,seqType);
out=[RDM(:);Sw(:);effVox];
