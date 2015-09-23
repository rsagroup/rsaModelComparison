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
        T=load(fullfile(regDir,'distances_sepPerm.mat')); 
        regions=unique(T.region); 
        for r=regions' 
            for h=[1 2] 
                indx = find(T.region==r & T.hemis==h); 
                D=getrow(T,indx); 
                for s=1:length(D.subj)
                    Sigma(:,:,s)=reshape(D.Sigma(s,:),8,8); 
                end; 
                Model.fcn = @sh1_getRDMmodelTau1;
                Model.numComp   = 4; % Number of linear components 
                Model.numPrior  = 4; % Number of prior variances on parameters
                Model.numNonlin = 4; % Number of nonlinear parameters 
                Model.nonlinP0  = [-2 1 0 0]; % Starting value of nonlinear(mixing) parameters 
                Model.constantParams = {chunk_set(D.subj),...       % Model terms 
                                        [1 2 4 6],...               % Chunk set 
                                        'sqEuclidean'...            % Distance term               
                                        };  
                [omega(indx,:),logEvidence(indx,:),lt]=rsa_fitModelHierarchEB(Model,D.RDM,Sigma,9,D.effVox);
                logtheta(indx,:)=repmat(lt,length(indx),1); 
            end; 
        end; 
        T.omega=omega; 
        T.logtheta=logtheta; 
        varargout={T}; 

    otherwise
        disp('there is no such case.')
end;

