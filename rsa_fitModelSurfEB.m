function rsa_fitModelSurfEB(Model,distFiles,sigmFiles,numPart,nvoxFiles,whiteFiles,pialFiles,outFiles,varargin);
% 
% Fits and empricial Bayes model on the surface 
% 
import rsa.util.*;
import rsa.*;

% Options: Not used right now 
Opt = []; 
Opt.minimizeLength = 10; 
Opt=rsa.getUserOptions(varargin,Opt);

numSubj = length(distFiles); 
if length(sigmFiles)~=numSubj
    error('number of sigmFiles need to be the same as distFiles'); 
end; 

% Map the files 
distV = spm_vol(distFiles); 
sigmV = spm_vol(sigmFiles); 
nvoxV = spm_vol(nvoxFiles); 

% generate the mid surfaces 
for i=1:numSubj 
    P=caret_load(pialFiles{i}); 
    PIAL(:,:,i)=P.data; 
    W=caret_load(whiteFiles{i}); 
    WHITE(:,:,i)=W.data;
end; 
MID = (PIAL+WHITE)/2; 
clear PIAL WHITE; 

% Now read blocks of data and fit hem
numVertex = size(MID,1); 
numDist  = length(distV{1}); 
numSigma = length(sigmV{1});
numReg   = size(Model.RDM,1); % Number of regressors 

stepSize=2000;
RESomega=nan(numVertex,numSubj,4); 
RESnsubj=nan(numVertex,1); 
numDist=length(distV{1}); 
numSigm=length(sigmV{1}); 
for j=1:stepSize:numVertex
    fprintf('reading %d/%d....',j,numVertex); tic();
    
    idx_init    = j;
    idx_end     = min(j+min(stepSize-1,numVertex),numVertex);    
    indx = [idx_init:idx_end];
    
    numPoints = length(indx); 
    DIST=nan(numPoints,numSubj,numDist); 
    SIGM=nan(numPoints,numSubj,numSigma);
    NVOX=nan(numPoints,numSubj);
    for s=1:numSubj
        coords=MID(indx,:,s); 
        [x,y,z]=spmj_affine_transform(coords(:,1),coords(:,2),coords(:,3),inv(distV{s}(1).mat));
        for n=1:numDist
            DIST(:,s,n)=spm_sample_vol(distV{s}(n),x,y,z,0); 
        end; 
        [x,y,z]=spmj_affine_transform(coords(:,1),coords(:,2),coords(:,3),inv(sigmV{s}(1).mat));
        for n=1:numSigma
            SIGM(:,s,n)=spm_sample_vol(sigmV{s}(n),x,y,z,0); 
        end; 
        [x,y,z]=spmj_affine_transform(coords(:,1),coords(:,2),coords(:,3),inv(nvoxV{s}(1).mat));
        NVOX(:,s)=spm_sample_vol(nvoxV{s},x,y,z,0); 
   end; 
   DIST = permute(DIST,[2 3 1]); 
   SIGM = permute(SIGM,[2 3 1]); 
   NVOX = permute(NVOX,[2 1]); 
   fprintf('fitting....'); 
   for n=1:length(indx);
       check=sum(abs(DIST(:,:,n)),2);
       in=~isnan(check) & check~=0; 
       RESnsubj(indx(n),1)=sum(in); 
       if (RESnsubj(indx(n),1)>8)
        M=Model; 
         M.RDM=Model.RDM(:,:,in); 
        A=rsa_fitModelHierarchEB(M,DIST(in,:,n),SIGM(in,:,n),numPart,NVOX(in,n),'minimizeLength',Opt.minimizeLength);
        RESomega(indx(n),in,:)=permute(A,[3 1 2]);
       end; 
   end; 
   fprintf('done. '); toc();
end; 

% Now write the results to metric files 

for i=1:numReg 
    M=caret_struct('metric','data',RESomega(:,:,i));
    caret_save(outFiles{i},M); 
end; 
M=caret_struct('metric','data',RESnsubj);
caret_save(outFiles{numReg+1},M); 
