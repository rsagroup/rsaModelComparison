function Model=sh1_getRDMmodelTau1(logtau,chunkset,components,distfun);
% function Model=sh1_getRDMmodelTau1(logtau,chunkset,type,distfun);
% Returns a RDM model structure using a number of different components
% INPUT:
%   theta: Log of Decay time constant y=exp(-1/exp(theta))
%          time of 1 is the total sequence
%          time of 1/4 is the time constant of a chunk
%          time of 1/11 is a time constant of a press
%   chunkset:
%          A:1 chunkset B:2
%   components:
%          onedigit       (1)
%          twodigit       (2)
%          threedigit     (3)
%          chunk          (4)
%          superchunk     (5)
%          sequence       (6)
%   distfun:
%          Distance function for activtion (cityblock, sqEucledian)
%          or normSqEucledian 
% OUTPUT:
%   Model.name:         Name of the components (Kx1) cell array
%   Model.RDM:          Predicted distances under that component and tau
%   Model.dRMDdtheta:   Derivative of the component in respect to each of
%                       the thetas
% Joern Diedrichsen, Atsushi Yokoi
persistent onedigit twodigit threedigit chunk superchunk sequence ChunkFing Con Features modelNames;
numSeq    = 8;
numChunk  = 8;
seqLength = 11;
numDist   = numSeq*(numSeq-1)/2;


% -----------------------------------------------
% First time into the function: execute this code
if (isempty(ChunkFing))
    modelNames = {'onedigit','twodigit','threedigit','chunk','superchunk','sequence'};
    normalisationFact = [1 1.8 2 1.4 1 0.8];
    
    ChunkFing  = {{[1 3];[5 2 4];[2 3 2];[5 1 4];
        [3 5];[4 2 1];[2 5 2];[1 4 3]};
        {[1 4];[3 2 5];[2 4 2];[1 3 5];
        [5 1];[4 2 3];[2 5 2];[4 1 3]}};
    chunk = {[1 2 3 4;4 3 2 1;5 6 7 8;8 7 6 5;...
        1 2 7 8;8 7 2 1;5 6 3 4;4 3 6 5];...
        [4 3 2 1;5 6 7 8;8 7 6 5;1 2 3 4;
        4 3 6 5;1 2 7 8;8 7 2 1;5 6 3 4];};
    
    % Make single fingers: For each sequence set onedigit is a
    % 8 x 11 matrix of indcators for fingers
    for s=1:2
        for i=1:numSeq;
            A=ChunkFing{s}(chunk{s}(i,:));
            AA = [];
            for j = 1:numel(A)
                AA = [AA, A{j}];
            end
            onedigit{s}(i,:)=AA;
        end;
    end;
    
    % Make single fingers: For each sequence set twodigit is a
    % 8 x 10 matrix of indcators for finger transitions
    for s=1:2
        for i=1:8
            for p=1:10
                twodigit{s}(i,p)=(onedigit{s}(i,p)-1)*5+onedigit{s}(i,p+1); % Unique 2-chunk
            end;
        end;
    end;
    
    % Make 3finger transition
    for s=1:2
        for i=1:8
            for p=1:10
                threedigit{s}(i,p)=(onedigit{s}(i,p)-1)*25+(onedigit{s}(i,p+1)-1)*5+onedigit{s}(i,p+1); % Unique 2-chunk
            end;
        end;
    end;
    
    % Make Superchunk
    for s=1:2       % Sequence sets
        for i = 1:8     % Sequence
            for p = 1:2
                superchunk{s}(i,p) = (chunk{s}(i,p*2-1)-1)*numChunk+chunk{s}(i,p*2); % Unique 2-chunk
            end;
        end;
    end;
    
    % Make sequence
    for s=1:2       % Sequence sets
        sequence{s} = [1:8]';
    end;
    
    % Contrast matrix
    Con=indicatorMatrix('allpairs',[1:8]);
    
    % Now make pre-make Feature matrices for all terms
    for s=1:2 % Loop over chunk sets
        for m=1:length(modelNames)  % Different model terms
            variable     = eval([modelNames{m} '{s};']);
            numTimep     = size(variable,2);     % Number of time points
            maxType      = max(variable(:));   % Maximal number of types
            
            for t=1:numTimep
                indx1        = [1:numSeq]';                    % Index for the sequence
                indx2        = variable(:,t);                  % Index for the variable
                Features{s,m}(:,:,t) = full(sparse(indx1,indx2,1,numSeq,maxType))./normalisationFact(m);
            end;
        end;
    end;
end;

numComp   = length(components);
numParam  = length(logtau);
numSubj   = length(chunkset);
chunksets = unique(chunkset);

for i=1:numComp
    for c = chunksets   % precalculate all necessary chunksets 
        Feature       = Features{c,components(i)};
        numTimep      = size(Feature,3);                                    % Number of time points
        t             = [0:numTimep-1]/numTimep;
        tempWeight    = exp(-t*exp(-logtau(i)));                            % Calculate weighting factor        
        dTempWeight   = tempWeight.*t.*exp(-logtau(i));                     % This is a derivative wrt logtau
        wFeature      = bsxfun(@times,Feature,permute(tempWeight,[1 3 2]));          % Multiply every time slice with the weight
        Pattern       = sum(wFeature,3); % Integrate over time
        dwFeature      = bsxfun(@times,Feature,permute(dTempWeight,[1 3 2]));          % Multiply every time slice with the weight
        dPattern       = sum(dwFeature,3); % Integrate over time
        
        % Calculate distance - the following is much faster than pdist
        Model.name{i,1} = modelNames{components(i)};
        switch (distfun)
            case 'cityblock'
                Model.RDM(i,:)=sum(abs(Con*Pattern),2)';
            case 'sqEuclidean'
                Model.RDM(i,:,c)          = sum((Con*Pattern).^2,2)';
                Model.dRDMdTheta(:,:,i,c) = zeros(numParam,numDist);                
                Model.dRDMdTheta(i,:,i,c) = sum(2*(Con*dPattern).*(Con*Pattern),2); % Fix: corrected derivative
            case 'normSqEuclidean'
                Model.RDM(i,:,c)          = sum((Con*Pattern).^2,2)';
                Model.RDM(i,:,c)          = bsxfun(@rdivide,Model.RDM(i,:,c),sqrt(mean(Model.RDM(i,:,c).^2)));
            otherwise
                error('wrong parameter for distfun');
        end;
    end;
end;

% Expand the models for the different subjects 
Model.RDM = Model.RDM(:,:,chunkset); 
if strcmp(distfun,'sqEuclidean')
    Model.dRDMdTheta = Model.dRDMdTheta(:,:,:,chunkset);    
end
