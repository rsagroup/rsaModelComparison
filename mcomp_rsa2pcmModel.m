function Mpcm=mcomp_rsa2pcmModel(Mrsa) 
% Transforms a set of models from rsa to pcm format 
for i=1:length(Mrsa) 
    Mpcm{i}.name = Mrsa(i).name; 
    Mpcm{i}.type = 'fixed'; 
    Mpcm{i}.numGparams=0; 
    RDM = squareform(Mrsa(i).RDM); 
    numCond = size(RDM,1); 
    H = eye(numCond)-ones(numCond)/numCond; 
    Mpcm{i}.Gc = -0.5 * H*RDM*H';
    Mpcm{i}.RDM=Mrsa(i).RDM; 
end; 
