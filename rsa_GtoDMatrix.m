function A=rsa_GtoDMatrix(K);
C=pcm_indicatorMatrix('allpairs',[1:K]); 
numDist = size(C,1); 
for i=1:numDist
    weightG=C(i,:)'*C(i,:); 
    A(i,:)=weightG(:)'; 
end;
