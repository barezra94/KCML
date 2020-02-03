function crossValidationSamples = mycvpartition(Y,K)
% Generate K random cross validatation samples from training data
% Y: label matrix
% K: number of cross validation samples
% crossValidationSamples is a matrix of the same length as Y where every 
% column is the indeces of one cross validation sample the index 
%
% Copyright (c) Heba Sailem 2018
%

posIdx = Y==1;
negIdx = find(Y~=1);
%negSubSampleIdx=floor(linspace(0,length(negIdx),K+1));
crossValidationSamples = logical(zeros(length(Y),K));
for i=1:K
    crNeg = logical(zeros(length(Y),1));
    rng(i);
    negSamples = randperm(length(negIdx), [sum(posIdx)])';
    crNeg(negIdx(negSamples)) = 1;
    crossValidationSamples(:,i) = posIdx | crNeg; 
end
