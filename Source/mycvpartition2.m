function crossValidationSamples = mycvpartition2(Y,K)
% Generate K random cross validatation samples into training and test
% samples
% Input:
% Y: label matrix
% K: number of cross validation samples

% Output:
% crossValidationSamples: is a matrix of the same length as Y where every 
% column is the index of one cross validation sample.
%
% Copyright (c) Heba Sailem 2019
%
posIdx = find(Y==1);
numPos = length(posIdx);
sampleSize = round((numPos./K)*(K-1));
negIdx = find(Y~=1);
crossValidationSamples = logical(zeros(length(Y),K));
rng(1);
cvIndices = crossvalind('kfold',numPos,K);

for i=1:K
    rng(i);
    negSamples = randperm(length(negIdx), sampleSize)';
    crossValidationSamples(negIdx(negSamples),i) = 1;
    crossValidationSamples(posIdx(cvIndices~=i),i) = 1;
end

