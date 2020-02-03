function [Sample]=sampleFunctionalAnnotation(GOTerm, annotations, data, geneSymbols,negSampleSize,negativeGeneList)
% Genrate random samples from the data matrix for a given GOTerm
%
% Copyright (c) Heba Sailem 2019

opts.RandomSeed=1;
% Retrieve annotations for the queried term
idx=find(cell2mat(annotations(:,2))==GOTerm);
positiveGeneList = geneSymbols(ismember(geneSymbols,annotations(idx,1)));

% Ensure that no positive genes in the negative list
negativeGeneList = negativeGeneList(~ismember(negativeGeneList,positiveGeneList));



% Use 70% of positive samples for training
percent_train = 0.7;
npos = length(positiveGeneList);

coef = round(negSampleSize./length(positiveGeneList));%obtain large list of negative samples
nneg = min([length(positiveGeneList)*coef;length(negativeGeneList)]);


%select the negative sample based on significance 

ssize1 = floor(npos*percent_train);
ssize2 = floor(nneg*percent_train);

% Set the random seed to guarantee reproducibility
rng(opts.RandomSeed);
shuffled_idx1 = randperm(npos);
rng(opts.RandomSeed);
shuffled_idx2 = randperm(length(negativeGeneList));

train_gn1 = positiveGeneList(shuffled_idx1(1:ssize1));
test_gn1 = positiveGeneList(shuffled_idx1(ssize1+1:npos));

if nneg<=length(negativeGeneList)
    train_gn2 = negativeGeneList(1:ssize2);
    test_gn2 = negativeGeneList(ssize2+1:nneg);
else
    negSample=[negativeGeneList;rem_ls2(shuffled_idx2(1:nneg-length(negativeGeneList)))];
    train_gn2 = negSample(1:ssize2);
    test_gn2 = negSample(ssize2+1:nneg);
end

genelist.pos=positiveGeneList;
genelist.neg=[train_gn1;test_gn2];

train_idx1 = find(ismember(geneSymbols,train_gn1));
train_idx2 = find(ismember(geneSymbols,train_gn2));
dataTrainIdx = [train_idx1; train_idx2];

test_idx1 = find(ismember(geneSymbols,test_gn1));
test_idx2 = find(ismember(geneSymbols,test_gn2));
dataTestIdx = [test_idx1; test_idx2];

dataTrainG1 = data(train_idx1,:);
dataTrainG2 = data(train_idx2,:);
dataTrain = [dataTrainG1; dataTrainG2];

yTrain1=repmat(1,size(dataTrainG1,1),1);
yTrain2=repmat(-1,size(dataTrainG2,1),1);
yTrain=[yTrain1; yTrain2];

dataTestG1 = data(test_idx1,:);
dataTestG2 = data(test_idx2,:);
dataTest = [dataTestG1; dataTestG2];

yTest1 = repmat(1,size(dataTestG1,1),1);
yTest2 = repmat(-1,size(dataTestG2,1),1);
yTest = [yTest1; yTest2];

Sample.x1 = dataTrain;
Sample.x1Idx = dataTrainIdx;
Sample.x2 = dataTest;
Sample.x2Idx = dataTestIdx;
Sample.y1 = yTrain;
Sample.y2 = yTest;