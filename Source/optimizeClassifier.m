function [Mdl]=optimizeClassifier(X,Y,x2,y2,crossValidationSamples,sigma)
% 
% Input:
% - X: training data
% - x2: test data
% - Y: the labels for training samples (binary [-1,1])
% - y2: the labels for test samples (binary [-1,1])
% - crossValidationSamples: is a matrix of the same length as Y where every 
%   column is the index of one cross validation sample.
% - sigma: an array of gamma parameters for optimising SVM RBF kernel
%
% Output:
% Mdl: a structure of trained model, andmodel statistics
%
% Copyright (c) Heba Sailem 2018
%


rng('default') % For reproducibility

for j=1:length(sigma)
    parfor i=1:size(crossValidationSamples,2)
        mdl{i,j} = fitcsvm(X(crossValidationSamples(:,i),:), Y(crossValidationSamples(:,i)),...
            'KernelFunction','rbf',...
            'KernelScale',sigma(j));
        
        pred1=predict(mdl{i,j},X);
        [Fscore1(i,j), acc1(i,j), recall1(i,j), precision1(i,j),FPR1(i,j)]=getClassificationMetrics(Y,pred1);
        
        pred2=predict(mdl{i,j},x2);
        [Fscore2(i,j), acc2(i,j), recall2(i,j), precision2(i,j),FPR2(i,j)]=getClassificationMetrics(y2,pred2);
      
    end
end

acc = (acc1+acc2)./2;
FPR = (FPR1+FPR2)./2;
Fscore = (Fscore1+Fscore2)./2;
recall = (recall1+recall2)./2;
precision = (precision1+precision2)./2;

[mdlPerform mdlIdx] = max(Fscore2); % select best model based on Score
[bestGamma gammaIdx] = max(mdlPerform);
xIdx = mdlIdx(gammaIdx);
yIdx = gammaIdx;

bestMdl = mdl{xIdx,yIdx};

Mdl.mdl = bestMdl;
Mdl.FPR = FPR(xIdx,yIdx);
Mdl.z = sigma(yIdx);
Mdl.performance = Fscore(xIdx,yIdx);
Mdl.sample = crossValidationSamples(:,xIdx);
Mdl.acc = [acc1(xIdx,yIdx), acc2(xIdx,yIdx), acc(xIdx,yIdx)];
Mdl.Fscore = [Fscore1(xIdx,yIdx), Fscore2(xIdx,yIdx), Fscore(xIdx,yIdx)];
Mdl.recall = [recall1(xIdx,yIdx), recall2(xIdx,yIdx), recall(xIdx,yIdx)];
Mdl.precision = [precision1(xIdx,yIdx), precision2(xIdx,yIdx), precision(xIdx,yIdx)];

