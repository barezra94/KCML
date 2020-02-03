function [keepin, bestMdl] = selectFeatures(x1,y1,x2,y2,crossValidationSamples, sigma)
%
% Input:
% - x1: training data
% - x2: test data
% - y1: the labels for training samples (binary [-1,1])
% - y2: the labels for test samples (binary [-1,1])
% - crossValidationSamples: is a matrix of the same length as Y where every 
%   column is the index of one cross validation sample.
% - sigma: an array of gamma parameters for optimising SVM RBF kernel
%
% Output:
% - keepin: a logical vector of the index of selected features
% - bestMdl: a structure of trained model, andmodel statistics
% Copyright (c) Heba Sailem 2018
%

%Add features one at time and keep if it decrease the error
bestCrit = 0;

%Forward pass
keepin = logical(zeros(size(x1,2),1));
keepin(1) = 1;

bestMdl = struct();
for i=2:size(x1,2)
    feat_idx = keepin;
    feat_idx(i) = 1;    
    Mdl = optimizeClassifier(x1(:,feat_idx), y1, x2(:,feat_idx), y2, crossValidationSamples, sigma);
    if bestCrit < Mdl.performance % add only features that improve the model
        keepin(i) = 1;
        bestCrit = Mdl.performance;
        bestMdl = Mdl;
        if sum(keepin)>=100
            return;
        end
    end
end

