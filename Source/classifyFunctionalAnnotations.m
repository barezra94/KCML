function [selfeats, Mdl] = classifyFunctionalAnnotations(data, featurenames, crossValidationSamples, opts)
% Input:
% - data: a structure with the following parameters
%       - x1: training data
%       - x2: test data
%       - y1: the labels for training samples (binary [-1,1])
%       - y2: the labels for test samples (binary [-1,1])
%       - x1Idx: the indices for training samples
%       - x2Idx: the indices for test samples
% - featureNames: the name of the features in data matrix.
% - crossValidationSamples:  is a matrix of the same length as x1 where every 
%   column is the index of one cross validation sample.
% - opts: a structure with the following parameters
%           - sigma:  an array of gamma parameters for optimising SVM RBF kernel
%           - maxNumFeats: the maximum number of features to be selected 
%             for a classifier
%
% Output:
% - selFeatus: the list of selected features
% - Mdl: the trained model, and model statistics
%
% Copyright (c) Heba Sailem 2018
%

if nargin <4
    opts.maxNumFeat = size(data.x1,2);
    opts.sigma = 1;
end

%data
    x1 = data.x1;
    x2 = data.x2;

    %
    g1_idx1 = data.y1==-1;
    g1_idx2 = data.y2==-1;

    g2_idx1 = data.y1==1;
    g2_idx2 = data.y2==1;
    p=[];
    
    for i=1:size(x1,2)
        [h p(i)] = kstest2([x1(g1_idx1,i);x2(g1_idx2,i)],[x1(g2_idx1,i);x2(g2_idx2,i)]);
    end

    [s si] = sort(p);
    [sel, Mdl] = selectFeatures(x1(:,si(1:opts.maxNumFeat)),data.y1,x2(:,si(1:opts.maxNumFeat)),data.y2,crossValidationSamples, opts.sigma);
    sfeaturenames = featurenames(si(1:opts.maxNumFeat));
    selfeats = sfeaturenames(sel);
end
