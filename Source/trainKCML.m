function trainKCML(dataFile, annotationFile, savePath, opts)
%
% The main functionality of KCML
%   
%   1) dataFile:
%       the path and filename of the data file.  dataFile should have
%       the following arrays:
%       data: 
%         a double matrix where each row is a vector of gene perturbation
%         features.
%       gene:
%         An array of gene names or symbols. The number of rows in the gene
%         array is the same as the number of rows in the data matrix   
%       featureNames: 
%         An array of the feature names (columns) of the 'data' matrix.
%
%    2) annoationsFile: 
%       the path and file name of annotation file. annoationsFile should have
%       the following arrays:
%       
%       functionalTerms: 
%         An array of of the annotations that will be classified by KCML 
%       annotations: 
%         a cell array of genes (first column) and their annotations
%       negativeAnnotations: 
%         similar as annotations but specifying the
%         terms that should not be included in the negative list. 
%         This might be due to prior belief that they might be
%         related to the respective functional annotation.
%     3) savePath: the directory for saving generated sample and trained
%        models
%     4) opts: a structure with the following parameters
%           - sigma: an array of gamma parameters for optimising SVM RBF kernel
%           - nCrossValidation: for feature selection and optimising model
%           parameters.
%           - maxNumFeats: the maximum number of features to be selected 
%             for a classifier
%           
% Copyright (c) Heba Sailem 2018
fprintf('Start KCML training\n');
load(dataFile);
load(annotationFile);

% Set KCML parameters here
negativeSampleSize = 5000;
modelPath = fullfile(savePath,'Models/');
samplePath = fullfile(savePath,'Samples/');
maxNumFeats = size(data,2);
nModels = length(functionalTerms);

% Gamma parameter for the RBF kernel
sigma = opts.sigma;


% Directory for saving the trained classifiers
if ~exist(modelPath,'file')
    mkdir(modelPath);
end


% Generate random samples of positive and negative gene profiles for each GO term
% Samples are saved for tracability (this includes the cross validation
% samples.)
% The sample is composed of the index of gene perturbation that should be
% used for trainining (x1) or testing (x2) as well their labels (y1, and y2)
% y=1 if the gene is annotated with the respective functional annotation
% and y=-1 otherwise.

if ~exist(fullfile(savePath,'Samples/Sample1.mat'),'file')
    fprintf('Genrating training samples\n');
    generateSamples(samplePath, dataFile, annotationFile, negativeSampleSize,opts.nCrossValidation);
end

% Apply KCML to generated samples

for ii=1:length(functionalTerms)
    fprintf('Training classifier %d out of %d \n', ii, length(functionalTerms));
    load(strcat(samplePath,'Sample',num2str(functionalTerms(ii))));
    dataSample = dIdx;
    dataSample.x1 = data(dataSample.x1Idx,:);
    dataSample.x2 = data(dataSample.x2Idx,:);
    [sel_feats, mdl] = classifyFunctionalAnnotations(dataSample, featureNames, crossValidationSamples, opts);
    save(strcat(modelPath,'Mdl',num2str(functionalTerms(ii))),'sel_feats','mdl');
end