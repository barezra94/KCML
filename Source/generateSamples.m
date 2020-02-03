function generateSamples(samplesPath,datafile,annotationfile,negativeSampleSize,numberCrossValidations,termIdx)
% Generate and save samples for training and testing KCML
%
% Copyright (c) Heba Sailem 2018
load(datafile);
load(annotationfile);

if nargin<6
    termIdx = 1:length(functionalTerms);
end

if ~exist(samplesPath,'file')
    mkdir(samplesPath)
end

annotmat = cell2mat(annotations(:,2));

for ii=1:length(termIdx)
    posIdx = ismember(gene,annotations(annotmat==functionalTerms(ii),1));
    negativeGeneList = gene(~posIdx); % Here one can exclude genes in related terms
    
    dSample = sampleFunctionalAnnotation(functionalTerms(ii),annotations,data,gene,negativeSampleSize,negativeGeneList);    
    
    % For training sample, obtain many cross validation sample
    crossValidationSamples = mycvpartition(dSample.y1,numberCrossValidations);
    dIdx.x1Idx = dSample.x1Idx;
    dIdx.x2Idx = dSample.x2Idx;
    dIdx.y1 = dSample.y1;
    dIdx.y2 = dSample.y2;
    
    % Only the index of samples are saved
    save(strcat(samplesPath,'/Sample',num2str(functionalTerms(ii))),'dIdx','crossValidationSamples');
end




