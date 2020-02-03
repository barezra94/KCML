function [Mdl, pred, conf] = validateClassifier(mdlPath,selectedMdlIdx,...
    data, gene, annotations, featureNames, functionalTerms, mdlPrefix, K)
%
% An optional step of cross validation to reshuffle training and test samples.
%  Input:
%       - mdlPath: the directory where trained models are saved
%       - selectedMdlIdx: the index of models that pass selection criteria
%       - data: double matrix of the entire data where columns represent phenotypic
%        features and rows represent genetic perturbations
%       - gene: the gene symbol for each row in data matrix
%       - annotations: a matrix where first column indicate gene symbol
%       annotations and second column is the annotation id
%       - featureNames: the name of the features in data matrix 
%        (should have the same second dimension as data matrix)
%       - functionalTerms: the list of functional terms ids
%       - functionalTerms: the list of functional term names corresponding
%       to those in functionalTerms list
%       - mdlPrefix: the prefix that is used when models are saved
%       -K: number of corss validation
%
%  Output
%       - Mdl: an array of mdl objects containing trained models for all
%         functional terms and training performance
%       - Feats: the selected features for all functional terms
%       - stats: the performance of all models on training and test samples
%       - pred: a matrix of KCML predictions. Each column represent the 
%         prediction that a gene perturbation belong to a certain functional term  
%       - conf: a matrix of KCML prediction confidence. Each column represent the 
%         prediction confidence that a gene perturbation belong to a certain functional term 
%       - prank: a matrix of KCML prediction ranks. Each column represent the 
%         prediction rank where rank 1 has the highest confidence
%
% Copyright (c) Heba Sailem 2019

% Get the names of existing models
fl = dir(strcat(mdlPath,mdlPrefix,'*'));
fl = struct2cell(fl)';
fl = fl(:,1);


% Get model number from the file name
modelNo = strrep(fl,'.mat','');
modelNo = strrep(modelNo,mdlPrefix,'');

fl = strcat(mdlPath,fl);
pred = zeros(size(data,1),length(fl));
conf = [];
prank_nonExp = zeros(length(gene),length(functionalTerms));
for i=1:length(fl)
    mdl_idx = str2num(cell2mat(modelNo(i)));
    if find(mdl_idx==selectedMdlIdx)
        load(cell2mat(fl(i)));

        if ~isfield(mdl,'recall') % Models are empty if classification failed
            Feats{mdl_idx} = {};
            recall(mdl_idx,:) = [0, 0, 0];
            precession(mdl_idx,:) = [0, 0, 0];
            Fscore(mdl_idx,:) = [0, 0, 0];
            n(mdl_idx,1) = 0; nTP(mdl_idx,1) = 0; nPP(mdl_idx,1) = 0; percPP(mdl_idx,1) = 0;
            myRecall(mdl_idx,1) = 0; auc(mdl_idx,1) = 0;
            
        else
            try
                Mdl{mdl_idx} = mdl;
                Feats{mdl_idx} = sel_feats;
                recall(mdl_idx,:) = [mdl.recall(1), mdl.recall(2),mdl.recall(3)];
            end
            precession(mdl_idx,:) = [mdl.precision(1), mdl.precision(2),mdl.precision(3)];
            Fscore(mdl_idx,:) = [mdl.Fscore(1), mdl.Fscore(2),mdl.Fscore(3)];
            
            [feat_idx loc] = ismember(featureNames,Feats{mdl_idx});
            [~,sif] = sort(loc(feat_idx));
            idx = find(feat_idx);
            idx = idx(sif);
            
            posGenes = annotations(cell2mat(annotations(:,2))==functionalTerms(mdl_idx),1);
            l = double(ismember(gene,posGenes));
            l(l==0)=-1;
            
            if K<1%perform leave one out
                K=sum(l>0);
            end
            
            crossValidationSamples = mycvpartition2(l,K);
            
            [prePred(:,mdl_idx), mdlConf] = predict(Mdl{mdl_idx}.mdl,data(:,idx));
            prePositiveIdx = prePred(:,mdl_idx)>0;
            for j=1:K
                cvmdl{j} = fitcsvm(data(crossValidationSamples(:,j),idx),l(crossValidationSamples(:,j)),...
                    'KernelFunction','rbf',...
                    'KernelScale',Mdl{mdl_idx}.z);
                [subpred(:,j), cc] = predict(cvmdl{j},data(:,idx));
                
                try
                    subconf(:,j) = cc(:,2);
                catch
                    subconf(:,j) = cc(:,1);
                end
            end
            Mdl{mdl_idx}.cvmdl = cvmdl;
            pred(:,mdl_idx) = sum(subpred,2)./K;
            
            % SVM confidence
            conf(:,mdl_idx) = mean(subconf,2);
            
        end
    end
end

