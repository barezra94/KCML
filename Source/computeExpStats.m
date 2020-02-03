function [Mdl, Feats, stats, pred, conf, prank]=computeExpStats(mdlPath,...
    data, gene, annotations,featureNames,functionalTerms,mdlPrefix)
% Estimate KCML performance, and select successful classifiers.
%  Input:
%       - mdlPath: the directory where trained models are saved
%       - data: double matrix of the entire data where columns represent phenotypic
%        features and rows represent genetic perturbations
%       - gene: the gene symbol for each row in data matrix
%       - annotations: a matrix where first column indicate gene symbol
%       annotations and second column is the annotation id
%       - featureNames: the name of the features in data matrix 
%        (should have the same second dimension)
%       - functionalTerms: the list of functional terms ids
%       - functionalTerms: the list of functional term names corresponding
%       to those in functionalTerms list
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
% Copyright (c) Heba Sailem 2018
%

% Get the names of existing models
fl = dir(strcat(mdlPath,mdlPrefix,'*'));
fl = struct2cell(fl)';
fl = fl(:,1);


% Get model number from the file name
modelNo = strrep(fl,'.mat','');
modelNo = strrep(modelNo,mdlPrefix,'');

fl = strcat(mdlPath,fl);
pred = [];
conf = [];

for i=1:length(fl)
    mdl_idx = str2num(cell2mat(modelNo(i)));
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

        
        [pred(:,mdl_idx), mdlConf] = predict(Mdl{mdl_idx}.mdl,data(:,idx));
        try
            conf(:,mdl_idx) = mdlConf(:,2);
        catch
            conf(:,mdl_idx) = mdlConf(:,1);
        end
        
        
        [~,si] = sort(conf(:,mdl_idx),'Descend');
        [~,r] = sort(si);
        prank(:,mdl_idx) = r;
        
        pos_idx = pred(:,mdl_idx)>0;
        posGenes = annotations(cell2mat(annotations(:,2))==mdl_idx,1);
        l = double(ismember(gene,posGenes));
        [~,~,~,auc(mdl_idx,1)] = perfcurve(l,pred(:,mdl_idx),'1');
        
        n(mdl_idx,1) = sum(l);%considering number of genes in the dataset
        nTP(mdl_idx,1) = sum(ismember(gene(pos_idx),posGenes));
        
        %num predicted positive
        nPP(mdl_idx,1) = sum(pos_idx)-nTP(mdl_idx,1);
        percPP(mdl_idx,1) = nPP(mdl_idx,1)./size(data,1);
        myRecall(mdl_idx,1) = sum(ismember(gene(pos_idx),posGenes))/length(posGenes);
    end
end

stats = array2table([n, nTP, nPP, percPP, myRecall(:,1), recall(:,1),recall(:,2), recall(:,3), ...
    precession(:,1), precession(:,2), precession(:,3),...
    Fscore(:,1), Fscore(:,2), (Fscore(:,1)+Fscore(:,2))/2, auc],...
    'VariableNames',{'N','numberTruePositives','numberPredictedPositives',...
    'percentPredictedPositive', 'overallRecall','trainingRecall','testRecall',....
    'subRecall','trainingPrecession','testPrecession','overallPrecession',...
    'trainingFscore','testFscore','overallFscore','AUC'});

