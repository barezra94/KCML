function [Fscore, accuracy, recall, precession, FPR, AUC]=getClassificationMetrics(Y,pred,predScores)
% Compute various classification metrics for a label vector Y and
% predcition vector pred and prediciton confidence predScores
%
% Copyright (c) Heba Sailem 2018
if nargin<3
    predScores = pred;
end
pred = sign(pred);
TP = sum(pred==Y & Y==1);
TN = sum(pred==Y & Y==-1);
FP = sum(pred~=Y & pred==1);
FN = sum(pred~=Y & pred==-1);

accuracy = sum(pred==Y)/length(Y);
recall = TP/(TP+FN);
FPR = FP/(TN+FP);
precession = TP/(TP+FP);
Fscore = (2*recall*precession)./(recall+precession);
%[X,Y,T,auc]=perfcurve(l,predScores,'1');
[~, ~, ~, AUC]=perfcurve(Y, predScores, '1');
