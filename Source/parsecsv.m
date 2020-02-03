% parse csv file

% load annotations
annotationsTable = readtable(strcat(annotationFile,'.csv'),'PreserveVariableNames',true,'Delimiter',',');
functionalTermNames = unique(annotationsTable.FunctionalTermName);
functionalTerms = [1:length(functionalTermNames)]'; % assign ids to each functional term

termId=[];
for i = 1:length(functionalTerms)
    idx = strcmpi(functionalTermNames(i),annotationsTable.FunctionalTermName);
    termId(idx,1) = functionalTerms(i);
end

annotations = [ annotationsTable.GeneSymbol, num2cell(termId) ];
save(annotationFile,'annotations','functionalTermNames','functionalTerms')

% load data
dataTable = readtable(strcat(dataFile,'.csv'),'PreserveVariableNames',true,'Delimiter',',');
gene = dataTable.GeneSymbol;
featureNames = dataTable.Properties.VariableNames;
data = table2array(dataTable(:,2:end));
save(dataFile,'data','featureNames','gene');
