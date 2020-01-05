clear
clc
dataPath = './data/';
DATAfile = [dataPath,'/FusionNetwork.txt'];

[PPIA,PPIB,source] = textread(DATAfile, '%n %n %s');
TempNode = [PPIA;PPIB];                                  
UniqueNode = unique(TempNode); 
NumberOfProteins = size(UniqueNode,1);
[~,edge(:,1)] = ismember(PPIA,UniqueNode); 
[~,edge(:,2)] = ismember(PPIB,UniqueNode);
AdjMatrix = zeros(NumberOfProteins); 
for i = 1:length(PPIA) 
    AdjMatrix(edge(i,1),edge(i,2)) = 1;
    AdjMatrix(edge(i,2),edge(i,1)) = 1;
end
%% ****************************GO Score******************************
[GeneID,GO] = xlsread([dataPath,'/GO.xlsx']); 

Genes = unique(GeneID);
NumberOfUniqueGO = size(Genes,1);

GOScore = zeros(size(UniqueNode,1),2);
GOScore(:,1) = UniqueNode;
for i = 1:size(Genes,1)
    Frequency = find(GeneID(:,1)==Genes(i,1));
    NumberOfAnnotation = size(Frequency,1);
    GOScore_Temp = NumberOfAnnotation/NumberOfUniqueGO;
    [~,~,iUniqueNode] = intersect(Genes(i,1),UniqueNode); 
    GOScore(iUniqueNode,2) = GOScore_Temp;
end
Rank_GO = LP(AdjMatrix,0.6,GOScore(:,2));
RankGO = zeros(size(UniqueNode,1),2);
RankGO(:,1) = UniqueNode; 
RankGO(:,2) = Rank_GO;
%% ******************KnownDiseaseGene Score**************************
[Genes_expression_number,~] = xlsread([dataPath,'/KnownDiseaseGenes.xlsx']);  

KnownDiseaseGene = Genes_expression_number(:,1);
Score = Genes_expression_number(:,5);
KDGScore = zeros(size(UniqueNode,1),2);
KDGScore(:,1) = UniqueNode;
for j = 1:size(KnownDiseaseGene,1)
    [~,~,jUniqueNode] = intersect(KnownDiseaseGene(j,1),UniqueNode); 
    KDGScore(jUniqueNode,2) = Score(j,1);   
end
Rank_KDG = LP(AdjMatrix,0.6,KDGScore(:,2));
RankKDG = zeros(size(UniqueNode,1),2);
RankKDG(:,1) = UniqueNode;
RankKDG(:,2) = Rank_KDG;   
%% ****************************PCA and K-means****************************************************
load([dataPath,'/GenesExpressionData.mat'])
GenesID = GenesExpressionData(2:end,1);
label = GenesExpressionData(1,2:end);  
GenesExpressionValue = GenesExpressionData(2:end,2:end);  
[GeneExpressionDataNormalized,~,~] = zscore(GenesExpressionValue);
GeneExpressionDataNormalizedTrans = GeneExpressionDataNormalized';
[COEFF,SCORE,latent] = pca(GeneExpressionDataNormalizedTrans);  
pcaData1 = SCORE(:,1:2);   
GoodSamples = pcaData1(label==0,:);
PoorSamples = pcaData1(label==1,:);
[Idx,Ctrs,SumD] = kmeans(pcaData1,2);
%% **************************t-test Score*********************************************************
NumberOfSample = size(Idx,1);
ClusterNumber = max(Idx);
ClusterLabelInRawData = cell(ClusterNumber,1);
ClusterGenesExpressionValue = cell(ClusterNumber,1);
ScoreOfAllGene = cell(ClusterNumber,1);
rs = 0.6;
NumOfBiomarker = 100;
SortGene = cell(ClusterNumber,1);
for i = 1:ClusterNumber
    ClusterPostion = find(Idx==i);
    ClusterLabelInRawData{i,1} = label(ClusterPostion);
    ClusterGenesExpressionValue{i,1} = GeneExpressionDataNormalized(:,ClusterPostion);
    StatisticScore = Obtain_StatisticScore(ClusterGenesExpressionValue{i,1},ClusterLabelInRawData{i,1},UniqueNode,GenesID);
    Rank_Statistic = LP(AdjMatrix,rs,StatisticScore(:,2));  
    ScoreOfAllGene{i,1} = [UniqueNode,Rank_Statistic]; 
    [SortGene{i,1},~] = sortrows(ScoreOfAllGene{i,1},-2);
end
if length(SortGene{1,1}) == 1
    Biomarkers = unique(SortGene{1,1},'stable');
else 
    Sort = zeros(length(SortGene{1,1}),2);
    SortTemp = [];
    for Num = 1:length(SortGene{1,1})
        iSortGene = 0;
        for j = 2:ClusterNumber
            [~,~,iSortGeneTemp] = intersect(SortGene{1,1}(Num,1),SortGene{j,1}(:,1));
            iSortGene = iSortGene + iSortGeneTemp;
        end
        RankScore = (Num+iSortGene)/2;
        SortTemp = [SortTemp;RankScore];
    end
    Sort(:,1) = SortGene{1,1}(:,1);
    Sort(:,2) = SortTemp;
    SortNew = sortrows(Sort,2);
    Biomarkers = SortNew(:,1); 
end
FinalMarker = Biomarkers(1:NumOfBiomarker,:);
save Biomarker.mat FinalMarker
