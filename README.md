# FUNMarker
A FUsion Network-based method (FUNMarker) to identify prognostic and heterogeneous biomarkers for accurate survival prediction of breast cancer outcomes

1) Main.m is the main script to call FUNMarker.

2) There are some Matlab scripts for each step of FUNMarker analysis, and called in Main.m
   1. Obtain_StatisticScore.m
   2. LP.m
    
3) The input datasets include:
   1. FusionNetwork.txt
   2. GO.xlsx
   3. KnownDiseaseGenes.xlsx
   4. GenesExpressionData.mat
    
4) The final biomarker list is saved as Biomarker.mat

5) As a demo, users can directly run Main.m in Matlab. This package has been tested in different computer environments as: Window 7 or above; Matlab 2016 or above.
