if strcmp(DataSet,'RNAMix1') 
    
    load('RNAMix1Labels.mat'); load('RNAMix1Linnorm.mat'); XLin = X; load('RNAMix1LinnormND.mat'); XLinND = X; load('RNAMix1BasicND.mat'); XND = X; load('RNAMix1Basic.mat'); epsilon = 2.6; epsilon_lin = 1; epsilon_umap = 0.64; %epsilon_umap = 0.7;
    
end

if strcmp(DataSet,'RNAMix2') 

    load('RNAMix2Labels.mat'); load('RNAMix2Linnorm.mat'); XLin = X; load('RNAMix2LinnormND.mat'); XLinND = X; load('RNAMix2BasicND.mat'); XND = X; load('RNAMix2Basic.mat'); epsilon = 1.6; epsilon_lin = .9; epsilon_umap = 0.7; 
    
end

if strcmp(DataSet,'CellMix') 
    
    load('CellMixLabels.mat'); load('CellMixLinnorm.mat'); XLin = X; load('CellMixLinnormND.mat'); XLinND = X; load('CellMixSCTND.mat'); XND = X; load('CellMixSCT.mat'); epsilon = 7; epsilon_lin = 4.05; epsilon_umap = 3.5; 
    % Anna's fisrt epsilon_umap = 3.5;
    % epsilon = 7 for Basic
end

if strcmp(DataSet,'CellMixSng') 
    
    load('CellMixSngLabels.mat'); load('CellMixSngLinnorm.mat'); XLin = X; load('CellMixSngLinnormND.mat'); XLinND = X; load('CellMixSngSCTND.mat'); XND = X; load('CellMixSngSCT.mat'); epsilon = 9; epsilon_lin = 7; epsilon_umap = 5; 
    %Anna's fisrt epsilon_umap = 5;
    % epsilon = 7 for Basic
end

if strcmp(DataSet,'Beta')
    
    load('BetaLabels.mat'); load('BetaLinnorm.mat'); XLin = X; load('BetaLinnormND.mat'); XLinND = X; load('BetaBasicND.mat'); XND = X; load('BetaBasic.mat'); epsilon = 2.45; epsilon_lin = 3.8; epsilon_umap = 0.575;
    
end

if strcmp(DataSet,'Beta2')
    
    load('Beta2Labels.mat'); load('Beta2Linnorm.mat'); XLin = X; load('Beta2LinnormND.mat'); XLinND = X; load('Beta2BasicND.mat'); XND = X; load('Beta2Basic.mat'); epsilon = .0415; epsilon_lin = .045; epsilon_umap = 0.49 ;  %0.52
    
end

if strcmp(DataSet,'BaronPanc')
    
    load('BaronPancLabels.mat'); load('BaronPancLinnorm.mat'); XLin = X; load('BaronPancLinnormND.mat'); XLinND = X; load('BaronPancSCTND.mat'); XND = X; load('BaronPancSCT.mat'); epsilon = 0.0283; epsilon_lin = 14.5; epsilon_umap = 2.5;
    % epsilon = 5.5 for Basic
end

if strcmp(DataSet,'TMPanc')
    
    load('TMPancLabels.mat'); load('TMPancLinnorm.mat'); XLin = X; load('TMPancLinnormND.mat'); XLinND = X; load('TMPancBasicND.mat'); XND = X; load('TMPancBasic.mat'); epsilon = 6; epsilon_lin = 10.3; epsilon_umap = 0.46; 
    
end

if strcmp(DataSet,'TMLung')
    
    load('TMLungLabels.mat'); load('TMLungLinnorm.mat'); XLin = X; load('TMLungLinnormND.mat'); XLinND = X; load('TMLungBasicND.mat'); XND = X; load('TMLungBasic.mat'); epsilon = 10; epsilon_lin = 18; epsilon_umap =0.52; %0.6
    
end

if strcmp(DataSet,'PBMC3k')
    
    load('PBMC3klabels.mat'); load('PBMC3kLinnorm.mat'); XLin = X; load('PBMC3kLinnormND.mat'); XLinND = X; load('PBMC3kSCTND.mat'); XND = X; load('PBMC3kSCT.mat'); epsilon = 2.5; epsilon_lin = 5; epsilon_umap = 0.6;
   
end
if strcmp(DataSet,'PBMC3k_SingleR')
    
    load('PBMC3kSingleRLabels.mat'); load('PBMC3kSingleRLinnorm.mat'); XLin = X; load('PBMC3kSingleRLinnormND.mat'); XLinND = X; load('PBMC3kSingleRSCTND.mat'); XND = X; load('PBMC3kSingleRSCT.mat'); epsilon = 2.5; epsilon_lin = 5; epsilon_umap = 0.5;
 
end
if strcmp(DataSet,'PBMC4k_main')
    
    load('PBMC4kmainlabels_merged_Tcells.mat'); load('PBMC4kMAINLinnorm.mat'); XLin = X; load('PBMC4kMAINLinnormND.mat'); XLinND = X; load('PBMC4kMAINBasicND.mat'); XND = X; load('PBMC4kMAINBasic.mat'); epsilon = 5; epsilon_lin = 5.5; epsilon_umap = 5.7;
                                                                                                                                                                                                            %epsilon = 2.2; epsilon_lin = 3.6; epsilon_umap = 0.65                        
end
if strcmp(DataSet,'PBMC4k_extented')
    
    load('PBMC4kextralabels.mat'); load('PBMC4kEXTRALinnorm.mat'); XLin = X; load('PBMC4kEXTRALinnormND.mat'); XLinND = X; load('PBMC4kEXTRASCTND.mat'); XND = X; load('PBMC4kEXTRASCT.mat'); epsilon = 1.5; epsilon_lin = 2.435;epsilon_umap = 0.325;
    
end