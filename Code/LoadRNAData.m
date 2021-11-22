if strcmp(DataSet,'RNAMix1')
    
    if strcmp(Normalization,'Basic')
        
        load('RNAmix1original.mat');Data=RNAmix1original';load('rnamix1originallabels.mat');Labels=rnamix1originallabels;
        Data = diag(1./sum(Data,2))*Data*10000;Data = log(Data+1); DataOpts.VarCutoff=.3;
        
    elseif strcmp(Normalization,'Linnorm')
        
        load('RNAmixSortSeq.mat'); Data = RNAmixSortSeq'; Labels = RNAmixSortSeqLabels; mix1_idx=find(Labels==1); Labels(mix1_idx) = 2; Labels=Labels-1;
        
    elseif strcmp(Normalization,'Seurat')
        
        load('rnamix1seuratscaled.mat');Data=rnamix1seuratscaled';load('rnamix1originallabels.mat');Labels=rnamix1originallabels;
    elseif strcmp(Normalization,'SCT')
        
        load('rnamix1SCT.mat');Data=rnamix1SCT';load('rnamix1originallabels.mat');Labels=rnamix1originallabels;
        
    end
    
end

if strcmp(DataSet,'RNAMix2')
    
    if strcmp(Normalization,'Basic')
        
        load('RNAmix2original.mat');Data=RNAmix2original';load('rnamix2originallabels.mat');Labels=rnamix2originallabels;
        Data = diag(1./sum(Data,2))*Data*10000;Data = log(Data+1);DataOpts.VarCutoff=.2;
        
    elseif strcmp(Normalization,'Linnorm')
        
        load('RNAmixCELSeq2.mat'); Data = RNAmixCELSeq2'; Labels = RNAmixCELSeq2Labels; mix1_idx=find(Labels==1); Labels(mix1_idx) = 2; Labels=Labels-1;
        
    elseif strcmp(Normalization,'Seurat')
        
        load('rnamix2seuratscaled.mat');Data=rnamix2seuratscaled';load('rnamix2originallabels.mat');Labels=rnamix2originallabels;
    elseif strcmp(Normalization,'SCT')
        
        load('rnamix2SCT.mat');Data=rnamix2SCT';load('rnamix2originallabels.mat');Labels=rnamix2originallabels;
      
    end
    
end

if strcmp(DataSet,'CellMix')
    
    if strcmp(Normalization,'Basic')
        
        load('sc10x5clnopreprocessing.mat');load('truelabelssc10x5clnopreprocessing.mat');Data=sc10x5clnopreprocessing'; Labels= truelabelssc10x5clnopreprocessing;
        Data = diag(1./sum(Data,2))*Data*10000;Data = log(Data+1); DataOpts.VarCutoff=1;
        
    elseif strcmp(Normalization,'Linnorm')
        
        load('sc10x5clLinnorm.mat');load('truelabelssc10x5clnopreprocessing.mat');Data=sc10x5clLinnorm'; Labels= truelabelssc10x5clnopreprocessing;
        
    elseif strcmp(Normalization,'Seurat')
        
        load('Cellmixseuratscaled.mat');load('truelabelssc10x5clnopreprocessing.mat');Data=Cellmixseuratscaled'; Labels= truelabelssc10x5clnopreprocessing;
    
    elseif strcmp(Normalization,'SCT')
        
        load('CellmixSCT.mat');load('truelabelssc10x5clnopreprocessing.mat');Data=CellmixSCT'; Labels= truelabelssc10x5clnopreprocessing;
        
    end
    
end

if strcmp(DataSet,'CellMixSng')
    
    if strcmp(Normalization,'Basic')
        
        load('Cellmixsng.mat');load('cellmixsnglabels.mat');Data=cellmixsng'; Labels= cellmixsnglabels;
        Data = diag(1./sum(Data,2))*Data*10000;Data = log(Data+1); DataOpts.VarCutoff=1;
        
    elseif strcmp(Normalization,'Linnorm')
        
        load('cellmixsngLinnorm.mat');load('cellmixsnglabels.mat');Data=cellmixsngLinnorm'; Labels=cellmixsnglabels;
        
    elseif strcmp(Normalization,'Seurat')
        
        load('cellmixsngseuratscaled.mat');load('cellmixsnglabels.mat');Data=cellmixsngseuratscaled'; Labels=cellmixsnglabels;
    
    elseif strcmp(Normalization,'SCT')
        
        load('cellmixsngSCT.mat');load('cellmixsnglabels.mat');Data=CellmixsngSCT'; Labels=cellmixsnglabels;
        
    end
    
end

if strcmp(DataSet,'Beta')
    
    if strcmp(Normalization,'Basic')
        
        load('beta_saver_no_trans_34_Gamma_10_100.mat'); Data = betasavernotrans34Gamma10100; Labels = betacellgroups34;
        Data = diag(1./sum(Data,2))*Data*10000;Data = log(Data+1); DataOpts.VarCutoff = 0.4;
        
    elseif strcmp(Normalization,'Linnorm')
        
        load('betasaver34Gamma10100Linnorm.mat');load('betacellgroups.mat');Data=betasaver34Gamma10100Linnorm';Labels=betacellgroups;
        
    elseif strcmp(Normalization,'Seurat')
        
        load('Betasaver34Gamma10100seuratscaled.mat'); Data = Betasaver34Gamma10100seuratscaled';load('betacellgroups.mat'); Labels = betacellgroups;
     elseif strcmp(Normalization,'SCT')
        
        load('betasaver34Gamma10100SCT.mat'); Data = betasaver34Gamma10100SCT';load('betacellgroups.mat'); Labels = betacellgroups;
     
    end
    
end

if strcmp(DataSet,'Beta2')
    
    if strcmp(Normalization,'Basic')
        
        load('beta3410filteredsaver.mat'); load('betacellgroups3410afterfiltering.mat'); Data = beta3410filteredsaver'; Labels = betacellgroups3410afterfiltering;
        Data = diag(1./sum(Data,2))*Data*10000;Data = log(Data+1); DataOpts.VarCutoff = 0.4;
        
    elseif strcmp(Normalization,'Linnorm')
        
        load('beta3410filteredsaverLinnorm.mat');load('betacellgroups3410afterfiltering.mat');Data=beta3410filteredsaverLinnorm'; Labels=betacellgroups3410afterfiltering;
        
    elseif strcmp(Normalization,'Seurat')
        
        load('Betafilteredsaver34Gamma10100seuratscaled.mat'); load('betacellgroups3410afterfiltering.mat'); Data = Betafilteredsaver34Gamma10100seuratscaled';load('betacellgroups.mat'); Labels=betacellgroups3410afterfiltering;
    elseif strcmp(Normalization,'SCT')
        
        load('BetafilteredSCT.mat'); load('betacellgroups3410afterfiltering.mat'); Data =BetafilteredSCT';Labels=betacellgroups3410afterfiltering;
        
    end
    
end

if strcmp(DataSet,'BaronPanc')
    
    MinClusterSampleSize = 60;
    
    if strcmp(Normalization,'Basic')
        
        load('pancreaticsavernonormalization.mat');load('pancreaticnumlabels.mat');Data=pancreaticsavernonormalization';Labels=pancreaticnumlabels;
        Data = diag(1./sum(Data,2))*Data*10000; Data =log(Data+1); DataOpts.VarCutoff = 1;
        
    elseif strcmp(Normalization,'Linnorm')
        
        load('pancreaticsaverLinnorm.mat'); load('pancreaticnumlabels.mat'); Data=pancreaticsaverLinnorm'; Labels=pancreaticnumlabels;      
        
    elseif strcmp(Normalization,'Seurat')
        
        load('BaronsPancreaticseuratscaled.mat'); load('pancreaticnumlabels.mat'); Data=BaronsPancreaticseuratscaled'; Labels=pancreaticnumlabels;
     elseif strcmp(Normalization,'SCT')
        
        load('BaronPancreaticSCT.mat'); load('pancreaticnumlabels.mat'); Data=BaronPancreaticSCT'; Labels=pancreaticnumlabels;
        
    end
    
end

if strcmp(DataSet,'TMPanc')
    
    if strcmp(Normalization,'Basic')
        
        load('pancreatictabulamurislabels.mat');load('pancreastabulamurissaver.mat'); Data=pancreastabulamurissaver';Labels=pancreatictabulamurislabels;
        Data = diag(1./sum(Data,2))*Data*10000;Data = log(Data+1); DataOpts.VarCutoff = 1;
        
    elseif strcmp(Normalization,'Linnorm')
        
        load('pancreatictabulamurislabels.mat');load('pancreastabulamurissaverLinnorm.mat'); Data=pancreastabulamurissaverLinnorm';Labels=pancreatictabulamurislabels;
        
    elseif strcmp(Normalization,'Seurat')
        
        load('pancreatictabulamurislabels.mat');load('pancreastabulamurissaverseuratscaled.mat'); Data=pancreastabulamurissaverseuratscaled';Labels=pancreatictabulamurislabels;
    elseif strcmp(Normalization,'SCT')
        
        load('pancreatictabulamurislabels.mat');load('pancreasTMSCT.mat'); Data=pancreasTMSCT';Labels=pancreatictabulamurislabels;
        
    end
    
end

if strcmp(DataSet,'TMLung')
    
    if strcmp(Normalization,'Basic')
        
        load('lungtabulamurissaver.mat');load('lungtabulamurislabels.mat'); Data=lungtabulamurissaver';Labels=lungtabulamurislabels;
        Data = diag(1./sum(Data,2))*Data*10000;Data = log(Data+1); DataOpts.VarCutoff = 2;
        
    elseif strcmp(Normalization,'Linnorm')
        
        load('lungtabulamurissaverLinnorm.mat');load('lungtabulamurislabels.mat'); Data=lungtabulamurissaverLinnorm';Labels=lungtabulamurislabels;
        
    elseif strcmp(Normalization,'Seurat')
        
        load('lungtabulamurissaverseuratscaled.mat');load('lungtabulamurislabels.mat'); Data=lungtabulamurissaverseuratscaled';Labels=lungtabulamurislabels;
        
    elseif strcmp(Normalization,'SCT')
        
        load('lungTMSCTnosaver.mat');load('lungtabulamurislabels.mat'); Data=lungTMSCTnosaver';Labels=lungtabulamurislabels;
        
    end
    
end

if strcmp(DataSet,'PBMC3k')
    
  MinClusterSampleSize = 20;   
    
    if strcmp(Normalization,'Basic')
        
        load('pbmc3k_qcsaver.mat');load('pbmc3k_labels.mat'); Data=pbmc3kqcsaver';Labels=pbmc3klabels;
        Data = diag(1./sum(Data,2))*Data*10000;Data = log(Data+1); DataOpts.VarCutoff = 2;
        
    elseif strcmp(Normalization,'Linnorm')
        
        load('pbmc3k_Linnorm.mat');load('pbmc3k_labels.mat'); Data=pbmc3kLinnorm';Labels=pbmc3klabels;
        
    elseif strcmp(Normalization,'Seurat')
        
        load('pbmc3k_saverseuratscaled.mat');load('pbmc3k_labels.mat'); Data=pbmc3ksaverseuratscaled';Labels=pbmc3klabels;
        
    elseif strcmp(Normalization,'SCT')
        
        load('pbmc3k_SCT.mat');load('pbmc3k_labels.mat'); Data=pbmc3kSCT';Labels=pbmc3klabels;
        
    end
    
end

if strcmp(DataSet,'PBMC3k_SingleR')
    
  MinClusterSampleSize = 50; %Clusters sizes: [1490 665 334 166 31]; remove smallest
    
    if strcmp(Normalization,'Basic')
        
        load('pbmc3k_qcsaver.mat');load('pbmcsingleRnumcelltypes.mat'); Data=pbmc3kqcsaver';Labels=pbmcsingleRnumcelltypes;
        Data = diag(1./sum(Data,2))*Data*10000;Data = log(Data+1); DataOpts.VarCutoff = 2;
        
    elseif strcmp(Normalization,'Linnorm')
        
        load('pbmc3k_Linnorm.mat');load('pbmcsingleRnumcelltypes.mat'); Data=pbmc3kLinnorm';Labels=pbmcsingleRnumcelltypes;
        
    elseif strcmp(Normalization,'Seurat')
        
        load('pbmc3k_saverseuratscaled.mat');load('pbmcsingleRnumcelltypes.mat'); Data=pbmc3ksaverseuratscaled';Labels=pbmcsingleRnumcelltypes;
        
    elseif strcmp(Normalization,'SCT')
        
        load('pbmc3k_SCT.mat');load('pbmcsingleRnumcelltypes.mat'); Data=pbmc3kSCT';Labels=pbmcsingleRnumcelltypes;
        
    end
    
end

if strcmp(DataSet,'PBMC4k_main')
    
  MinClusterSampleSize = 20;   
    
    if strcmp(Normalization,'Basic')
        
        load('pbmc4k_qcsaver.mat');load('pbmc4k_mainlabels.mat'); Data=pbmc4k_qcsaver';Labels=pbmc4k_mainlabels;
        Data = diag(1./sum(Data,2))*Data*10000;Data = log(Data+1); DataOpts.VarCutoff = 2;
        
    elseif strcmp(Normalization,'Linnorm')
        
        load('pbmc4k_Linnorm.mat');load('pbmc4k_mainlabels.mat'); Data=pbmc4k_Linnorm';Labels=pbmc4k_mainlabels;
        
    elseif strcmp(Normalization,'Seurat')
        
        load('pbmc4k_saverseuratscaled.mat');load('pbmc4k_mainlabels.mat'); Data=pbmc4k_saverseuratscaled';Labels=pbmc4k_mainlabels;
        
    elseif strcmp(Normalization,'SCT')
        
        load('pbmc4k_SCT.mat');load('pbmc4k_mainlabels.mat'); Data=pbmc4k_SCT';Labels=pbmc4k_mainlabels;
        
    end
    
end

if strcmp(DataSet,'PBMC4k_extented')
    
  MinClusterSampleSize = 20;   
    
    if strcmp(Normalization,'Basic')
        
        load('pbmc4k_qcsaver.mat');load('pbmc4k_extralabels.mat'); Data=pbmc4k_qcsaver';Labels=pbmc4k_extralabels;
        Data = diag(1./sum(Data,2))*Data*10000;Data = log(Data+1); DataOpts.VarCutoff = 2;
        
    elseif strcmp(Normalization,'Linnorm')
        
        load('pbmc4k_Linnorm.mat');load('pbmc4k_extralabels.mat'); Data=pbmc4k_Linnorm';Labels=pbmc4k_extralabels;
        
    elseif strcmp(Normalization,'Seurat')
        
        load('pbmc4k_saverseuratscaled.mat');load('pbmc4k_extralabels.mat'); Data=pbmc4k_saverseuratscaled';Labels=pbmc4k_extralabels;
        
    elseif strcmp(Normalization,'SCT')
        
        load('pbmc4k_SCT.mat');load('pbmc4k_extralabels.mat'); Data=pbmc4k_SCT';Labels=pbmc4k_extralabels;
        
    end
    
end
