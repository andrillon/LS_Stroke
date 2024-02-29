%% Paths, Loading & Variables
 clc, clear all
 if isempty(findstr(pwd,'Daniel'))==0
     addpath(genpath('/fs04/so34/Daniel/Script/Dependencies/'));
     rmpath(genpath('/fs04/so34/Daniel/Functions/eeglab2021.0/'));
     addpath(genpath('/fs04/so34/Daniel/Functions/LSCPtools/'));
     path_temp=('/fs04/so34/Daniel/Data/SWdetection/');
     path_raw=('/fs04/so34/Daniel/Data/Raw/');
     
     load('/fs04/so34/Daniel/Script/Dependencies/chanlocs_Oxford.mat');
     load('/fs04/so34/Daniel/Script/Dependencies/chanlocs_Monash.mat');
     load([path_temp 'SW_individualThreshold_categorical_subGroups.mat']);
     SW_table=readtable([path_temp 'SW_individualThreshold_acrossBlocks.csv']);
     A=load([path_temp 'SW_individualThreshold_categorical_subGroups.mat']);
     commonChans=A.commonChans;
     demo_table=readtable([path_raw 'Lycette_StrokeMeta_260423.xlsx']);
 else
%           addpath(genpath('/fs04/so34/Daniel/Script/Dependencies/'));
     addpath(('/Users/thandrillon/WorkGit/projects/ext/fieldtrip'));
     ft_defaults;
     addpath(genpath('/Users/thandrillon/WorkGit/LSCPtools/'));
     path_temp=('/Users/thandrillon/WorkGit/projects/inprogress/LS_Stroke/');
     
     load('/Users/thandrillon/WorkGit/projects/inprogress/LS_Stroke/chanlocs_Oxford.mat');
     load('/Users/thandrillon/WorkGit/projects/inprogress/LS_Stroke/chanlocs_Monash.mat');
     SW_table=readtable([path_temp 'SW_individualThreshold.csv']);
 end

ROI_SWDens=[]; all_SWDens=[];
GroupIDs=unique(SW_table.subGroupID);
stroke_SubIDs=unique(SW_table.SubID(ismember(SW_table.subGroupID,GroupIDs(2)) | ismember(SW_table.subGroupID,GroupIDs(3))));
lstr_SubIDs=unique(SW_table.SubID(ismember(SW_table.subGroupID,GroupIDs(2))));
rstr_SubIDs=unique(SW_table.SubID(ismember(SW_table.subGroupID,GroupIDs(3))));
neglect_current = categorical({'S010','S017','S027','S032','S018','S201','S202'});
NegGroup = categorical({'control','no_neglect','neglect','left_hemi'});
allsubj=unique(SW_table.SubID);

SW_table.SubID=categorical(SW_table.SubID);
SW_table.subGroupID=categorical(SW_table.subGroupID);
SW_table.Elec=categorical(SW_table.Elec);

% Define neglect / no neglect groups for later
for pps = 1:size(allsubj)
        if ismember(allsubj(pps),lstr_SubIDs)
            SW_table.NegGroup(ismember(SW_table.SubID,allsubj(pps))) = NegGroup(4);
        elseif ismember(allsubj(pps),rstr_SubIDs)
            if ismember(allsubj(pps),neglect_current)
                SW_table.NegGroup(ismember(SW_table.SubID,allsubj(pps))) = NegGroup(3);
            else
                SW_table.NegGroup(ismember(SW_table.SubID,allsubj(pps))) = NegGroup(2);
            end
        else
            SW_table.NegGroup(ismember(SW_table.SubID,allsubj(pps))) = NegGroup(1);
        end
end

SW_table.NegGroup = categorical(SW_table.NegGroup);

% Loading Demographics & Stroke Variables
allsubj=unique(SW_table.SubID);
for pps = 1:length(allsubj)
    SW_table.Age(ismember(SW_table.SubID,allsubj(pps)))=demo_table.Age(ismember(demo_table.ID,allsubj(pps)));
    SW_table.LesionVol(ismember(SW_table.SubID,allsubj(pps)))=demo_table.StrokeVolume(ismember(demo_table.ID,allsubj(pps)));
    SW_table.TSS(ismember(SW_table.SubID,allsubj(pps)))=demo_table.StrokeTimeYrs(ismember(demo_table.ID,allsubj(pps)));
    SW_table.IQ(ismember(SW_table.SubID,allsubj(pps)))=demo_table.PremIQ(ismember(demo_table.ID,allsubj(pps)));
    SW_table.CRI(ismember(SW_table.SubID,allsubj(pps)))=demo_table.CRI_Tot(ismember(demo_table.ID,allsubj(pps)));
    SW_table.PCRS(ismember(SW_table.SubID,allsubj(pps)))=demo_table.PCRS(ismember(demo_table.ID,allsubj(pps)));
    SW_table.FSS(ismember(SW_table.SubID,allsubj(pps)))=demo_table.FSS(ismember(demo_table.ID,allsubj(pps)));
    SW_table.CFQ(ismember(SW_table.SubID,allsubj(pps)))=demo_table.CFQ(ismember(demo_table.ID,allsubj(pps)));
    SW_table.MoCA(ismember(SW_table.SubID,allsubj(pps)))=demo_table.MoCA(ismember(demo_table.ID,allsubj(pps)));
end

% Load Cluster Output


% Channel Layout
cfg = [];
cfg.layout = 'biosemi64.lay';
layout = ft_prepare_layout(cfg);

cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel = layout.label;
a = match_str(layout.label,commonChans);
cfg.channel = layout.label(a);
cfg.center = 'yes';
layout = ft_prepare_layout(cfg);

%% ROI Definition

% Calculate mean SW density per block in the ROI for the stroke groups
for pps = 1:size(stroke_SubIDs)
    for nB = 1:9
        if ismember(stroke_SubIDs(pps),lstr_SubIDs)
            SW_table.ROI_SWDens(ismember(SW_table.SubID,stroke_SubIDs(pps)) & SW_table.Block==nB) = nanmean(SW_table.SW_density(ismember(SW_table.SubID,stroke_SubIDs(pps)) & ismember(SW_table.Elec,ClustersByDrugs{2,1}) & SW_table.Block==nB));
        elseif ismember(stroke_SubIDs(pps),rstr_SubIDs)
            SW_table.ROI_SWDens(ismember(SW_table.SubID,stroke_SubIDs(pps)) & SW_table.Block==nB) = nanmean(SW_table.SW_density(ismember(SW_table.SubID,stroke_SubIDs(pps)) & ismember(SW_table.Elec,ClustersByDrugs{2,2}) & SW_table.Block==nB));
        end
    end
    ROI_mean_SWDens(pps) = nanmean(SW_table.ROI_SWDens(ismember(SW_table.SubID,stroke_SubIDs(pps)))); % Overall mean per participant
    TSS(pps) = demo_table.StrokeTimeYrs(ismember(demo_table.ID,stroke_SubIDs(pps))); % Time since stroke
    LesionVol(pps) = demo_table.StrokeVolume(ismember(demo_table.ID,stroke_SubIDs(pps))); % Lesion volume
end

% Calculating BRDM stats - RT, hit rate per block per participant
allsubj=unique(SW_table.SubID); RT_all = []; sides = {'Left','Right'}; nCh=1;
SW_perTrial = readtable([path_temp 'SW_and_Behav_perTrial.csv']);
SW_perTrial.SubID=categorical(SW_perTrial.SubID);
SW_perTrial.GroupID=categorical(SW_perTrial.GroupID);
SW_perTrial.Elec=categorical(SW_perTrial.Elec);
SW_perTrial.StimSide=categorical(SW_perTrial.StimSide);

for pps = 1:size(allsubj)
    for nB = 1:9
        % RT by block - all, left & right targets
        SW_table.RT_block(ismember(SW_table.SubID,allsubj(pps)) & SW_table.Block==nB) = nanmean(SW_perTrial.RT(ismember(SW_perTrial.SubID,allsubj(pps)) & ismember(SW_perTrial.Elec,layout.label{nCh}) & SW_perTrial.Block==nB & SW_perTrial.RT>=150 & SW_perTrial.RT<=3000));
        temp_left = nanmean(SW_perTrial.RT(ismember(SW_perTrial.SubID,allsubj(pps)) & ismember(SW_perTrial.Elec,layout.label{nCh}) & SW_perTrial.Block==nB & SW_perTrial.RT>=150 & SW_perTrial.RT<=3000 & ismember(SW_perTrial.StimSide,sides{1})));
        SW_table.RTleft_block(ismember(SW_table.SubID,allsubj(pps)) & SW_table.Block==nB) = temp_left;
        temp_right = nanmean(SW_perTrial.RT(ismember(SW_perTrial.SubID,allsubj(pps)) & ismember(SW_perTrial.Elec,layout.label{nCh}) & SW_perTrial.Block==nB & SW_perTrial.RT>=150 & SW_perTrial.RT<=3000 & ismember(SW_perTrial.StimSide,sides{2})));;
        SW_table.RTright_block(ismember(SW_table.SubID,allsubj(pps)) & SW_table.Block==nB) = temp_right;
        SW_table.RTasym_block(ismember(SW_table.SubID,allsubj(pps)) & SW_table.Block==nB) = (temp_left-temp_right)/((temp_left+temp_right)/2);
        % hits by block - all, left & right targets
        temp_hits = sum(SW_perTrial.Response(ismember(SW_perTrial.SubID,allsubj(pps)) & ismember(SW_perTrial.Elec,layout.label{nCh}) & SW_perTrial.Block==nB & SW_perTrial.RT>=150 & SW_perTrial.RT<=3000 & SW_perTrial.FixBreak~=1) == 1);
        SW_table.hits_block(ismember(SW_table.SubID,allsubj(pps)) & SW_table.Block==nB) = temp_hits;
        temp_hits_left = sum(SW_perTrial.Response(ismember(SW_perTrial.SubID,allsubj(pps)) & ismember(SW_perTrial.Elec,layout.label{nCh}) & SW_perTrial.Block==nB & SW_perTrial.RT>=150 & SW_perTrial.RT<=3000 & SW_perTrial.FixBreak~=1 & ismember(SW_perTrial.StimSide,sides{1})) == 1);
        SW_table.hits_left_block(ismember(SW_table.SubID,allsubj(pps)) & SW_table.Block==nB) = temp_hits_left;
        temp_hits_right = sum(SW_perTrial.Response(ismember(SW_perTrial.SubID,allsubj(pps)) & ismember(SW_perTrial.Elec,layout.label{nCh}) & SW_perTrial.Block==nB & SW_perTrial.RT>=150 & SW_perTrial.RT<=3000 & SW_perTrial.FixBreak~=1 & ismember(SW_perTrial.StimSide,sides{2})) == 1);
        SW_table.hits_right_block(ismember(SW_table.SubID,allsubj(pps)) & SW_table.Block==nB) = temp_hits_right;
        temp_validt = temp_hits+sum(SW_perTrial.Response(ismember(SW_perTrial.SubID,allsubj(pps)) & ismember(SW_perTrial.Elec,layout.label{nCh}) & SW_perTrial.Block==nB & SW_perTrial.FixBreak~=1) == 0);
        SW_table.validt_block(ismember(SW_table.SubID,allsubj(pps)) & SW_table.Block==nB) = temp_validt;
        temp_validt_left = temp_hits_left+sum(SW_perTrial.Response(ismember(SW_perTrial.SubID,allsubj(pps)) & ismember(SW_perTrial.Elec,layout.label{nCh}) & SW_perTrial.Block==nB & SW_perTrial.FixBreak~=1& ismember(SW_perTrial.StimSide,sides{1})) == 0);
        SW_table.validt_block_left(ismember(SW_table.SubID,allsubj(pps)) & SW_table.Block==nB) = temp_validt_left;
        temp_validt_right = temp_hits_right+sum(SW_perTrial.Response(ismember(SW_perTrial.SubID,allsubj(pps)) & ismember(SW_perTrial.Elec,layout.label{nCh}) & SW_perTrial.Block==nB & SW_perTrial.FixBreak~=1& ismember(SW_perTrial.StimSide,sides{2})) == 0);
        SW_table.validt_block_right(ismember(SW_table.SubID,allsubj(pps)) & SW_table.Block==nB) = temp_validt_right;
        % Hit rate by block - all, left & right targets
        SW_table.hitrate_block(ismember(SW_table.SubID,allsubj(pps)) & SW_table.Block==nB) = (temp_hits/temp_validt)*100;
        temp_hitrate_left = (temp_hits_left/temp_validt_left)*100;
        SW_table.hitrate_block_left(ismember(SW_table.SubID,allsubj(pps)) & SW_table.Block==nB) = temp_hitrate_left;
        temp_hitrate_right = (temp_hits_right/temp_validt_right)*100;
        SW_table.hitrate_block_right(ismember(SW_table.SubID,allsubj(pps)) & SW_table.Block==nB) = temp_hitrate_right;
        SW_table.hitrate_asym_block(ismember(SW_table.SubID,allsubj(pps)) & SW_table.Block==nB) = (temp_hitrate_left-temp_hitrate_right)/((temp_hitrate_left+temp_hitrate_right)/2);
    end
    all_SWDens(pps) = nanmean(SW_table.SW_density(ismember(SW_table.SubID,allsubj(pps)) & ismember(SW_table.Elec,layout.label{nCh}) & SW_table.Block<=9));
    RT_all(pps) = nanmean(SW_perTrial.RT(ismember(SW_perTrial.SubID,allsubj(pps)) & ismember(SW_perTrial.Elec,layout.label{nCh}) & SW_perTrial.Block<=9 & SW_perTrial.RT>=150 & SW_perTrial.RT<=3000));
    groups_all(pps) = unique(SW_table.subGroupID(ismember(SW_table.SubID,allsubj(pps))));
    SW_table.RT_all(ismember(SW_table.SubID,allsubj(pps))) = RT_all(pps);
    temp_left = nanmean(SW_perTrial.RT(ismember(SW_perTrial.SubID,allsubj(pps)) & ismember(SW_perTrial.Elec,layout.label{nCh}) & SW_perTrial.Block<=9 & SW_perTrial.RT>=150 & SW_perTrial.RT<=3000 & ismember(SW_perTrial.StimSide,sides{1})));
    SW_table.RTleft_all(ismember(SW_table.SubID,allsubj(pps))) = temp_left;
    temp_right = nanmean(SW_perTrial.RT(ismember(SW_perTrial.SubID,allsubj(pps)) & ismember(SW_perTrial.Elec,layout.label{nCh}) & SW_perTrial.Block<=9 & SW_perTrial.RT>=150 & SW_perTrial.RT<=3000 & ismember(SW_perTrial.StimSide,sides{2})));
    SW_table.RTright_all(ismember(SW_table.SubID,allsubj(pps))) = temp_right;
    RTasym_all(pps) = (temp_left-temp_right)/((temp_left+temp_right)/2);
    SW_table.RTasym_all(ismember(SW_table.SubID,allsubj(pps))) = RTasym_all(pps);
    temp_hits = sum(SW_perTrial.Response(ismember(SW_perTrial.SubID,allsubj(pps)) & ismember(SW_perTrial.Elec,layout.label{nCh}) & SW_perTrial.RT>=150 & SW_perTrial.RT<=3000 & SW_perTrial.FixBreak~=1) == 1);
    SW_table.hits_all(ismember(SW_table.SubID,allsubj(pps))) = temp_hits;
    temp_hits_left = sum(SW_perTrial.Response(ismember(SW_perTrial.SubID,allsubj(pps)) & ismember(SW_perTrial.Elec,layout.label{nCh}) & SW_perTrial.RT>=150 & SW_perTrial.RT<=3000 & SW_perTrial.FixBreak~=1 & ismember(SW_perTrial.StimSide,sides{1})) == 1);
    SW_table.hits_all_left(ismember(SW_table.SubID,allsubj(pps))) = temp_hits_left;
    temp_hits_right = sum(SW_perTrial.Response(ismember(SW_perTrial.SubID,allsubj(pps)) & ismember(SW_perTrial.Elec,layout.label{nCh}) & SW_perTrial.RT>=150 & SW_perTrial.RT<=3000 & SW_perTrial.FixBreak~=1 & ismember(SW_perTrial.StimSide,sides{2})) == 1);
    SW_table.hits_all_right(ismember(SW_table.SubID,allsubj(pps))) = temp_hits_right;
    temp_validt = temp_hits+sum(SW_perTrial.Response(ismember(SW_perTrial.SubID,allsubj(pps)) & ismember(SW_perTrial.Elec,layout.label{nCh}) & SW_perTrial.FixBreak~=1) == 0);
    SW_table.validt_all(ismember(SW_table.SubID,allsubj(pps))) = temp_validt;
    temp_validt_left = temp_hits_left+sum(SW_perTrial.Response(ismember(SW_perTrial.SubID,allsubj(pps)) & ismember(SW_perTrial.Elec,layout.label{nCh}) & SW_perTrial.FixBreak~=1 & ismember(SW_perTrial.StimSide,sides{1})) == 0);
    SW_table.validt_left(ismember(SW_table.SubID,allsubj(pps))) = temp_validt_left;
    temp_validt_right = temp_hits_right+sum(SW_perTrial.Response(ismember(SW_perTrial.SubID,allsubj(pps)) & ismember(SW_perTrial.Elec,layout.label{nCh}) & SW_perTrial.FixBreak~=1 & ismember(SW_perTrial.StimSide,sides{2})) == 0);
    SW_table.validt_right(ismember(SW_table.SubID,allsubj(pps))) = temp_validt_right;
    hitrate_all(pps) = (temp_hits/temp_validt)*100;
    SW_table.hitrate_all(ismember(SW_table.SubID,allsubj(pps))) = hitrate_all(pps);
    temp_hitrate_left = (temp_hits_left/temp_validt_left)*100;
    SW_table.hitrate_all_left(ismember(SW_table.SubID,allsubj(pps))) = temp_hitrate_left;
    temp_hitrate_right = (temp_hits_right/temp_validt_right)*100;
    SW_table.hitrate_all_right(ismember(SW_table.SubID,allsubj(pps))) = temp_hitrate_right;
    hitrateasym_all(pps) = (temp_hitrate_left-temp_hitrate_right)/((temp_hitrate_left+temp_hitrate_right)/2);
    SW_table.hitrate_all_asym(ismember(SW_table.SubID,allsubj(pps))) = hitrateasym_all(pps);
end

%% Density x BRDM Behaviour LMMs

% Hit rate x ROI SW Density
sub_table=SW_table(SW_table.Elec==layout.label{1} & SW_table.Block <=9 & (ismember(SW_table.subGroupID,GroupIDs(2)) | ismember(SW_table.subGroupID,GroupIDs(3))),:);
mdl=fitlme(sub_table,'hitrate_block~ROI_SWDens+Block+(1|SubID)')
figure; scatter(sub_table.hitrate_block,sub_table.ROI_SWDens,'filled'); lsline;

% RT x ROI SW Density
sub_table=SW_table(SW_table.Elec==layout.label{1} & SW_table.Block <=9 & (ismember(SW_table.subGroupID,GroupIDs(2)) | ismember(SW_table.subGroupID,GroupIDs(3))),:);
mdl=fitlme(sub_table,'RT_block~ROI_SWDens+Block+(1|SubID)')
figure; scatter(sub_table.RT_block,sub_table.ROI_SWDens,'filled'); lsline;

% Hit rate Asymmetry x ROI SW Density
sub_table=SW_table(SW_table.Elec==layout.label{1} & SW_table.Block <=9 & (ismember(SW_table.subGroupID,GroupIDs(2)) | ismember(SW_table.subGroupID,GroupIDs(3))),:);
mdl=fitlme(sub_table,'hitrate_asym_block~ROI_SWDens+Block+(1|SubID)')
figure; scatter(sub_table.hitrate_asym_block,sub_table.ROI_SWDens,'filled'); lsline;

% RT Asymmetry x ROI SW Density
sub_table=SW_table(SW_table.Elec==layout.label{1} & SW_table.Block <=9 & (ismember(SW_table.subGroupID,GroupIDs(2)) | ismember(SW_table.subGroupID,GroupIDs(3))),:);
mdl=fitlme(sub_table,'RTasym_block~ROI_SWDens+Block+(1|SubID)')
figure; scatter(sub_table.RTasym_block,sub_table.ROI_SWDens,'filled'); lsline;

% Scalp SW Density x RT - Partial Correlation
[rho, pval]=partialcorr(RT_all',all_SWDens',grp2idx(groups_all))
figure; scatter(RT_all,all_SWDens,'filled'); lsline;

% Scalp SW Density x RT Asymmetry - Partial Correlation
[rho, pval]=partialcorr(RTasym_all',all_SWDens',grp2idx(groups_all))
figure; scatter(RTasym_all,all_SWDens,'filled'); lsline;

% Scalp SW Density x Hit Rate - Partial Correlation
[rho, pval]=partialcorr(hitrate_all',all_SWDens',grp2idx(groups_all))
figure; scatter(hitrate_all,all_SWDens,'filled'); lsline;

% Scalp SW Density x Hit Rate Asymmetry - Partial Correlation
[rho, pval]=partialcorr(hitrateasym_all',all_SWDens',grp2idx(groups_all))
figure; scatter(hitrateasym_all,all_SWDens,'filled'); lsline;

%% Stroke Characteristic Analyses

% Lesion Volume x ROI SW Density LMM
sub_table=SW_table(SW_table.Elec==layout.label{1} & SW_table.Block <=9 & (ismember(SW_table.subGroupID,GroupIDs(2)) | ismember(SW_table.subGroupID,GroupIDs(3))),:);
mdl=fitlme(sub_table,'ROI_SWDens~LesionVol*Block+(1|SubID)')
figure; scatter(sub_table.ROI_SWDens,sub_table.LesionVol,'filled'); lsline;

% Time Since Stroke x ROI SW Density LMM
sub_table=SW_table(SW_table.Elec==layout.label{1} & SW_table.Block <=9 & (ismember(SW_table.subGroupID,GroupIDs(2)) | ismember(SW_table.subGroupID,GroupIDs(3))),:);
mdl=fitlme(sub_table,'ROI_SWDens~TSS*Block+(1|SubID)')
figure; scatter(sub_table.ROI_SWDens,sub_table.TSS,'filled'); lsline;

% TSS*Block Interaction Post-Hocs
for n=1:9
    sub_table=SW_table(SW_table.Elec==layout.label{1} & SW_table.Block == n & (ismember(SW_table.subGroupID,GroupIDs(2)) | ismember(SW_table.subGroupID,GroupIDs(3))),:);
    mdl=fitlme(sub_table,'ROI_SWDens~TSS+(1|SubID)')
    stats = anova(mdl)
    pvals(n) = stats.pValue(2)
end
[corrected_p, ~] = bonf_holm(pvals,0.05)

% Neglect vs non-neglect sub-analysis
sub_table=SW_table(SW_table.Elec==layout.label{1} & SW_table.Block <=9 & ismember(SW_table.NegGroup,1),:);
mdl=fitlme(sub_table,'RTasym_block~ROI_SWDens+Block+(1|SubID)')
figure; scatter(sub_table.RTasym_block,sub_table.ROI_SWDens,'filled'); lsline;

%% Functional Questionnaire x SW Density Partial Correlations
% ROI SW Density x PCRS - Partial Correlation
[rho, pval]=partialcorr(PCRS',all_SWDens',grp2idx(groups_all))
figure; scatter(PCRS,all_SWDens,'filled'); lsline;

% ROI SW Density x FSS - Partial Correlation
[rho, pval]=partialcorr(FSS',all_SWDens',grp2idx(groups_all))
figure; scatter(FSS,all_SWDens,'filled'); lsline;

% ROI SW Density x CFQ - Partial Correlation
[rho, pval]=partialcorr(CFQ',all_SWDens',grp2idx(groups_all))
figure; scatter(CFQ,all_SWDens,'filled'); lsline;

    
%% Exploratory Analysis - Focal Strokes
Groups = unique(SW_table.subGroupID); data_to_plot=[];
focal = categorical({'S109','S114','S004','S025','S030','S033','S201','S010','S017'});


for nG=1:3
    for nCh=1:size(layout.label)-2
        data_to_plot(nG,nCh)=nanmean(SW_table.SW_density(SW_table.Elec==layout.label{nCh}&SW_table.Block<=9&ismember(SW_table.subGroupID,Groups(nG))&(ismember(SW_table.SubID,focal) | ismember(SW_table.subGroupID,Groups(1))))); % Get nanmean of variable for rows of group, blocks 1:9, nCh
    end
end

% Separate loop for plot so that I can define c axis similarly across groups

for nG=1:3
    figure;
    simpleTopoPlot_ft(data_to_plot(nG,:), layout,'on',[],0,1); % Plot a topo for each group
    colorbar;
    caxis([min(min(data_to_plot)) max(max(data_to_plot))])
    title([Groups(nG) ': SW Density (Focal Strokes)']);
    
end