%% Two groups
%{
clc, clear all

addpath(genpath('/fs04/so34/Daniel/Script/Dependencies/'));
addpath(genpath('/fs04/so34/Daniel/Functions/LSCPtools/'));
path_temp=('/fs04/so34/Daniel/Data/SWdetection/');

load('/fs04/so34/Daniel/Script/Dependencies/chanlocs_Oxford.mat');
load('/fs04/so34/Daniel/Script/Dependencies/chanlocs_Monash.mat');
load([path_temp 'SW_individualThreshold_withICA4.mat']);

% addpath(genpath('/Users/manab/Desktop/Functions/eeglab2021.0/'));
% eeglab;
% addpath(genpath('/Users/manab/Desktop/Functions/fieldtrip/'));
% 
% load('SW_table_twoGroups.mat')
% load('commonChans.mat')
% SW_table = SW_table_twoGroups;

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

Group_effect=[];
for nCh=1:length(layout.label)-2
    sub_table=SW_table(SW_table.Elec==layout.label{nCh},:);
    mdl=fitlme(sub_table,'SW_density~1+GroupID+(1|SubID)');
    Group_effect(nCh,1)=mdl.Coefficients.tStat(find_trials(mdl.Coefficients.Name,'GroupID'));
    Group_effect(nCh,2)=mdl.Coefficients.pValue(find_trials(mdl.Coefficients.Name,'GroupID'));
end

cmap2=cbrewer('div','RdBu',64); cmap2=flipud(cmap2);

figure;
simpleTopoPlot_ft(Group_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(Group_effect(:,2)<fdr(Group_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(Group_effect(:,1))));
title('Group Effect');
%title('Misses', 'FontSize', 16)

%%
SW_density=[];
SW_dslope=[];
SW_amplitude=[];
for nCh=1:length(layout.label)-2
    sub_table=SW_table(SW_table.Elec==layout.label{nCh},:);
    SW_density(nCh,1)=nanmean(sub_table.SW_density(sub_table.GroupID=='healthy_old'));
    SW_density(nCh,2)=nanmean(sub_table.SW_density(sub_table.GroupID=='neglect'));
    
    SW_amplitude(nCh,1)=nanmean(sub_table.SW_amplitude(sub_table.GroupID=='healthy_old'));
    SW_amplitude(nCh,2)=nanmean(sub_table.SW_amplitude(sub_table.GroupID=='neglect'));
    
    SW_dslope(nCh,1)=nanmean(sub_table.SW_downslope(sub_table.GroupID=='healthy_old'));
    SW_dslope(nCh,2)=nanmean(sub_table.SW_downslope(sub_table.GroupID=='neglect'));
end
figure;
% subplot(3,2,1);
simpleTopoPlot_ft(SW_density(:,1), layout,'on',[],0,1);
% subplot(3,2,2);
figure;
simpleTopoPlot_ft(SW_density(:,2), layout,'on',[],0,1);

subplot(3,2,3);
simpleTopoPlot_ft(SW_amplitude(:,1), layout,'on',[],0,1);
subplot(3,2,4);
simpleTopoPlot_ft(SW_amplitude(:,2), layout,'on',[],0,1);

subplot(3,2,5);
simpleTopoPlot_ft(SW_dslope(:,1), layout,'on',[],0,1);
subplot(3,2,6);
simpleTopoPlot_ft(SW_dslope(:,2), layout,'on',[],0,1);
%}
%% Three Groups
 clc, clear all
 if isempty(findstr(pwd,'Daniel'))==0
     addpath(genpath('/fs04/so34/Daniel/Script/Dependencies/'));
     addpath(genpath('/fs04/so34/Daniel/Functions/LSCPtools/'));
     path_temp=('/fs04/so34/Daniel/Data/SWdetection/');
     
     load('/fs04/so34/Daniel/Script/Dependencies/chanlocs_Oxford.mat');
     load('/fs04/so34/Daniel/Script/Dependencies/chanlocs_Monash.mat');
     load([path_temp 'SW_individualThreshold_withICA4_subGroups.mat']);
 else
%           addpath(genpath('/fs04/so34/Daniel/Script/Dependencies/'));
     addpath(genpath('/Users/thandrillon/WorkGit/LSCPtools/'));
     path_temp=('/Users/thandrillon/WorkGit/projects/inprogress/LS_Stroke/');
     
     load('/Users/thandrillon/WorkGit/projects/inprogress/LS_Stroke/chanlocs_Oxford.mat');
     load('/Users/thandrillon/WorkGit/projects/inprogress/LS_Stroke/chanlocs_Monash.mat');
     load([path_temp 'SW_individualThreshold_withICA4_subGroups.mat']);
 end
% Reorder for L vs R comparison
SW_table_reorder=SW_table;
SW_table_reorder.subGroupID=reordercats(SW_table_reorder.subGroupID,{'left_stroke','right_stroke','healthy_old'});

% Mark's suggested analysis
% SW_table2=SW_table; % Copy the table
% SW_table2(~ismember(SW_table2.Elec,commonChans),:)=[]; % Remove rows where the electrode is not common across sites
% 
% subGroupID2=zeros(height(SW_table2),1); % Create blank subGroupID variable
% subGroupID2(ismember(SW_table2.subGroupID,'healthy_old'))=1; % Set rows with healthy older to 1
% subGroupID2(ismember(SW_table2.subGroupID,'left_stroke'))=2; % Set rows with left stroke to 2
% subGroupID2(ismember(SW_table2.subGroupID,'right_stroke'))=3; % Set rows with right stroke to 3
% 
% Elec2=zeros(height(SW_table2),1); % Create blank electrode variable
% for elec=1:59 % For all electrodes
%     currelec=commonChans(elec); % Choose one electrode from common chans
%     Elec2(ismember(SW_table2.Elec,currelec))=elec; % Set all rows with that electrode label to the same number
% end
% 
% SubID2=zeros(height(SW_table2),1); % Create blank SubID variable
% allsubs=unique(SW_table2.SubID); % Get all SubIDs
% for sub=1:57 % For all subjects
%     currsub=allsubs(sub); % Choose one SubID from list
%     SubID2(ismember(SW_table2.SubID,currsub))=sub; % Set all rows with that participant label to the same number
% end
% 
% SW_table2.SubID2=SubID2; % Add new sub IDs to table
% SW_table2.SubID2=nominal(SW_table2.SubID2); % Convert new sub IDs to nominal
% SW_table2.subGroupID2=subGroupID2; % Add new group IDs to table
% SW_table2.subGroupID2=nominal(SW_table2.subGroupID2); % Convert new group IDs to nominal
% SW_table2.Elec2=Elec2; % Add new sub IDs to table
% SW_table2.Elec2=nominal(SW_table2.Elec2); % Convert new sub IDs to nominal

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

% Density
LeftStroke_effect=[]; RightStroke_effect=[];
% mdl=fitlme(SW_table,'SW_density~subGroupID*Elec + (1|SubID)'); % Mark
for nCh=1:length(layout.label)-2
    sub_table=SW_table(SW_table.Elec==layout.label{nCh},:);
    mdl=fitlme(sub_table,'SW_density~1+subGroupID+(1|SubID)');
    LeftStroke_effect(nCh,1)=mdl.Coefficients.tStat(find_trials(mdl.Coefficients.Name,'left'));
    LeftStroke_effect(nCh,2)=mdl.Coefficients.pValue(find_trials(mdl.Coefficients.Name,'left'));
    RightStroke_effect(nCh,1)=mdl.Coefficients.tStat(find_trials(mdl.Coefficients.Name,'right'));
    RightStroke_effect(nCh,2)=mdl.Coefficients.pValue(find_trials(mdl.Coefficients.Name,'right'));
end

% And then left vs right stroke
for nCh=1:length(layout.label)-2
    sub_table=SW_table_reorder(SW_table_reorder.Elec==layout.label{nCh},:);
    mdl=fitlme(sub_table,'SW_density~1+subGroupID+(1|SubID)');
    LxR_effect(nCh,1)=mdl.Coefficients.tStat(find_trials(mdl.Coefficients.Name,'right'));
    LxR_effect(nCh,2)=mdl.Coefficients.pValue(find_trials(mdl.Coefficients.Name,'right'));
end

cmap2=cbrewer('div','RdBu',64); cmap2=flipud(cmap2);

figure;
simpleTopoPlot_ft(LeftStroke_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(LeftStroke_effect(:,2)<0.05), 'pointsymbol','o','pointcolor','green','pointsize',36,'box','no','label','no')
ft_plot_lay_me(layout, 'chanindx', find(LeftStroke_effect(:,2)<fdr(LeftStroke_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(LeftStroke_effect(:,1))))
title('Healthy x Left Stroke ~ SW Density');

figure;
simpleTopoPlot_ft(RightStroke_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(RightStroke_effect(:,2)<0.05), 'pointsymbol','o','pointcolor','green','pointsize',36,'box','no','label','no')
ft_plot_lay_me(layout, 'chanindx', find(RightStroke_effect(:,2)<fdr(RightStroke_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(RightStroke_effect(:,1))))
title('Healthy x Right Stroke ~ SW Density');

figure;
simpleTopoPlot_ft(LxR_effect(:,1), layout,'on',[],0,1); hold on;
% ft_plot_lay_me(layout, 'chanindx', find(LxR_effect(:,2)<0.05), 'pointsymbol','o','pointcolor','green','pointsize',36,'box','no','label','no')
% ft_plot_lay_me(layout, 'chanindx', find(LxR_effect(:,2)<fdr(LxR_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(LxR_effect(:,1))))
title('Left x Right Stroke ~ SW Density');

%% Cluster Permutations
clusteralpha = 0.05;
montecarloalpha = 0.05
nperm = 1000; % Need to edit for however many is possible
tailFlag = 0;

cfg_neighb = [];
cfg_neighb.method = 'template';
cfg_neighb.layout = layout;
cfg_neighb.channel = layout.label;
neighbours = ft_prepare_neighbours(cfg_neighb);

[SW_clus] = get_clusterperm_lme_lsstroke(SW_est,clusteralpha,montecarloalpha,nperm,neighbours);
% What is the input file (SW_est)? Is it the LME outputs? I.e., run LME at each elec and put in SW_est
% Structure = Elec x ? x t value x p value x comparison/group?

%% Amplitude
%{
LeftStroke_effect=[]; RightStroke_effect=[]; LxR_effect=[];
for nCh=1:length(layout.label)-2
    sub_table=SW_table(SW_table.Elec==layout.label{nCh},:);
    mdl=fitlme(sub_table,'SW_amplitude~1+subGroupID+(1|SubID)');
    LeftStroke_effect(nCh,1)=mdl.Coefficients.tStat(find_trials(mdl.Coefficients.Name,'left'));
    LeftStroke_effect(nCh,2)=mdl.Coefficients.pValue(find_trials(mdl.Coefficients.Name,'left'));
    RightStroke_effect(nCh,1)=mdl.Coefficients.tStat(find_trials(mdl.Coefficients.Name,'right'));
    RightStroke_effect(nCh,2)=mdl.Coefficients.pValue(find_trials(mdl.Coefficients.Name,'right'));
end

% And then left vs right stroke
for nCh=1:length(layout.label)-2
    sub_table=SW_table_reorder(SW_table_reorder.Elec==layout.label{nCh},:);
    mdl=fitlme(sub_table,'SW_amplitude~1+subGroupID+(1|SubID)');
    LxR_effect(nCh,1)=mdl.Coefficients.tStat(find_trials(mdl.Coefficients.Name,'right'));
    LxR_effect(nCh,2)=mdl.Coefficients.pValue(find_trials(mdl.Coefficients.Name,'right'));
end

cmap2=cbrewer('div','RdBu',64); cmap2=flipud(cmap2);

figure;
simpleTopoPlot_ft(LeftStroke_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(LeftStroke_effect(:,2)<0.05), 'pointsymbol','o','pointcolor','green','pointsize',36,'box','no','label','no')
% ft_plot_lay_me(layout, 'chanindx', find(LeftStroke_effect(:,2)<fdr(LeftStroke_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(LeftStroke_effect(:,1))))
title('Healthy Older x Left Stroke ~ SW Amplitude');

figure;
simpleTopoPlot_ft(RightStroke_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(RightStroke_effect(:,2)<0.05), 'pointsymbol','o','pointcolor','green','pointsize',36,'box','no','label','no')
% ft_plot_lay_me(layout, 'chanindx', find(RightStroke_effect(:,2)<fdr(RightStroke_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(RightStroke_effect(:,1))))
title('Healthy Older x Right Stroke ~ SW Amplitude');

figure;
simpleTopoPlot_ft(LxR_effect(:,1), layout,'on',[],0,1); hold on;
ft_plot_lay_me(layout, 'chanindx', find(LxR_effect(:,2)<0.05), 'pointsymbol','o','pointcolor','green','pointsize',36,'box','no','label','no')
% ft_plot_lay_me(layout, 'chanindx', find(LxR_effect(:,2)<fdr(LxR_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(LxR_effect(:,1))))
title('Left x Right Stroke ~ SW Amplitude');

%% Frequency
LeftStroke_effect=[]; RightStroke_effect=[];
for nCh=1:length(layout.label)-2
    sub_table=SW_table(SW_table.Elec==layout.label{nCh},:);
    mdl=fitlme(sub_table,'SW_frequency~1+subGroupID+(1|SubID)');
    LeftStroke_effect(nCh,1)=mdl.Coefficients.tStat(find_trials(mdl.Coefficients.Name,'left'));
    LeftStroke_effect(nCh,2)=mdl.Coefficients.pValue(find_trials(mdl.Coefficients.Name,'left'));
    RightStroke_effect(nCh,1)=mdl.Coefficients.tStat(find_trials(mdl.Coefficients.Name,'right'));
    RightStroke_effect(nCh,2)=mdl.Coefficients.pValue(find_trials(mdl.Coefficients.Name,'right'));
end

% And then left vs right stroke
for nCh=1:length(layout.label)-2
    sub_table=SW_table_reorder(SW_table_reorder.Elec==layout.label{nCh},:);
    mdl=fitlme(sub_table,'SW_frequency~1+subGroupID+(1|SubID)');
    LxR_effect(nCh,1)=mdl.Coefficients.tStat(find_trials(mdl.Coefficients.Name,'right'));
    LxR_effect(nCh,2)=mdl.Coefficients.pValue(find_trials(mdl.Coefficients.Name,'right'));
end

cmap2=cbrewer('div','RdBu',64); cmap2=flipud(cmap2);

figure;
simpleTopoPlot_ft(LeftStroke_effect(:,1), layout,'on',[],0,1);
% ft_plot_lay_me(layout, 'chanindx', find(LeftStroke_effect(:,2)<0.05), 'pointsymbol','o','pointcolor','green','pointsize',36,'box','no','label','no')
% ft_plot_lay_me(layout, 'chanindx', find(LeftStroke_effect(:,2)<fdr(LeftStroke_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(LeftStroke_effect(:,1))))
title('Healthy Older x Left Stroke ~ SW Frequency');

figure;
simpleTopoPlot_ft(RightStroke_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(RightStroke_effect(:,2)<0.05), 'pointsymbol','o','pointcolor','green','pointsize',36,'box','no','label','no')
% ft_plot_lay_me(layout, 'chanindx', find(RightStroke_effect(:,2)<fdr(RightStroke_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(RightStroke_effect(:,1))))
title('Healthy Older x Right Stroke ~ SW Frequency');

figure;
simpleTopoPlot_ft(LxR_effect(:,1), layout,'on',[],0,1); hold on;
% ft_plot_lay_me(layout, 'chanindx', find(LxR_effect(:,2)<0.05), 'pointsymbol','o','pointcolor','green','pointsize',36,'box','no','label','no')
% ft_plot_lay_me(layout, 'chanindx', find(LxR_effect(:,2)<fdr(LxR_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(LxR_effect(:,1))))
title('Left x Right Stroke ~ SW Frequency');

%% Upslope
LeftStroke_effect=[]; RightStroke_effect=[];
for nCh=1:length(layout.label)-2
    sub_table=SW_table(SW_table.Elec==layout.label{nCh},:);
    mdl=fitlme(sub_table,'SW_upslope~1+subGroupID+(1|SubID)');
    LeftStroke_effect(nCh,1)=mdl.Coefficients.tStat(find_trials(mdl.Coefficients.Name,'left'));
    LeftStroke_effect(nCh,2)=mdl.Coefficients.pValue(find_trials(mdl.Coefficients.Name,'left'));
    RightStroke_effect(nCh,1)=mdl.Coefficients.tStat(find_trials(mdl.Coefficients.Name,'right'));
    RightStroke_effect(nCh,2)=mdl.Coefficients.pValue(find_trials(mdl.Coefficients.Name,'right'));
end

% And then left vs right stroke
for nCh=1:length(layout.label)-2
    sub_table=SW_table_reorder(SW_table_reorder.Elec==layout.label{nCh},:);
    mdl=fitlme(sub_table,'SW_upslope~1+subGroupID+(1|SubID)');
    LxR_effect(nCh,1)=mdl.Coefficients.tStat(find_trials(mdl.Coefficients.Name,'right'));
    LxR_effect(nCh,2)=mdl.Coefficients.pValue(find_trials(mdl.Coefficients.Name,'right'));
end

cmap2=cbrewer('div','RdBu',64); cmap2=flipud(cmap2);

figure;
simpleTopoPlot_ft(LeftStroke_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(LeftStroke_effect(:,2)<0.05), 'pointsymbol','o','pointcolor','green','pointsize',36,'box','no','label','no')
% ft_plot_lay_me(layout, 'chanindx', find(LeftStroke_effect(:,2)<fdr(LeftStroke_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(LeftStroke_effect(:,1))))
title('Healthy Older x Left Stroke ~ SW Upward Slope');

figure;
simpleTopoPlot_ft(RightStroke_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(RightStroke_effect(:,2)<0.05), 'pointsymbol','o','pointcolor','green','pointsize',36,'box','no','label','no')
% ft_plot_lay_me(layout, 'chanindx', find(RightStroke_effect(:,2)<fdr(RightStroke_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(RightStroke_effect(:,1))))
title('Healthy Older x Right Stroke ~ SW Upward Slope');

figure;
simpleTopoPlot_ft(LxR_effect(:,1), layout,'on',[],0,1); hold on;
% ft_plot_lay_me(layout, 'chanindx', find(LxR_effect(:,2)<0.05), 'pointsymbol','o','pointcolor','green','pointsize',36,'box','no','label','no')
% ft_plot_lay_me(layout, 'chanindx', find(LxR_effect(:,2)<fdr(LxR_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(LxR_effect(:,1))))
title('Left x Right Stroke ~ SW Upward Slope');

%% Downslope
LeftStroke_effect=[]; RightStroke_effect=[];
for nCh=1:length(layout.label)-2
    sub_table=SW_table(SW_table.Elec==layout.label{nCh},:);
    mdl=fitlme(sub_table,'SW_downslope~1+subGroupID+(1|SubID)');
    LeftStroke_effect(nCh,1)=mdl.Coefficients.tStat(find_trials(mdl.Coefficients.Name,'left'));
    LeftStroke_effect(nCh,2)=mdl.Coefficients.pValue(find_trials(mdl.Coefficients.Name,'left'));
    RightStroke_effect(nCh,1)=mdl.Coefficients.tStat(find_trials(mdl.Coefficients.Name,'right'));
    RightStroke_effect(nCh,2)=mdl.Coefficients.pValue(find_trials(mdl.Coefficients.Name,'right'));
end

% And then left vs right stroke
for nCh=1:length(layout.label)-2
    sub_table=SW_table_reorder(SW_table_reorder.Elec==layout.label{nCh},:);
    mdl=fitlme(sub_table,'SW_downslope~1+subGroupID+(1|SubID)');
    LxR_effect(nCh,1)=mdl.Coefficients.tStat(find_trials(mdl.Coefficients.Name,'right'));
    LxR_effect(nCh,2)=mdl.Coefficients.pValue(find_trials(mdl.Coefficients.Name,'right'));
end

cmap2=cbrewer('div','RdBu',64); cmap2=flipud(cmap2);

figure;
simpleTopoPlot_ft(LeftStroke_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(LeftStroke_effect(:,2)<0.05), 'pointsymbol','o','pointcolor','green','pointsize',36,'box','no','label','no')
% ft_plot_lay_me(layout, 'chanindx', find(LeftStroke_effect(:,2)<fdr(LeftStroke_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(LeftStroke_effect(:,1))))
title('Healthy Older x Left Stroke ~ SW Downward Slope');

figure;
simpleTopoPlot_ft(RightStroke_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(RightStroke_effect(:,2)<0.05), 'pointsymbol','o','pointcolor','green','pointsize',36,'box','no','label','no')
% ft_plot_lay_me(layout, 'chanindx', find(RightStroke_effect(:,2)<fdr(RightStroke_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(RightStroke_effect(:,1))))
title('Healthy Older x Right Stroke ~ SW Downward Slope');

figure;
simpleTopoPlot_ft(LxR_effect(:,1), layout,'on',[],0,1); hold on;
% ft_plot_lay_me(layout, 'chanindx', find(LxR_effect(:,2)<0.05), 'pointsymbol','o','pointcolor','green','pointsize',36,'box','no','label','no')
% ft_plot_lay_me(layout, 'chanindx', find(LxR_effect(:,2)<fdr(LxR_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(LxR_effect(:,1))))
title('Left x Right Stroke ~ SW Downward Slope');
%}

%% Topos
%{
SW_density=[];
SW_dslope=[];
SW_amplitude=[];
for nCh=1:length(layout.label)-2
    sub_table=SW_table(SW_table.Elec==layout.label{nCh},:);
    SW_density(nCh,1)=nanmean(sub_table.SW_density(sub_table.subGroupID=='healthy_old'));
    SW_density(nCh,2)=nanmean(sub_table.SW_density(sub_table.subGroupID=='right_stroke'));
    SW_density(nCh,3)=nanmean(sub_table.SW_density(sub_table.subGroupID=='left_stroke'));

    SW_amplitude(nCh,1)=nanmean(sub_table.SW_amplitude(sub_table.subGroupID=='healthy_old'));
    SW_amplitude(nCh,2)=nanmean(sub_table.SW_amplitude(sub_table.subGroupID=='right_stroke'));
    SW_amplitude(nCh,3)=nanmean(sub_table.SW_amplitude(sub_table.subGroupID=='left_stroke'));

    SW_dslope(nCh,1)=nanmean(sub_table.SW_downslope(sub_table.subGroupID=='healthy_old'));
    SW_dslope(nCh,2)=nanmean(sub_table.SW_downslope(sub_table.subGroupID=='right_stroke'));
    SW_dslope(nCh,3)=nanmean(sub_table.SW_downslope(sub_table.subGroupID=='left_stroke'));
end
figure;
% subplot(3,3,1);
simpleTopoPlot_ft(SW_density(:,1), layout,'on',[],0,1);
% subplot(3,3,2);
figure;
simpleTopoPlot_ft(SW_density(:,2), layout,'on',[],0,1);
% subplot(3,3,3);
figure;
simpleTopoPlot_ft(SW_density(:,3), layout,'on',[],0,1);

subplot(3,3,4);
simpleTopoPlot_ft(SW_amplitude(:,1), layout,'on',[],0,1);
subplot(3,3,5);
simpleTopoPlot_ft(SW_amplitude(:,2), layout,'on',[],0,1);
subplot(3,3,6);
simpleTopoPlot_ft(SW_amplitude(:,3), layout,'on',[],0,1);

subplot(3,3,7);
simpleTopoPlot_ft(SW_dslope(:,1), layout,'on',[],0,1);
subplot(3,3,8);
simpleTopoPlot_ft(SW_dslope(:,2), layout,'on',[],0,1);
subplot(3,3,9);
simpleTopoPlot_ft(SW_amplitude(:,3), layout,'on',[],0,1);
%}