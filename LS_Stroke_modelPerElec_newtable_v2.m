
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
     addpath(('/Users/thandrillon/WorkGit/projects/ext/fieldtrip'));
     ft_defaults;
     addpath(genpath('/Users/thandrillon/WorkGit/LSCPtools/'));
     path_temp=('/Users/thandrillon/WorkGit/projects/inprogress/LS_Stroke/');
     
     load('/Users/thandrillon/WorkGit/projects/inprogress/LS_Stroke/chanlocs_Oxford.mat');
     load('/Users/thandrillon/WorkGit/projects/inprogress/LS_Stroke/chanlocs_Monash.mat');
     SW_table=readtable([path_temp 'SW_individualThreshold.csv']);
     A=load([path_temp 'SW_individualThreshold_withICA4_subGroups.mat']);
     commonChans=A.commonChans;
end
% Reorder for L vs R comparison
SW_table.SubID=categorical(SW_table.SubID);
SW_table.GroupID=categorical(SW_table.GroupID);
SW_table.Elec=categorical(SW_table.Elec);
SW_table(SW_table.Block==0,:)=[];

uniqueSubIDs=unique(SW_table.SubID);
SW_table.subGroupID=nan(size(SW_table,1),1);
SW_table.subGroupID=categorical(SW_table.subGroupID);
for nSub=1:length(uniqueSubIDs)
    SW_table.subGroupID(SW_table.SubID==uniqueSubIDs(nSub))=unique(A.SW_table.subGroupID(A.SW_table.SubID==uniqueSubIDs(nSub)));
end
SW_table_reorder=SW_table;
SW_table_reorder.subGroupID=reordercats(SW_table_reorder.subGroupID,{'left_stroke','right_stroke','healthy_old'});

%%
normSW_table=SW_table_reorder;
normSW_table(normSW_table.GroupID=='healthy_old',:)=[];
normSW_table.GroupID=removecats(normSW_table.GroupID);
normSW_table.subGroupID=removecats(normSW_table.subGroupID);
uniqueSubIDs=unique(normSW_table.SubID);
uniqueBlockIDs=unique(normSW_table.Block);
for nSub=1:length(uniqueSubIDs)
    for nBlock=1:length(uniqueBlockIDs)
        sub_normSW_table=normSW_table(normSW_table.SubID==uniqueSubIDs(nSub) & normSW_table.Block==uniqueBlockIDs(nBlock),:);
        if ~isempty(sub_normSW_table)
            sub_normSW_table.SW_density=zscore(sub_normSW_table.SW_density);
            normSW_table(normSW_table.SubID==uniqueSubIDs(nSub) & normSW_table.Block==uniqueBlockIDs(nBlock),:)=sub_normSW_table;
        end
    end
end
%%
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

%%
% Density
LeftStroke_effect=[]; RightStroke_effect=[];
% mdl=fitlme(SW_table,'SW_density~subGroupID*Elec + (1|SubID)'); % Mark
SWdens_est=cell(1,2);
totperm=1;
for nCh=1:length(layout.label)-2
    sub_table=SW_table(SW_table.Elec==layout.label{nCh},:);
    mdl=fitlme(sub_table,'SW_density~1+subGroupID+(1|SubID)');
    if nCh==1
        out_pred_perm=[];
        [real_out, cont_out, perm_out, cont_perm_out, out_pred_perm]=lme_perm_lsstroke(sub_table,'subGroupID','SW_density~1+Block+pred+(1+Block|SubID)',totperm);
    else
        [real_out, cont_out, perm_out, cont_perm_out, next_out_pred_perm]=lme_perm_lsstroke(sub_table,'subGroupID','SW_density~1+Block+pred+(1+Block|SubID)',totperm,out_pred_perm);
    end
    SWdens_est{1}=[SWdens_est{1} ; [nCh*ones(3,1) real_out]];
    SWdens_est{2}=[SWdens_est{2} ; [nCh*ones(totperm*3,1) perm_out]];
end

%% Filter clusters
clus_alpha=0.05;
montecarlo_alpha=0.05;

cfg_neighb=[];
cfg_neighb.method = 'template';
cfg_neighb.layout='biosemi64.lay';
cfg.channel=unique(SW_table.Elec);
neighbours = ft_prepare_neighbours(cfg_neighb);
neighbours(~ismember({neighbours.label},unique(SW_table.Elec)))=[];
[SWdens_clus]=get_clusterperm_lme_lsstroke(SWdens_est,clus_alpha,montecarlo_alpha,totperm,neighbours,1);

%%
orderedConds={'left_stroke vs healthy_old','right_stroke vs healthy_old','right_stroke vs left_stroke'};

cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);
limNumClus=1;
limMax=8;
ClustersByDrugs=cell(2,3);
Cond2={'LStr vs CTL','RStr vs CTL','RStr vs LStr'};
for nCond=1:3
    figure; format_fig; set(gcf,'Position',[1     1   420   420]);
    [ha pos]=tight_subplot(1,1,0.02,0.05,0.05);
    hs=subplot(1,1,1); format_fig;
    set(hs,'Position',pos);
    
    temp_topo=SWdens_est{1}(SWdens_est{1}(:,5)==nCond,3);
    temp_topo2=zeros(size(temp_topo));
    temp_topo3=zeros(size(temp_topo));
    temp_clus=SWdens_clus{nCond};
    
    if ~isempty(temp_clus)
        for nclus=1:length(temp_clus)
            if length(match_str(layout.label,temp_clus{nclus}{2}))<limNumClus
                continue;
            end
            %             ft_plot_lay_me(layout, 'chanindx',match_str(layout.label,temp_clus{nclus}{2}),'pointsymbol','o','pointcolor','r','pointsize',64,'box','no','label','yes')
            temp_topo2(match_str(layout.label,temp_clus{nclus}{2}))=temp_topo(match_str(layout.label,temp_clus{nclus}{2}));
            temp_topo3(match_str(layout.label,temp_clus{nclus}{2}))=1;
            fprintf('... ... found %s cluster (%g) of %g electrodes (tval cluster=%g, Pmc=%g)\n',temp_clus{nclus}{1},nclus,length(temp_clus{nclus}{2}),temp_clus{nclus}{3},temp_clus{nclus}{4})
            if strcmp(temp_clus{nclus}{1},'pos')
                ClustersByDrugs{2,nCond}=[ClustersByDrugs{2,nCond} ; temp_clus{nclus}{2}];
            elseif strcmp(temp_clus{nclus}{1},'neg')
                ClustersByDrugs{1,nCond}=[ClustersByDrugs{1,nCond} ; temp_clus{nclus}{2}];
            end
        end
    end
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    %     ft_plot_lay_me(layout, 'chanindx',1:length(layout.label)-2,'pointsymbol','o','pointcolor',[1 1 1]*0.7,'pointsize',6,'box','no','label','no')
    format_fig;
    caxis([-1 1]*limMax)
    if ~isempty(temp_clus)
        for nclus=1:length(temp_clus)
            if length(match_str(layout.label,temp_clus{nclus}{2}))<limNumClus
                continue;
            end
            ft_plot_lay_me(layout, 'chanindx',match_str(layout.label,temp_clus{nclus}{2}),'pointsymbol','o','pointcolor','k','pointsize',96,'box','no','label','no')
        end
    end
    if nCond==3
        hb=colorbar('Position',[0.94    0.6    0.05    0.33]);
        
    end
    colormap(cmap2);
    
    title(Cond2{nCond})
end

%%
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


%%
% And then left vs right stroke
TopoLeft=[];
TopoRight=[];
for nCh=1:length(layout.label)-2
    sub_table=normSW_table(normSW_table.Elec==layout.label{nCh},:);
    mdl=fitlme(sub_table,'SW_density~1+subGroupID+(1|SubID)');
    LxR_effect(nCh,1)=mdl.Coefficients.tStat(find_trials(mdl.Coefficients.Name,'right'));
    LxR_effect(nCh,2)=mdl.Coefficients.pValue(find_trials(mdl.Coefficients.Name,'right'));

    TopoLeft(nCh)=mean(normSW_table.SW_density(normSW_table.Elec==layout.label{nCh} & normSW_table.subGroupID=='left_stroke'));
    TopoRight(nCh)=mean(normSW_table.SW_density(normSW_table.Elec==layout.label{nCh} & normSW_table.subGroupID=='right_stroke'));
end

cmap2=cbrewer('div','RdBu',64); cmap2=flipud(cmap2);
figure;
subplot(1,3,1)
simpleTopoPlot_ft(TopoLeft, layout,'on',[],0,1); hold on;
colorbar; caxis([-1 1]*0.7)
title('Left SW density norm');
subplot(1,3,2)
simpleTopoPlot_ft(TopoRight, layout,'on',[],0,1); hold on;
colorbar; caxis([-1 1]*0.7)
title('Right SW density norm');
subplot(1,3,3)
% colorbar; caxis([-1 1]*0.7)
% 
% figure;
simpleTopoPlot_ft(LxR_effect(:,1), layout,'on',[],0,1); hold on;
ft_plot_lay_me(layout, 'chanindx', find(LxR_effect(:,2)<0.05), 'pointsymbol','o','pointcolor','green','pointsize',36,'box','no','label','no')
ft_plot_lay_me(layout, 'chanindx', find(LxR_effect(:,2)<fdr(LxR_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(LxR_effect(:,1))))
title('Left vs Right');

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