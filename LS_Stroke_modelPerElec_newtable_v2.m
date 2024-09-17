
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
     SW_table=readtable([path_temp 'SW_5uVThreshold_acrossBlocks_trim.csv']);
     A=load([path_temp 'SW_individualThreshold_withICA4_subGroups.mat']);
     commonChans=A.commonChans;
end
% Reorder for L vs R comparison
SW_table.SubID=categorical(SW_table.SubID);
SW_table.GroupID=categorical(SW_table.GroupID);
SW_table.Elec=categorical(SW_table.Elec);
SW_table.subGroupID=categorical(SW_table.subGroupID);

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
totperm=500;
for nCh=1:length(layout.label)-2
    sub_table=SW_table(SW_table.Elec==layout.label{nCh},:);
    if nCh==1
        out_pred_perm=[];
        [real_out, cont_out, perm_out, cont_perm_out, out_pred_perm]=lme_perm_lsstroke(sub_table,'subGroupID','SW_density~1+Block+pred+(1|SubID)',totperm);
    else
        [real_out, cont_out, perm_out, cont_perm_out, next_out_pred_perm]=lme_perm_lsstroke(sub_table,'subGroupID','SW_density~1+Block+pred+(1|SubID)',totperm,out_pred_perm);
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
cfg_neighb.channel=layout.label(1:end-2);
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


%% Effect of block
groupIDs=unique(SW_table.subGroupID);
figure; hold on;
for nCond=1:3
    temp_SWdensity=grpstats(SW_table.SW_density(SW_table.subGroupID==groupIDs(nCond)),SW_table.Block(SW_table.subGroupID==groupIDs(nCond)));
    plot(1:9,temp_SWdensity)
end
%%
SWdens_norm_est=cell(1,2);
for nCh=1:length(layout.label)-2
    sub_table=normSW_table(normSW_table.Elec==layout.label{nCh},:);
    if nCh==1
        out_pred_perm=[];
        [real_out, cont_out, perm_out, cont_perm_out, out_pred_perm]=lme_perm_lsstroke_rightleft(sub_table,'subGroupID','SW_density~1+Block+pred+(1|SubID)',totperm);
    else
        [real_out, cont_out, perm_out, cont_perm_out, next_out_pred_perm]=lme_perm_lsstroke_rightleft(sub_table,'subGroupID','SW_density~1+Block+pred+(1|SubID)',totperm,out_pred_perm);
    end
    SWdens_norm_est{1}=[SWdens_norm_est{1} ; [nCh real_out]];
    SWdens_norm_est{2}=[SWdens_norm_est{2} ; [nCh*ones(totperm,1) perm_out]];
end

%% Filter clusters
clus_alpha=0.1;
montecarlo_alpha=0.05;

cfg_neighb=[];
cfg_neighb.method = 'template';
cfg_neighb.layout='biosemi64.lay';
cfg_neighb.channel=layout.label(1:end-2);
neighbours = ft_prepare_neighbours(cfg_neighb);
neighbours(~ismember({neighbours.label},unique(SW_table.Elec)))=[];
[SWdens_norm_clus]=get_clusterperm_lme_lsstroke(SWdens_norm_est,clus_alpha,montecarlo_alpha,totperm,neighbours,1);


%
orderedConds={'right_stroke vs left_stroke norm.'};

cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);
limNumClus=1;
limMax=4;
ClustersByDrugs=cell(2,3);
Cond2={'RStr vs LStr norm'};
for nCond=1
    figure; format_fig; set(gcf,'Position',[1     1   420   420]);

    temp_topo=SWdens_norm_est{1}(SWdens_norm_est{1}(:,5)==nCond,3);
    temp_topo2=zeros(size(temp_topo));
    temp_topo3=zeros(size(temp_topo));
    temp_clus=SWdens_norm_clus{nCond};

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
    hb=colorbar('Position',[0.94    0.6    0.05    0.33]);
    colormap(cmap2);
    title(Cond2{nCond})
end
