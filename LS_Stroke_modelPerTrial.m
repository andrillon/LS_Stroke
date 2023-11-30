clc, clear all

addpath(genpath('/fs04/so34/Daniel/Script/Dependencies/'));
addpath(genpath('/fs04/so34/Daniel/Functions/LSCPtools/'));
path_temp=('/fs04/so34/Daniel/Data/SWdetection/');

load('/fs04/so34/Daniel/Script/Dependencies/chanlocs_Oxford.mat');
load('/fs04/so34/Daniel/Script/Dependencies/chanlocs_Monash.mat');
SW_table = readtable([path_temp 'SW_and_Behav_perTrial.csv']);
load([path_temp 'SW_individualThreshold_withICA4.mat'],'commonChans','order_Monash','order_Oxford');

%% Restructure table
SW_table2=SW_table; % Copy the table
SW_table2(~ismember(SW_table2.Elec,commonChans),:)=[]; % Remove rows where the electrode is not common across sites

% Group ID
GroupID2=zeros(height(SW_table2),1); % Create blank subGroupID variable
GroupID2(ismember(SW_table2.GroupID,'healthy_old'))=1; % Set rows with healthy older to 1
GroupID2(ismember(SW_table2.GroupID,'left_stroke'))=2; % Set rows with left stroke to 2
GroupID2(ismember(SW_table2.GroupID,'right_stroke'))=3; % Set rows with right stroke to 3

% Hemifield
% StimSide2=zeros(height(SW_table2),1); % Create blank subGroupID variable
% StimSide2(ismember(SW_table2.StimSide,'Left'))=1; % Set rows with left side to 1
% StimSide2(ismember(SW_table2.StimSide,'Right'))=2; % Set rows with right side to 2
% StimSide2(ismember(SW_table2.StimSide,'Undefined'))=3; % Set rows with undefined side to 3

% Elec2=zeros(height(SW_table2),1); % Create blank electrode variable
% for elec=1:59 % For all electrodes
%     currelec=commonChans(elec); % Choose one electrode from common chans
%     Elec2(ismember(SW_table2.Elec,currelec))=elec; % Set all rows with that electrode label to the same number
% end

% SubID
SubID2=zeros(height(SW_table2),1); % Create blank SubID variable
allsubs=unique(SW_table2.SubID); % Get all SubIDs
for sub=1:57 % For all subjects
    currsub=allsubs(sub); % Choose one SubID from list
    SubID2(ismember(SW_table2.SubID,currsub))=sub; % Set all rows with that participant label to the same number
end

% Add new variables to table
SW_table2.SubID2=SubID2; % Add new sub IDs to table
SW_table2.SubID2=categorical(SW_table2.SubID2); % Convert new sub IDs to cat
% SW_table2.StimSide2=StimSide2; % Add new stim sides to table
% SW_table2.StimSide2=categorical(SW_table2.StimSide2); % Convert new stim sides to cat
SW_table2.GroupID2=GroupID2; % Add new group IDs to table
SW_table2.GroupID2=categorical(SW_table2.GroupID2); % Convert new group IDs to cat
SW_table2.Elec=categorical(SW_table2.Elec);
% SW_table2.Elec2=Elec2; % Add new elec IDs to table
% SW_table2.Elec2=categorical(SW_table2.Elec2); % Convert new sub IDs to cat
SW_table2.SW = categorical(SW_table2.SW);

% Reorder table & remove old variables
SW_table2 = movevars(SW_table2,"SubID2",'Before',"SubID"); 
% SW_table2 = movevars(SW_table2,"StimSide2",'Before',"StimSide");
SW_table2 = movevars(SW_table2,"GroupID2",'Before',"GroupID");
% SW_table2 = movevars(SW_table2,"Elec2",'Before',"Elec");
SW_table2.SubID = [];
% SW_table2.StimSide = [];
SW_table2.GroupID = [];
% SW_table2.Elec = [];
clearvars SW_table;

%
% SW_table2(ismember(SW_table2.StimSide2,'3'),:) = [];
SW_table2(SW_table2.Response==0,:) = [];
SW_table2(SW_table2.RT>3000,:) = [];
SW_table2(SW_table2.RT<150,:) = [];

% Reorder for L vs R comparison
% SW_table_reorder=SW_table;
% SW_table_reorder.GroupID=reordercats(SW_table_reorder.GroupID,{'left_stroke','right_stroke','healthy_old'});

%% Modelling
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

% Need to do separately per group
groups=unique(SW_table2.GroupID2);
for grp = 1:length(groups)
    for nCh=1:length(layout.label)-2
        sub_table=SW_table2(SW_table2.Elec==layout.label{nCh} & SW_table2.GroupID2==groups(grp),:);
        mdl=fitlme(sub_table,'RT~SW*StimSide+Block+(1|SubID2)');
        output=anova(mdl);
        BlockEffect(nCh,1)=output{2,2};
        BlockEffect(nCh,2)=output{2,5};
        HemifieldEffect(nCh,1)=output{3,2};
        HemifieldEffect(nCh,2)=output{3,5};
        SWEffect(nCh,1)=output{4,2};
        SWEffect(nCh,2)=output{4,5};
        SWxHemifieldEffect(nCh,1)=output{5,2};
        SWxHemifieldEffect(nCh,2)=output{5,5};
    end
    BlockLME{grp}=BlockEffect;
    HemifieldLME{grp}=HemifieldEffect;
    SWLME{grp}=SWEffect;
    SWxHemifieldLME{grp}=SWxHemifieldEffect;
end

%% Topoplots - RT ~ SW x Hemifield
cmap2=cbrewer('div','RdBu',64); cmap2=flipud(cmap2);

figure;
simpleTopoPlot_ft(SWxHemifieldLME{1}(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(SWxHemifieldLME{1}(:,2)<0.05), 'pointsymbol','o','pointcolor','green','pointsize',36,'box','no','label','no')
% ft_plot_lay_me(layout, 'chanindx', find(SWxHemifieldLME{1}(:,2)<fdr(SWxHemifieldLME{1}(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(SWxHemifieldLME{1}(:,1))))
title('Healthy: RT ~ SW x Hemifield');

figure;
simpleTopoPlot_ft(SWxHemifieldLME{2}(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(SWxHemifieldLME{2}(:,2)<0.05), 'pointsymbol','o','pointcolor','green','pointsize',36,'box','no','label','no')
% ft_plot_lay_me(layout, 'chanindx', find(SWxHemifieldLME{2}(:,2)<fdr(SWxHemifieldLME{2}(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(SWxHemifieldLME{2}(:,1))))
title('Left Stroke: RT ~ SW x Hemifield');

figure;
simpleTopoPlot_ft(SWxHemifieldLME{3}(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(SWxHemifieldLME{3}(:,2)<0.05), 'pointsymbol','o','pointcolor','green','pointsize',36,'box','no','label','no')
% ft_plot_lay_me(layout, 'chanindx', find(SWxHemifieldLME{3}(:,2)<fdr(SWxHemifieldLME{3}(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(SWxHemifieldLME{3}(:,1))))
title('Right Stroke: RT ~ SW x Hemifield');

%% Topoplots - RT ~ SW
cmap2=cbrewer('div','RdBu',64); cmap2=flipud(cmap2);

figure;
simpleTopoPlot_ft(SWLME{1}(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(SWLME{1}(:,2)<0.05), 'pointsymbol','o','pointcolor','green','pointsize',36,'box','no','label','no')
% ft_plot_lay_me(layout, 'chanindx', find(SWLME{1}(:,2)<fdr(SWLME{1}(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(SWLME{1}(:,1))))
title('Healthy: RT ~ SW');

figure;
simpleTopoPlot_ft(SWLME{2}(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(SWLME{2}(:,2)<0.05), 'pointsymbol','o','pointcolor','green','pointsize',36,'box','no','label','no')
% ft_plot_lay_me(layout, 'chanindx', find(SWLME{2}(:,2)<fdr(SWLME{2}(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(SWLME{2}(:,1))))
title('Left Stroke: RT ~ SW');

figure;
simpleTopoPlot_ft(SWLME{3}(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(SWLME{3}(:,2)<0.05), 'pointsymbol','o','pointcolor','green','pointsize',36,'box','no','label','no')
ft_plot_lay_me(layout, 'chanindx', find(SWLME{3}(:,2)<fdr(SWLME{3}(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar;
colormap(cmap2);
caxis([-1 1]*max(abs(SWLME{3}(:,1))))
title('Right Stroke: RT ~ SW');