cfg=[];
cfg.elec='LS_Stroke_ElecPos.sfp';
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
caxis([-1 1]*max(abs(Group_effect(:,1))))
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
subplot(3,2,1);
simpleTopoPlot_ft(SW_density(:,1), layout,'on',[],0,1);
subplot(3,2,2);
simpleTopoPlot_ft(SW_density(:,2), layout,'on',[],0,1);

subplot(3,2,3);
simpleTopoPlot_ft(SW_amplitude(:,1), layout,'on',[],0,1);
subplot(3,2,4);
simpleTopoPlot_ft(SW_amplitude(:,2), layout,'on',[],0,1);

subplot(3,2,5);
simpleTopoPlot_ft(SW_dslope(:,1), layout,'on',[],0,1);
subplot(3,2,6);
simpleTopoPlot_ft(SW_dslope(:,2), layout,'on',[],0,1);