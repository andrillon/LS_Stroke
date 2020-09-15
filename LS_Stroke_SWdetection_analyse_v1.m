%%
clear all
close all

run LS_Stroke_localdef.m
% Participants HN989 and HN978 completed nine blocks of the Dots task at 90% coherence,
% with the naming convention "[Participant ID] _ [Block #]". Participants A101 and A102
% completed five blocks of the Dots task at 90% coherence, followed by five blocks of
% the Dots task at 25% coherence, with the naming convention "[Participant ID] [Coherence] [Block #]".
eeg_files=dir([path_data filesep '*' filesep '*.eeg']);

addpath(path_fieldtrip);
ft_defaults;

addpath(genpath(path_LSCPtools));

%%
all_slow_Waves=[];
table_SW=[];
for nF=1:length(eeg_files)
    
    filename=eeg_files(nF).name;
    fprintf('... ... load SW detection for %s %g/%g\n',filename,nF,length(eeg_files))
    load([path_data filesep 'SWdetection'  filesep 'SW_clean_90thr_P2P_allE_' filename(1:end-4)]); %,'slow_Waves','labels','Fs','paramSW');
    all_slow_Waves=[all_slow_Waves ; slow_Waves];
    
    header=ft_read_header([eeg_files(nF).folder filesep eeg_files(nF).name]);
    nout=hist(slow_Waves(:,3),size(labels,1));
    nout=nout./(header.nSamples/header.Fs/60);
    table_SW=[table_SW ; repmat(slow_Waves(1,1:2),size(labels,1),1) (1:size(labels,1))' nout'];
end

%%
% sens = ft_read_sens('LS_Stroke_ElecPos.sfp');
cfg=[];
cfg.elec='LS_Stroke_ElecPos.sfp';
layout = ft_prepare_layout(cfg);
cmap=cbrewer('seq','YlOrRd',256); % select a sequential colorscale from yellow to red (64)

temp_topo=grpstats(table_SW(:,4),table_SW(:,3));
simpleTopoPlot_ft(temp_topo, layout,'on', [], 0,1);
caxis([2 14])
colormap(cmap);
title('SW density (W/min)')
%%
figure;
for nBlock=1:5
    subplot(1,5,nBlock);
    temp_topo=grpstats(table_SW(table_SW(:,2)==nBlock,4),table_SW(table_SW(:,2)==nBlock,3));
    simpleTopoPlot_ft(temp_topo, layout,'on', [], 0,1);
    caxis([2 14])
    colormap(cmap);
    title({'SW density (W/min)',sprintf('Block %g',nBlock)})
end