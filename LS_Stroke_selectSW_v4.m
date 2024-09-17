%%
close all
clear all
if isempty(findstr(pwd,'thandrillon'))==0
    path_data = '/Users/thandrillon/Data/StrokeData/EEG';
    path_save = '/Users/thandrillon/Data/StrokeData/SWdetection';
    path_LSCPtools='/Users/thandrillon/WorkGit/LSCPtools/';
    path_eeglab='/Users/thandrillon/WorkGit/projects/ext/eeglab/';
    path_TESA='/Users/thandrillon/WorkGit/projects/ext/TESA/';
    addpath(path_eeglab);
    addpath(path_TESA);
    addpath(genpath('/Users/thandrillon/WorkLocal/local/FastICA_25/'))
     path_temp=('/Users/thandrillon/WorkGit/projects/inprogress/LS_Stroke/');
elseif isempty(findstr(pwd,'Daniel'))==0
    path_data = '/fs04/so34/Daniel/Data/Raw';
    path_save = '/fs04/so34/Daniel/Data/SWdetection';
    path_LSCPtools = '/fs04/so34/LocalSleep/Stroke/Scripts';
    addpath(genpath('/fs04/so34/Daniel/Functions/eeglab2021.0/'));
end
path_fieldtrip='/Users/thandrillon/WorkGit/projects/ext/fieldtrip/';
addpath(path_fieldtrip)
addpath(genpath(path_LSCPtools));
IDs = {'HN996'; 'HN968'; 'HN969'; 'HN970'; 'HN971'; 'HN972'; 'HN973'; 'HN974'; 'HN976'; 'HN977'; ...
    'HN978'; 'HN980'; 'HN981'; 'HN982'; 'HN983'; 'HN985'; 'HN986'; 'HN987'; 'HN988'; 'HN989'; ...
    'HN990'; 'HN992'; 'HN993'; 'HN994'; 'HN995'; 'HN998'; 'HN999'; 'S002'; 'S003'; 'S004'; ...
    'S010'; 'S012'; 'S013'; 'S014'; 'S016';'S017'; 'S018'; 'S020'; 'S025'; 'S026'; ...
    'S027'; 'S029'; 'S030'; 'S031'; 'S032'; 'S033'; 'S103'; 'S104'; 'S107';'S109'; ...
    'S111'; 'S112'; 'S114'; 'S115'; 'S201'; 'S202'; 'S207' };

%%
redo=0;

SW_table=array2table(zeros(0,11),'VariableNames',{'SubID','GroupID','Elec','Block','SW_density','SW_amplitude','SW_frequency','SW_downslope','SW_upslope','SW_threshold','BlockDuration'});
SW_table.SubID=categorical(SW_table.SubID);
SW_table.GroupID=categorical(SW_table.GroupID);
SW_table.Elec=categorical(SW_table.Elec);

%     allSW_table=array2table(zeros(0,9),'VariableNames',{'SubID','GroupID','Elec','Block','SW_density','SW_amplitude','SW_frequency','SW_downslope','SW_upslope'});
%     allSW_table.SubID=categorical(allSW_table.SubID);
%     allSW_table.GroupID=categorical(allSW_table.GroupID);
%     allSW_table.Elec=categorical(allSW_table.Elec);
CommonChannels=[];
for idx = 1:length(IDs)
    EEG = []; ALLEEG = [];

    file_names=dir([path_data filesep IDs{idx} filesep IDs{idx} '*.mat']);
    if isempty(file_names)
        file_names=dir([path_data filesep IDs{idx} '*.mat']);
        if isempty(file_names)
            continue;
        end
    end
    merged_all_Waves=[];
    Duration=[];
    BlockInfo=[];
    for nF=1:length(file_names)
        FileName=file_names(nF).name;
        separators=findstr(FileName,'_');
        SubID=FileName(1:separators(1)-1);
        BlockID=str2num(FileName(separators(1)+6:separators(2)-1));
        GroupID=FileName(separators(2)+1:separators(end)-1);
        if BlockID>9
           fprintf('... %s - %2.0f - %s SKIPPING\n',SubID,BlockID,GroupID)
         continue;
        end
        if exist([path_save filesep 'allSW_trim_' file_names(nF).name])==0
            fprintf('... %s - %2.0f - %s MISSING\n',SubID,BlockID,GroupID)
            continue;
        end
        load([path_save filesep 'allSW_trim_' file_names(nF).name(1:end-4)]);
        merged_all_Waves=[merged_all_Waves ; double(all_Waves)];
        load([file_names(nF).folder filesep file_names(nF).name]);
        chan_labels={EEG_120Hz.chanlocs.labels};
        data=EEG_120Hz.data;
         %start from first stimulus
        if isnumeric([EEG_120Hz.event.type])==1
            stimEvents=find([EEG_120Hz.event.type]==4);
        else
            stimEvents=match_str({EEG_120Hz.event.type},'4');
        end
        this_Sart=EEG_120Hz.event(stimEvents(1)).latency; 
        data=data(:,this_Sart:end);
        Fs=EEG_120Hz.srate;
        Duration(nF)=size(data,2)/Fs/60;
        fprintf('... %s - %2.0f - %s - %g - %g\n',SubID,BlockID,GroupID,length(chan_labels),Duration(nF))
        if isempty(CommonChannels)
            CommonChannels=chan_labels;
        else
            CommonChannels=intersect(CommonChannels,chan_labels);
        end
        BlockInfo=[BlockInfo ; [nF BlockID Duration(nF) unique(all_Waves(:,2))]];

    end

    %%% Select slow waves
    %%%% parameters to clean SW and select top-amplitude SW
    paramSW.prticle_Thr=90; % percentile of slow waves selected
    paramSW.LimFrqW=[1 7]; % frequency range: delta [1 4] or theta [4 7]
    paramSW.AmpCriterionIdx=4; % Amplitude criterion: 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.fixThr=[]; % to select an absolute and not relative threshold
    paramSW.art_ampl=150; % threshold to discard too high amplitude events
    paramSW.max_posampl=75; % threshold for artifact on positive peak (blinks typically)
    paramSW.max_Freq=7; % max frequency considered

    all_freq=1./(abs((merged_all_Waves(:,5)-merged_all_Waves(:,7)))./Fs);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq<paramSW.LimFrqW(1) | all_freq>paramSW.LimFrqW(2) | all_freq>paramSW.max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(merged_all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(merged_all_Waves(:,11)>paramSW.max_posampl | merged_all_Waves(:,14)>paramSW.art_ampl| abs(merged_all_Waves(:,15))>paramSW.art_ampl)*100)
    merged_all_Waves(all_freq<paramSW.LimFrqW(1) | all_freq>paramSW.LimFrqW(2) | all_freq>paramSW.max_Freq | merged_all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | merged_all_Waves(:,11)>paramSW.max_posampl| merged_all_Waves(:,14)>paramSW.art_ampl| abs(merged_all_Waves(:,15))>paramSW.art_ampl,:)=[];

    slow_Waves=[];
    thr_Wave=[];
    for nE=1:length(chan_labels)
        thisE_Waves=merged_all_Waves(merged_all_Waves(:,3)==nE,:);
        temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);

        if ~isempty(paramSW.fixThr)
            thr_Wave(nE)=paramSW.fixThr;
        else
            thr_Wave(nE)=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
        end
        slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>thr_Wave(nE),:)];
    end
    save([path_save filesep 'SW_' IDs{idx}],'slow_Waves','Fs','chan_labels');

    for nB=1:size(BlockInfo,1)
        slow_Waves_perE=[];

        for nE=1:length(chan_labels)
            slow_Waves_perE=[slow_Waves_perE ; [sum(slow_Waves(slow_Waves(:,2)==BlockInfo(nB,4),3)==nE)/BlockInfo(nB,3) nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==BlockInfo(nB,4),4)) nanmean(1./((slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==BlockInfo(nB,4),7)-slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==BlockInfo(nB,4),5))/Fs)) ...
                nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==BlockInfo(nB,4),12)) nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==BlockInfo(nB,4),13))]];
        end

        table_length=size(SW_table,1);
        SW_table.SubID(table_length+(1:length(chan_labels)))=repmat({SubID},length(chan_labels),1);
        SW_table.GroupID(table_length+(1:length(chan_labels)))=repmat({GroupID},length(chan_labels),1);
        SW_table.Block(table_length+(1:length(chan_labels)))=repmat(BlockInfo(nB,4),length(chan_labels),1);
        SW_table.Elec(table_length+(1:length(chan_labels)))=chan_labels;
        SW_table.SW_density(table_length+(1:length(chan_labels)))=slow_Waves_perE(:,1);
        SW_table.SW_amplitude(table_length+(1:length(chan_labels)))=slow_Waves_perE(:,2);
        SW_table.SW_frequency(table_length+(1:length(chan_labels)))=slow_Waves_perE(:,3);
        SW_table.SW_downslope(table_length+(1:length(chan_labels)))=slow_Waves_perE(:,4);
        SW_table.SW_upslope(table_length+(1:length(chan_labels)))=slow_Waves_perE(:,5);
        SW_table.SW_threshold(table_length+(1:length(chan_labels)))=thr_Wave';
        SW_table.BlockDuration(table_length+(1:length(chan_labels)))=repmat(BlockInfo(nB,3) ,length(chan_labels),1);
    end
    %         table_length=size(allSW_table,1);
    %         allSW_table.SubID(table_length+(1:length(chan_labels)))=repmat({SubID},length(chan_labels),1);
    %         allSW_table.GroupID(table_length+(1:length(chan_labels)))=repmat({GroupID},length(chan_labels),1);
    %         allSW_table.Block(table_length+(1:length(chan_labels)))=repmat(BlockID,length(chan_labels),1);
    %         allSW_table.Elec(table_length+(1:length(chan_labels)))=chan_labels;
    %         allSW_table.SW_density(table_length+(1:length(chan_labels)))=allSlow_Waves_perE(:,1);
    %         allSW_table.SW_amplitude(table_length+(1:length(chan_labels)))=allSlow_Waves_perE(:,2);
    %         allSW_table.SW_frequency(table_length+(1:length(chan_labels)))=allSlow_Waves_perE(:,3);
    %         allSW_table.SW_downslope(table_length+(1:length(chan_labels)))=allSlow_Waves_perE(:,4);
    %         allSW_table.SW_upslope(table_length+(1:length(chan_labels)))=allSlow_Waves_perE(:,5);

end
%%
     A=load([path_temp 'SW_individualThreshold_withICA4_subGroups.mat']);
uniqueSubIDs=unique(SW_table.SubID);
SW_table.subGroupID=nan(size(SW_table,1),1);
SW_table.subGroupID=categorical(SW_table.subGroupID);
for nSub=1:length(uniqueSubIDs)
    SW_table.subGroupID(SW_table.SubID==uniqueSubIDs(nSub))=unique(A.SW_table.subGroupID(A.SW_table.SubID==uniqueSubIDs(nSub)));
end

writetable(SW_table,[path_save filesep 'SW_individualThreshold_acrossBlocks_trim.csv'])
% writetable(allSW_table,[path_save filesep 'SW_noThreshold.csv'])

%%
cfg = [];
cfg.channel = CommonChannels;
cfg.center = 'yes';
cfg.layout = 'biosemi64.lay';
layout = ft_prepare_layout(cfg);

%%
uniqueGroups=unique(SW_table.subGroupID);
figure;
    temp_topo=[];
for nG=1:length(uniqueGroups)
    subplot(1,length(uniqueGroups)+1,nG)
    for nCh=1:length(layout.label)-2
        temp_topo(nG,nCh)=nanmean(SW_table.SW_density((SW_table.Elec==layout.label{nCh}) & (SW_table.subGroupID==uniqueGroups(nG))));
    end
    simpleTopoPlot_ft(temp_topo(nG,:)', layout,'labels',[],0,1);
    colorbar;
    title(uniqueGroups(nG))
end
for nG=1:length(uniqueGroups)
    subplot(1,length(uniqueGroups)+1,nG)
caxis([min(min(temp_topo)) max(max(temp_topo))]);
end

    subplot(1,length(uniqueGroups)+1,length(uniqueGroups)+1)
    simpleTopoPlot_ft(temp_topo(3,:)'-temp_topo(2,:)', layout,'labels',[],0,1);
    colorbar;
caxis([-1 1]*max(max(abs(temp_topo(3,:)'-temp_topo(2,:)'))));

%%
figure;
    temp_topo=[];
for nG=1:length(uniqueGroups)
    subplot(1,length(uniqueGroups)+1,nG)
    for nCh=1:length(layout.label)-2
        temp_topo(nG,nCh)=nanmean(SW_table.SW_threshold((SW_table.Elec==layout.label{nCh}) & (SW_table.subGroupID==uniqueGroups(nG))));
    end
    simpleTopoPlot_ft(temp_topo(nG,:)', layout,'labels',[],0,1);
    colorbar;
end
for nG=1:length(uniqueGroups)
    subplot(1,length(uniqueGroups)+1,nG)
caxis([min(min(temp_topo)) max(max(temp_topo))]);
end
    subplot(1,length(uniqueGroups)+1,length(uniqueGroups)+1)
    simpleTopoPlot_ft(temp_topo(3,:)'-temp_topo(2,:)', layout,'labels',[],0,1);
    colorbar;
caxis([-1 1]*max(max(abs(temp_topo(3,:)'-temp_topo(2,:)'))));
