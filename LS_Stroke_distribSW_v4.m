%%
% clos all
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
bins=0:2:100;
all_distrib_SW=[];
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

    for nB=1:size(BlockInfo,1)
        table_length=size(SW_table,1);
        SW_table.SubID(table_length+1)={SubID};
        SW_table.GroupID(table_length+1)={GroupID};
        SW_table.Block(table_length+1)=BlockInfo(nB,4);
        SW_table.BlockDuration(table_length+1)=BlockInfo(nB,3);

        this_amp=merged_all_Waves(:,4);
        amp_bins=histc(this_amp,bins);
        amp_bins=(amp_bins/BlockInfo(nB,3))/length(chan_labels)';

        all_distrib_SW=[all_distrib_SW ; amp_bins'];
    end
end

%%
figure;
uniqueGroupIDs={'healthy_old','neglect'};
cmap=cbrewer('qual','Set2',6);
hplot=[];
for nG=1:length(uniqueGroupIDs)
    temp_plot=all_distrib_SW(SW_table.GroupID==uniqueGroupIDs(nG),:);
    [pV hplot(nG)]=simpleTplot(bins,temp_plot,0,cmap(nG,:),[0],'-',0.5,1,0,1,4);
end
legend(hplot,uniqueGroupIDs)
format_fig;
xlabel('Amplitude all waves')
ylabel('SW density')
xlim([0 70])
