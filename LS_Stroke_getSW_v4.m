%%
close all
clear all

if isempty(findstr(pwd,'thandrillon'))==0
    path_data = '/Users/thandrillon/Data/StrokeData/EEG';
    path_save = '/Users/thandrillon/Data/StrokeData/SWdetection';
    path_LSCPtools='/Users/thandrillon/WorkGit/LSCPtools/';
    path_fieldtrip='/Users/thandrillon/WorkGit/projects/ext/fieldtrip/';
elseif isempty(findstr(pwd,'Daniel'))==0
    path_data = '/fs04/so34/Daniel/Data/Raw';
    path_save = '/fs04/so34/Daniel/SWdetection';
    path_LSCPtools = '/fs04/so34/LocalSleep/Stroke/Scripts';
end
addpath(path_fieldtrip)
ft_defaults;
addpath(genpath(path_LSCPtools));
IDs = {'HN996'; 'HN968'; 'HN969'; 'HN970'; 'HN971'; 'HN972'; 'HN973'; 'HN974'; 'HN976'; 'HN977'; ...
    'HN978'; 'HN980'; 'HN9881'; 'HN982'; 'HN983'; 'HN985'; 'HN986'; 'HN987'; 'HN988'; 'HN989'; ...
    'HN990'; 'HN992'; 'HN993'; 'HN994'; 'HN995'; 'HN998'; 'HN999'; 'S002'; 'S003'; 'S004'; ...
    'S009'; 'S010'; 'S012'; 'S013'; 'S014'; 'S016';'S017'; 'S018'; 'S020'; 'S021'; ...
    'S025'; 'S026'; 'S027'; 'S029'; 'S030'; 'S031'; 'S032'; 'S033'; 'S103'; 'S104'; ...
    'S107';'S109'; 'S111'; 'S112'; 'S114'; 'S115'; 'S201'; 'S202'; 'S207' };

%%%%%% Events description
%      4: central fixation
%      5: random motion
%     12: response
%     13: ???
%     28: fixation break
% [101 105 109]; % left patch, up motion
% [102 106 110]; % left patch, down motion
% [103 107 111]; % right patch, up motion
% [104 108 112]; % right patch, down motion

%%
redo=1;

SW_table=array2table(zeros(0,9),'VariableNames',{'SubID','GroupID','Elec','Block','SW_density','SW_amplitude','SW_frequency','SW_downslope','SW_upslope'});
SW_table.SubID=categorical(SW_table.SubID);
SW_table.GroupID=categorical(SW_table.GroupID);
SW_table.Elec=categorical(SW_table.Elec);

%     allSW_table=array2table(zeros(0,9),'VariableNames',{'SubID','GroupID','Elec','Block','SW_density','SW_amplitude','SW_frequency','SW_downslope','SW_upslope'});
%     allSW_table.SubID=categorical(allSW_table.SubID);
%     allSW_table.GroupID=categorical(allSW_table.GroupID);
%     allSW_table.Elec=categorical(allSW_table.Elec);

for idx = 1:length(IDs)
    file_names=dir([path_data filesep IDs{idx} '*.mat']);
    if isempty(file_names)
        continue;
    end
    if exist([path_save filesep 'allSW_Feb2024_' IDs{idx} '.mat'])==0 || redo==1
        all_Waves=[];
        for nF=1:length(file_names)
            %%% load data
            load([file_names(nF).folder filesep file_names(nF).name]);
            chan_labels={EEG_120Hz.chanlocs.labels};
            data=EEG_120Hz.data;
            Fs=EEG_120Hz.srate;
            Duration=size(data,2)/Fs;
            FileName=file_names(nF).name;
            separators=findstr(FileName,'_');
            SubID=FileName(1:separators(1)-1);
            BlockID=str2num(FileName(separators(2)-1));
            GroupID=FileName(separators(2)+1:separators(end)-1);

            %%% Preprocess
            data=data-repmat(mean(data(match_str(chan_labels,{'TP7','TP8'}),:),1),size(data,1),1);
            data=data-mean(data,2);

            %%% Detect all slow waves
            [twa_results]=twalldetectnew_TA_v2(data,Fs,0);
            for nE=1:size(data,1)
                all_Waves=[all_Waves ; [repmat([nF BlockID nE],length(abs(cell2mat(twa_results.channels(nE).maxnegpkamp))),1) abs(cell2mat(twa_results.channels(nE).maxnegpkamp))'+abs(cell2mat(twa_results.channels(nE).maxpospkamp))' ...
                    cell2mat(twa_results.channels(nE).negzx)' ...
                    cell2mat(twa_results.channels(nE).poszx)' ...
                    cell2mat(twa_results.channels(nE).wvend)' ...
                    cell2mat(twa_results.channels(nE).maxnegpk)' ...
                    cell2mat(twa_results.channels(nE).maxnegpkamp)' ...
                    cell2mat(twa_results.channels(nE).maxpospk)' ...
                    cell2mat(twa_results.channels(nE).maxpospkamp)' ...
                    cell2mat(twa_results.channels(nE).mxdnslp)' ...
                    cell2mat(twa_results.channels(nE).mxupslp)' ...
                    cell2mat(twa_results.channels(nE).maxampwn)' ...
                    cell2mat(twa_results.channels(nE).minampwn)' ...
                    ]];
                % Columns of the all_Waves matrix
                %  1: subject number (as in the for loop), you could replace with a subject ID if numerical
                %  2: 0 here but could be a block number
                %  3: electrode number (same order as the data)
                %  4: peak-to-peak amplitude
                %  5: start (in sample)
                %  6: half wave position
                %  7: end
                %  8: position of negative peak
                %  9: amplitude of negative peak
                % 10: position of positive peak
                % 11: amplitude of positive peak
                % 12: maximum downward slope
                % 13: maximum upward slope
                % 14: max amplitude on the slow wave window
                % 15: min amplitude on the slow wave window
            end
        end
        save([path_save filesep 'allSW_Feb2024_' IDs{idx}],'all_Waves','Fs','chan_labels');
    else
        load([path_save filesep 'allSW_Feb2024_' IDs{idx}]);
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

    all_Waves=double(all_Waves);
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./Fs);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq<paramSW.LimFrqW(1) | all_freq>paramSW.LimFrqW(2) | all_freq>paramSW.max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>paramSW.max_posampl | all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl)*100)
    all_Waves(all_freq<paramSW.LimFrqW(1) | all_freq>paramSW.LimFrqW(2) | all_freq>paramSW.max_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];

    slow_Waves=[];
    thr_Wave=[];
    for nE=1:length(chan_labels)
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);

        if ~isempty(paramSW.fixThr)
            thr_Wave(nE)=paramSW.fixThr;
        else
            thr_Wave(nE)=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
        end
        slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>thr_Wave(nE),:)];
    end
    save([path_save filesep 'SW_Feb2024_' IDs{idx}],'slow_Waves','Fs','chan_labels');

    for nF=1:length(file_names)
        %%% load data
        load([file_names(nF).folder filesep file_names(nF).name]);
        data=EEG_120Hz.data;
        Fs=EEG_120Hz.srate;
        Duration=size(data,2)/Fs/60;
        FileName=file_names(nF).name;
        separators=findstr(FileName,'_');
        SubID=FileName(1:separators(1)-1);
        BlockID=str2num(FileName(separators(2)-1));
        GroupID=FileName(separators(2)+1:separators(end)-1);

        slow_Waves_perE=[];
        %         allSlow_Waves_perE=[];
        for nE=1:length(chan_labels)
            slow_Waves_perE=[slow_Waves_perE ; [sum(slow_Waves(slow_Waves(:,2)==BlockID,3)==nE)/Duration nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==BlockID,4)) nanmean(1./((slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==BlockID,7)-slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==BlockID,5))/Fs)) ...
                nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==BlockID,12)) nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==BlockID,13))]];

            %             allSlow_Waves_perE=[allSlow_Waves_perE ; [sum(all_Waves(:,3)==nE)/Duration nanmean(all_Waves(all_Waves(:,3)==nE,4)) nanmean(1./((all_Waves(all_Waves(:,3)==nE,7)-all_Waves(all_Waves(:,3)==nE,5))/Fs)) ...
            %                 nanmean(all_Waves(all_Waves(:,3)==nE,12)) nanmean(all_Waves(all_Waves(:,3)==nE,13))]];
        end

        table_length=size(SW_table,1);
        SW_table.SubID(table_length+(1:length(chan_labels)))=repmat({SubID},length(chan_labels),1);
        SW_table.GroupID(table_length+(1:length(chan_labels)))=repmat({GroupID},length(chan_labels),1);
        SW_table.Block(table_length+(1:length(chan_labels)))=repmat(BlockID,length(chan_labels),1);
        SW_table.Elec(table_length+(1:length(chan_labels)))=chan_labels;
        SW_table.SW_density(table_length+(1:length(chan_labels)))=slow_Waves_perE(:,1);
        SW_table.SW_amplitude(table_length+(1:length(chan_labels)))=slow_Waves_perE(:,2);
        SW_table.SW_frequency(table_length+(1:length(chan_labels)))=slow_Waves_perE(:,3);
        SW_table.SW_downslope(table_length+(1:length(chan_labels)))=slow_Waves_perE(:,4);
        SW_table.SW_upslope(table_length+(1:length(chan_labels)))=slow_Waves_perE(:,5);
        SW_table.SW_threshold(table_length+(1:length(chan_labels)))=thr_Wave';
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
end
writetable(SW_table,[path_save filesep 'SW_individualThreshold_acrossBlocks_Feb2024.csv'])
% writetable(allSW_table,[path_save filesep 'SW_noThreshold.csv'])

%%
cfg = [];
cfg.channel = cellstr(unique(SW_table.Elec));
cfg.center = 'yes';
cfg.layout = 'biosemi64.lay';
layout = ft_prepare_layout(cfg);

%%
uniqueGroups=unique(SW_table.GroupID);
figure;
    temp_topo=[];
for nG=1:length(uniqueGroups)
    subplot(1,length(uniqueGroups),nG)
    temp_topo=[];
    for nCh=1:length(layout.label)-2
        temp_topo(nG,nCh)=nanmean(SW_table.SW_density((SW_table.Elec==layout.label{nCh}) & (SW_table.GroupID==uniqueGroups(nG))));
    end
    simpleTopoPlot_ft(temp_topo(nG,:)', layout,'labels',[],0,1);
    colorbar;
end
for nG=1:length(uniqueGroups)
    subplot(1,length(uniqueGroups),nG)
caxis([min(min(temp_topo)) max(max(temp_topo))]);
end

%%
figure;
    temp_topo=[];
for nG=1:length(uniqueGroups)
    subplot(1,length(uniqueGroups),nG)
    for nCh=1:length(layout.label)-2
        temp_topo(nG,nCh)=nanmean(SW_table.SW_threshold((SW_table.Elec==layout.label{nCh}) & (SW_table.GroupID==uniqueGroups(nG))));
    end
    simpleTopoPlot_ft(temp_topo(nG,:)', layout,'labels',[],0,1);
    colorbar;
end
for nG=1:length(uniqueGroups)
    subplot(1,length(uniqueGroups),nG)
caxis([min(min(temp_topo)) max(max(temp_topo))]);
end