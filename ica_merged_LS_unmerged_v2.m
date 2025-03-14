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
elseif isempty(findstr(pwd,'Daniel'))==0
    path_data = '/fs04/so34/Daniel/Data/Raw';
    path_save = '/fs04/so34/Daniel/Data/SWdetection';
    path_LSCPtools = '/fs04/so34/LocalSleep/Stroke/Scripts';
addpath(genpath('/fs04/so34/Daniel/Functions/eeglab2021.0/'));
end
eeglab;
addpath(genpath(path_LSCPtools));
IDs = {'HN996'; 'HN968'; 'HN969'; 'HN970'; 'HN971'; 'HN972'; 'HN973'; 'HN974'; 'HN976'; 'HN977'; ...
    'HN978'; 'HN980'; 'HN981'; 'HN982'; 'HN983'; 'HN985'; 'HN986'; 'HN987'; 'HN988'; 'HN989'; ...
    'HN990'; 'HN992'; 'HN993'; 'HN994'; 'HN995'; 'HN998'; 'HN999'; 'S002'; 'S003'; 'S004'; ...
    'S010'; 'S012'; 'S013'; 'S014'; 'S016';'S017'; 'S018'; 'S020'; 'S025'; 'S026'; ...
    'S027'; 'S029'; 'S030'; 'S031'; 'S032'; 'S033'; 'S103'; 'S104'; 'S107';'S109'; ...
    'S111'; 'S112'; 'S114'; 'S115'; 'S201'; 'S202'; 'S207' }; 
    
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
    EEG = []; ALLEEG = [];

    file_names=dir([path_data filesep IDs{idx} filesep IDs{idx} '*.mat']);
    if isempty(file_names)
        file_names=dir([path_data filesep IDs{idx} '*.mat']);
        if isempty(file_names)
            continue;
        end
    end
    for nF=1:length(file_names)
        %%% load data
        load([file_names(nF).folder filesep file_names(nF).name]);
        a = length(EEG_120Hz.event);
        newTypeValues = repmat(['block_', num2str(nF)], a, 1);
        for i = 1:a
            EEG_120Hz.event(i).block = newTypeValues(i,:);
        end
        % Save EEG for each block type for each participant
        [ALLEEG, EEG] = eeg_store(ALLEEG,EEG_120Hz);
        
    end
    
    % Merge all blocks for each subject
    EEG = pop_mergeset(  ALLEEG, 1:length(file_names), 0);
    
    EEG.chanlocs = rmfield(EEG.chanlocs, 'type');
    
    %MB
    EEG = pop_tesa_fastica( EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );
    
    EEG= pop_tesa_compselect( EEG,'compCheck','off','comps',[],'figSize','small','plotTimeX',[0 500],'plotFreqX',[1 100],'tmsMuscle',...
        'off','tmsMuscleThresh',8,'tmsMuscleWin',[11 30],'tmsMuscleFeedback','off','blink','on','blinkThresh',2.5,'blinkElecs',{'Fp1','Fp2'},...
        'blinkFeedback','off','move','on','moveThresh',2,'moveElecs',{'F7','F8'},'moveFeedback','off','muscle','off','muscleThresh',-0.3,...
        'muscleFreqWin',[30 100],'muscleFreqIn',[7 75],'muscleFeedback','off','elecNoise','on','elecNoiseThresh',4,'elecNoiseFeedback','off' );
    boundaries=[1 match_str({EEG.event.type},'boundary')' length({EEG.event.type})];
    all_samples=[EEG.event.latency];
    
    chan_labels={EEG_120Hz.chanlocs.labels};
    EEG1 = [];
    for nF=1:length(file_names)
        % Extract the data from each condition
        EEG1 = pop_select( EEG, 'point',[all_samples(boundaries(nF)) all_samples(boundaries(nF+1))]);
        
        
        data=EEG1.data;
        Fs=EEG1.srate;

        %start from first stimulus
        if isnumeric([EEG1.event.type])==1
            stimEvents=find([EEG1.event.type]==4);
        else
            stimEvents=match_str({EEG1.event.type},'4');
        end
        this_Sart=EEG1.event(stimEvents(1)).latency; 
        data=data(:,this_Sart:end);

        Duration=size(data,2)/Fs/60;
        FileName=file_names(nF).name;
        separators=findstr(FileName,'_');
        SubID=FileName(1:separators(1)-1);
        BlockID=str2num(FileName(separators(1)+6:separators(2)-1));
        GroupID=FileName(separators(2)+1:separators(end)-1);
        
        %%% Preprocess
        data=data-repmat(mean(data(match_str(chan_labels,{'TP7','TP8'}),:),1),size(data,1),1);
        data=data-mean(data,2);
        
        %%% Detect all slow waves
        if exist([path_save filesep 'allSW_trim_' file_names(nF).name])==0 || redo==1
            [twa_results]=twalldetectnew_TA_v2(data,Fs,0);
            all_Waves=[];
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
            save([path_save filesep 'allSW_trim_' file_names(nF).name '_newfilt'],'all_Waves','Fs','chan_labels');
        else
            load([path_save filesep 'allSW_trim_' file_names(nF).name '_newfilt']);
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
        for nE=1:size(data,1)
            thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
            temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
            
            if ~isempty(paramSW.fixThr)
                thr_Wave(nE)=paramSW.fixThr;
            else
                thr_Wave(nE)=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
            end
            slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>thr_Wave(nE),:)];
        end
        save([path_save filesep 'SW_trim_' file_names(nF).name],'slow_Waves','Fs','chan_labels');
        
        slow_Waves_perE=[];
%         allSlow_Waves_perE=[];
        for nE=1:length(chan_labels)
            slow_Waves_perE=[slow_Waves_perE ; [sum(slow_Waves(:,3)==nE)/Duration nanmean(slow_Waves(slow_Waves(:,3)==nE,4)) nanmean(1./((slow_Waves(slow_Waves(:,3)==nE,7)-slow_Waves(slow_Waves(:,3)==nE,5))/Fs)) ...
                nanmean(slow_Waves(slow_Waves(:,3)==nE,12)) nanmean(slow_Waves(slow_Waves(:,3)==nE,13))]];
            
%             allSlow_Waves_perE=[allSlow_Waves_perE ; [sum(all_Waves(:,3)==nE)/Duration nanmean(all_Waves(all_Waves(:,3)==nE,4)) nanmean(1./((all_Waves(all_Waves(:,3)==nE,7)-all_Waves(all_Waves(:,3)==nE,5))/Fs)) ...
%                 nanmean(all_Waves(all_Waves(:,3)==nE,12)) nanmean(all_Waves(all_Waves(:,3)==nE,13))]];
        end
        
        table_length=size(SW_table,1);
        SW_table.SubID(table_length+(1:length(chan_labels)))=repmat({SubID},length(chan_labels),1);
        SW_table.GroupID(table_length+(1:length(chan_labels)))=repmat({GroupID},length(chan_labels),1);
        SW_table.Block(table_length+(1:length(chan_labels)))=repmat(BlockID,length(chan_labels),1);
        SW_table.BlockDuration(table_length+(1:length(chan_labels)))=repmat(Duration,length(chan_labels),1);
        SW_table.Elec(table_length+(1:length(chan_labels)))=chan_labels;
        SW_table.SW_density(table_length+(1:length(chan_labels)))=slow_Waves_perE(:,1);
        SW_table.SW_amplitude(table_length+(1:length(chan_labels)))=slow_Waves_perE(:,2);
        SW_table.SW_frequency(table_length+(1:length(chan_labels)))=slow_Waves_perE(:,3);
        SW_table.SW_downslope(table_length+(1:length(chan_labels)))=slow_Waves_perE(:,4);
        SW_table.SW_upslope(table_length+(1:length(chan_labels)))=slow_Waves_perE(:,5);
        
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
writetable(SW_table,[path_save filesep 'SW_trim_individualThreshold_newfilt.csv'])
% writetable(allSW_table,[path_save filesep 'SW_noThreshold.csv'])