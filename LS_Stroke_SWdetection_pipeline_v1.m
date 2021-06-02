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

mkdir([path_data filesep 'SWdetection']);

%%
sw_thr=[];
for nF=1:length(eeg_files)
    %%% 1. load EEG data, events and header
    %     data=ft_read_data([eeg_files(nF).folder filesep eeg_files(nF).name]);
    event=ft_read_event([eeg_files(nF).folder filesep eeg_files(nF).name]);
    header=ft_read_headeLSr([eeg_files(nF).folder filesep eeg_files(nF).name]);
    
    filename=eeg_files(nF).name;
    if exist([path_data filesep 'SWdetection' filesep 'SW_all_' filename(1:end-4) '.mat'])==0
        fprintf('... ... detect all SW for %s %g/%g\n',filename,nF,length(eeg_files))
        %%% 2. pre-process continuous EEG data for SW detection (re-reference to
        % mastoids)
        M1_idx=find(~cellfun(@isempty,regexp(header.label,'TP9')));
        M2_idx=find(~cellfun(@isempty,regexp(header.label,'TP10')));
        
        cfg=[];
        cfg.dataset             = [eeg_files(nF).folder filesep eeg_files(nF).name];
        cfg.reref      = 'yes';
        cfg.refchannel = [M1_idx M2_idx];
        
        cfg.hpfilter       = 'yes';        % enable high-pass filtering
        cfg.hpfilttype     = 'but';
        cfg.hpfiltord         = 4;
        cfg.hpfreq         = 0.1;
        
        cfg.lpfilter       = 'yes';        % enable high-pass filtering
        cfg.lpfilttype     = 'but';
        cfg.lpfiltord         = 4;
        cfg.lpfreq         = 40;
        
        cfg.dftfilter      = 'yes';        % enable notch filtering to eliminate power line noise
        cfg.dftfreq        = [50];
        
        dataset                   = ft_preprocessing(cfg); % read raw data
        data=dataset.trial{1}; %rows: channels / columns: sample
        
        
        %%% 3. Detect all slow waves candidates
        [twa_results]=twalldetectnew_TA_v2(data,dataset.fsample,0);
        if filename(1)=='A'
            SubjID=str2num(filename(2:4));
            BlockID=str2num(filename(7));
        else
            SubjID=str2num(filename(3:5));
            BlockID=str2num(filename(7));
        end
        all_Waves=[];
        for nE=1:size(data,1)
            all_Waves=[all_Waves ; [repmat([SubjID BlockID nE],length(abs(cell2mat(twa_results.channels(nE).maxnegpkamp))),1) abs(cell2mat(twa_results.channels(nE).maxnegpkamp))'+abs(cell2mat(twa_results.channels(nE).maxpospkamp))' ...
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
        end
        all_Waves=double(all_Waves);
        labels=dataset.label;
        Fs=dataset.fsample;
        save([path_data filesep 'SWdetection' filesep 'SW_all_' filename(1:end-4)],'all_Waves','labels','Fs');
    else
        load([path_data filesep 'SWdetection' filesep 'SW_all_' filename(1:end-4)]); %,'all_Waves','labels','Fs');
    end
    %%% Filter the top amplitude slow-waves per block
    fprintf('... ... clean SW detection for %s %g/%g\n',filename,nF,length(eeg_files))
    % parameters of cleaning SW detection
    paramSW.fixThr=[]; % if you want to use a fix threshold (eg 50) leave empty ([]) if you want to use the relative
    paramSW.prticle_Thr=90; % Choose percentile that you want to select: 80 or 90 or 95
    paramSW.LimFrqW=[1 4]; % CURRENTLY NOT USED: Freq range you want to select: [1 4] or [4 10] in Hz
    paramSW.AmpCriterionIdx=4; % Criterion to select waves on: 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P, peak-to-peak amplitude)
    paramSW.art_ampl=150; % Rejection criterion (max abs amplitude)
    paramSW.max_posampl=75; % Rejection criterion (max positive amplitude)
    paramSW.max_Freq=7; % Rejection criterion (max frequency of individual wave)
    paramSW.byElec=0; % 1: compute threshold by electrode, 0 across all electrodes
    
    % clean SW detection
    if isempty(all_Waves)
        continue;
    end
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./Fs);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>paramSW.max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>paramSW.max_posampl | all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl)*100)
    all_Waves(all_freq>paramSW.max_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];
    
    slow_Waves=[];
    for nE=1:header.nChans
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
        
        if ~isempty(paramSW.fixThr)
            thr_Wave=paramSW.fixThr;
        elseif paramSW.byElec==0
            thr_Wave=prctile(all_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
        elseif paramSW.byElec==1
            thr_Wave=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
        end
        sw_thr=[sw_thr ; [all_Waves(1,1) all_Waves(1,2) nE thr_Wave]];
        slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>thr_Wave,:)];
    end
    if ~isempty(paramSW.fixThr)
        save([path_data filesep 'SWdetection'  filesep 'SW_clean_fixThr_' filename(1:end-4)],'slow_Waves','labels','Fs','paramSW');
    elseif paramSW.byElec==0
        save([path_data filesep 'SWdetection'  filesep 'SW_clean_90thr_P2P_allE_' filename(1:end-4)],'slow_Waves','labels','Fs','paramSW');
    elseif paramSW.byElec==1
        save([path_data filesep 'SWdetection'  filesep 'SW_clean_90thr_P2P_byE_' filename(1:end-4)],'slow_Waves','labels','Fs','paramSW');
    end
    
end