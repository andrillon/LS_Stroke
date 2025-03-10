%%
close all
clear all

%%
if isempty(findstr(pwd,'thandrillon'))==0
    path_data='/Users/thandrillon/Data/StrokeDataEx/EEG';
    path_save='/Users/thandrillon/Data/StrokeDataEx/SWdetection';
    path_LSCPtools='/Users/thandrillon/WorkGit/LSCPtools/';
    addpath(genpath(path_LSCPtools));
elseif isempty(findstr(pwd,'yourusername'))==0
    path_data='';
    path_save='';
    path_LSCPtools='';
    addpath(genpath(path_LSCPtools));
end
file_names=dir([path_data filesep '*.mat']);

%%
redo=0;
for nF=1:length(file_names)
    %%% load data
    load([path_data filesep file_names(nF).name]);
    chan_labels={EEG_120Hz.chanlocs.labels};
    data=EEG_120Hz.data;
    Fs=EEG_120Hz.srate;
    
    %%% Preprocess
    data=data-repmat(mean(data(match_str(chan_labels,{'TP7','TP8'}),:),1),size(data,1),1);
    data=data-mean(data,2);
    
    %%% Detect all slow waves
    if exist([path_save filesep 'allSW_' file_names(nF).name])==0 || redo==1
        [twa_results]=twalldetectnew_TA_v2(data,Fs,0);
        all_Waves=[];
        for nE=1:size(data,1)
            all_Waves=[all_Waves ; [repmat([nF 0 nE],length(abs(cell2mat(twa_results.channels(nE).maxnegpkamp))),1) abs(cell2mat(twa_results.channels(nE).maxnegpkamp))'+abs(cell2mat(twa_results.channels(nE).maxpospkamp))' ...
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
        save([path_save filesep 'allSW_' file_names(nF).name],'all_Waves','Fs','chan_labels');
    else
        load([path_save filesep 'allSW_' file_names(nF).name]);
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
    
    thr_Wave=[];
    slow_Waves=[];
    ERP_SW=[];
    for nE=1:size(data,1)
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
        
        if ~isempty(paramSW.fixThr)
            thr_Wave(nE)=paramSW.fixThr;
        else
            thr_Wave(nE)=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
        end
        slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>thr_Wave(nE),:)];
        
        temp_ERP=[];
        temp_Waves=thisE_Waves(temp_p2p>thr_Wave(nE),:);
        for nW=1:size(temp_Waves,1)
            wave_onset=temp_Waves(nW,5);
            if min(wave_onset+(-1*Fs:1*Fs))<1 || max(wave_onset+(-1*Fs:1*Fs))>size(data,2)
                continue;
            end
            vec_EEG=data(nE,wave_onset+(-1*Fs:1*Fs));
            vec_EEG=vec_EEG-mean(vec_EEG(1:Fs/2));
            if max(abs(vec_EEG))>150
                continue;
            end
            temp_ERP=[temp_ERP ; vec_EEG];
        end
        if size(temp_ERP,1)>5
            ERP_SW(nE,:)=mean(temp_ERP,1);
        else
            ERP_SW(nE,:)=nan(1,length((-1*Fs:1*Fs)));
        end
    end
    save([path_save filesep 'SW_' file_names(nF).name],'slow_Waves','Fs','chan_labels');

    %%% Examine SW properties
    cfg=[];
    cfg.layout='easycapM1.mat';
    cfg.channel=chan_labels;
    layout=ft_prepare_layout(cfg);
    correspCh=[];
    for nCh=1:length(chan_labels)
        correspCh(nCh)=match_str(layout.label,chan_labels{nCh});
    end
    layout.label(1:length(chan_labels))=layout.label(correspCh);
    layout.pos(1:length(chan_labels),:)=layout.pos(correspCh,:);
    
    figure;
    num_SW_perElec=hist(slow_Waves(:,3),length(chan_labels));
    simpleTopoPlot_ft(num_SW_perElec, layout,'labels',[],0,1);
    colorbar;
    title('SW occurrence')
    
    ampl_SW_perElec=nan(1,length(chan_labels));
    for nE=1:length(chan_labels)
        ampl_SW_perElec(nE)=nanmean(slow_Waves(slow_Waves(:,3)==nE,4));
    end
    figure;
    simpleTopoPlot_ft(ampl_SW_perElec, layout,'labels',[],0,1);
    colorbar;
    title('SW amplitude')
    
    figure;
    plot((-1*Fs:1*Fs)/Fs,ERP_SW','Color','k');
    hold on;
    plot((-1*Fs:1*Fs)/Fs,nanmean(ERP_SW,1),'Color','r','LineWidth',2);
    xlim([-0.5 1])
    
end