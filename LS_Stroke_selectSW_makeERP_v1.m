%%
close all
clear all
if isempty(findstr(pwd,'thandrillon'))==0
    path_data = '/Users/thandrillon/Data/StrokeData/EEG';
    path_save = '/Users/thandrillon/Data/StrokeData/newSWdetection';
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
nFc=0;
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
    EEG = pop_tesa_fastica( EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );
    EEG= pop_tesa_compselect( EEG,'compCheck','off','comps',[],'figSize','small','plotTimeX',[0 500],'plotFreqX',[1 100],'tmsMuscle',...
        'off','tmsMuscleThresh',8,'tmsMuscleWin',[11 30],'tmsMuscleFeedback','off','blink','on','blinkThresh',2.5,'blinkElecs',{'Fp1','Fp2'},...
        'blinkFeedback','off','move','on','moveThresh',2,'moveElecs',{'F7','F8'},'moveFeedback','off','muscle','off','muscleThresh',-0.3,...
        'muscleFreqWin',[30 100],'muscleFreqIn',[7 75],'muscleFeedback','off','elecNoise','on','elecNoiseThresh',4,'elecNoiseFeedback','off' );
    boundaries=[1 match_str({EEG.event.type},'boundary')' length({EEG.event.type})];
    evt_lat=[EEG.event.latency];
    samples_boundaries=[1 evt_lat(match_str({EEG.event.type},'boundary')) size(EEG.data,2)];
    all_samples=[EEG.event.latency];
    chan_labels={EEG_120Hz.chanlocs.labels};
        Fs=EEG_120Hz.srate;

    Duration=[];
    BlockInfo=[];
    merged_all_Waves=[];
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
        if exist([path_save filesep 'allSW_' file_names(nF).name])==0
            fprintf('... %s - %2.0f - %s MISSING\n',SubID,BlockID,GroupID)
            continue;
        end
        load([path_save filesep 'allSW_' file_names(nF).name]);
        merged_all_Waves=[merged_all_Waves ; double(all_Waves)];
        load([file_names(nF).folder filesep file_names(nF).name]);
        chan_labels={EEG_120Hz.chanlocs.labels};
        data=EEG_120Hz.data;
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

    %%% Get ERP by block
    FileNamesSep=split(file_names(nF).name,'_');
    ERP_SW=cell(1,length(chan_labels));
    ERP_SW_avref=cell(1,length(chan_labels));
    for nB=1:size(BlockInfo,1)
        beg=samples_boundaries(BlockInfo(nB,1));
        ed=samples_boundaries(BlockInfo(nB,1)+1);
        if (beg-round(beg))~=0
            beg=floor(beg)+1;
        end
        if (ed-round(ed))~=0
            ed=floor(ed);
        end
        data=EEG.data(:,beg:ed);
        temp_data=data-repmat(mean(data(match_str(chan_labels,{'TP7','TP8'}),:),1),[length(chan_labels) 1]);
        temp_data2=data-repmat(mean(data,1),[length(chan_labels) 1]);
        for nE=1:length(chan_labels)
            temp_slow_Waves=slow_Waves(slow_Waves(:,2)==BlockInfo(nB,4) & slow_Waves(:,3)==nE,:);
            if isempty(temp_slow_Waves)
                continue;
            end

            for nW=1:size(temp_slow_Waves,1)
                if min((-0.3*Fs:Fs)+temp_slow_Waves(nW,5))>0 && max((-0.3*Fs:Fs)+temp_slow_Waves(nW,5))<length(temp_data)
                    temp=temp_data(nE,(-0.25*Fs:Fs)+temp_slow_Waves(nW,5));
                    temp=temp-mean(temp(1:0.25*Fs));

                    temp2=temp_data2(nE,(-0.25*Fs:Fs)+temp_slow_Waves(nW,5));
                    temp2=temp2-mean(temp2(1:0.25*Fs));

                    if max(abs(temp))<150
                        ERP_SW{nE}=[ERP_SW{nE} ; temp];
                         ERP_SW_avref{nE}=[ERP_SW_avref{nE} ; temp2];
                   end
                end
            end
        end
    end
    nFc=nFc+1;
    for nE=1:length(chan_labels)
        if size(ERP_SW{nE},1)<30
            mean_ERP_SW(nFc,nE,:)=nan(1,length((-0.25*Fs:Fs)));
            mean_ERP_SW_avref(nFc,nE,:)=nan(1,length((-0.25*Fs:Fs)));
        else
            mean_ERP_SW(nFc,nE,:)=nanmean(ERP_SW{nE},1);
             mean_ERP_SW_avref(nFc,nE,:)=nanmean(ERP_SW_avref{nE},1);
       end
    end
    mean_ERP_Cond{nFc}=GroupID;

end
%%
myChannels=chan_labels; %{'Fz','Cz','Pz','Oz'};
cmap=cbrewer('qual','Set2',3);
% figure('Position',[440     9   782   788]); hb=[];
for nP=1:length(myChannels)
%     subplot(2,2,nP);
figure;
    hold on;
    
    xTime=(-0.25*Fs:Fs)/Fs;
    temp_plot1=squeeze(mean_ERP_SW(match_str(mean_ERP_Cond,'healthy_old'),match_str(chan_labels,myChannels{nP}),:));
    [~, hb(1)]=simpleTplot(xTime,temp_plot1,0,cmap(1,:),0,'-',0.5,1,5,1,[]);

    temp_plot2=squeeze(mean_ERP_SW(match_str(mean_ERP_Cond,'neglect'),match_str(chan_labels,myChannels{nP}),:));
    [~, hb(2)]=simpleTplot(xTime,temp_plot2,0,cmap(2,:),0,'-',0.5,1,5,1,[]);
    xlim([-0.25 1])
    ylim([-10 5])
    format_fig;
    xlabel('time from onset (s)')
    ylabel('\muV')
    title(myChannels{nP})
%     if nP==1
        legend(hb,{'healthy','neglect'})
%     end
% data=[temp_plot1 ; temp_plot2];
% group=[zeros(size(temp_plot1,1),1) ; ones(size(temp_plot1,1),1)];
% [realpos realneg]=get_cluster_permutation_aov(data,group,0.05,0.05,1000,xTime,'full',[]);    
end

%%
% % Compute the mean voltage across all electrodes for each time point
% temp_plot=squeeze(std(mean_ERP_SW(match_str(mean_ERP_Cond,'healthy_old'),:,:),[],2));
% 
% 
% % Normalize GFP
% normalized_gfp = gfp / max(gfp);


figure; hb=[];
hold on;
xTime=(-0.25*Fs:Fs)/Fs;
temp_plot=squeeze(mean_ERP_SW_avref(match_str(mean_ERP_Cond,'healthy_old'),:,:));
mean_voltage = nanmean(temp_plot, 2);
% Compute the GFP for each time point
gfp = squeeze(sqrt(nanmean((temp_plot - repmat(mean_voltage,[1 size(temp_plot,2) 1])).^2, 2)));

[~, hb(1)]=simpleTplot(xTime,gfp,0,cmap(1,:),0,'-',0.5,1,[],1,[]);

temp_plot=squeeze(mean_ERP_SW_avref(match_str(mean_ERP_Cond,'neglect'),:,:));
mean_voltage = nanmean(temp_plot, 2);
% Compute the GFP for each time point
gfp = squeeze(sqrt(nanmean((temp_plot - repmat(mean_voltage,[1 size(temp_plot,2) 1])).^2, 2)));

[~, hb(2)]=simpleTplot(xTime,gfp,0,cmap(2,:),0,'-',0.5,1,[],1,[]);
xlim([-0.25 1])
% ylim([-10 5])
format_fig;
xlabel('time from onset (s)')
ylabel('\muV')
title('GFP')
    legend(hb,{'healthy','neglect'})
