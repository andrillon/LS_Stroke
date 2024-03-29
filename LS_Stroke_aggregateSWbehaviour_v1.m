%%
close all
clear all

if isempty(findstr(pwd,'thandrillon'))==0
    path_data = '/Users/thandrillon/Data/StrokeDataEx/EEG';
    path_save = '/Users/thandrillon/Data/StrokeDataEx/SWdetection';
    path_LSCPtools='/Users/thandrillon/WorkGit/LSCPtools/';
elseif isempty(findstr(pwd,'Daniel'))==0
    path_data = '/fs04/so34/Daniel/Data/Raw';
    path_save = '/fs04/so34/Daniel/Data/SWdetection';
    path_LSCPtools = '/fs04/so34/LocalSleep/Stroke/Scripts';
end
addpath(genpath(path_LSCPtools));
IDs = {'HN996'; 'HN968'; 'HN969'; 'HN970'; 'HN971'; 'HN972'; 'HN973'; 'HN974'; 'HN976'; 'HN977'; ...
       'HN978'; 'HN980'; 'HN981'; 'HN982'; 'HN983'; 'HN985'; 'HN986'; 'HN987'; 'HN988'; 'HN989'; ...
       'HN990'; 'HN992'; 'HN993'; 'HN994'; 'HN995'; 'HN998'; 'HN999'; 'S002'; 'S003'; 'S004'; ...
       'S010'; 'S012'; 'S013'; 'S014'; 'S016';'S017'; 'S018'; 'S020'; 'S025'; 'S026'; ...
       'S027'; 'S029'; 'S030'; 'S031'; 'S032'; 'S033'; 'S103'; 'S104'; 'S107';'S109'; ...
       'S111'; 'S112'; 'S114'; 'S115'; 'S201'; 'S202'; 'S207' };

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
Left_Patches=[101 105 109 102 106 110];
Right_Patches=[103 107 111 104 108 112];

%%
redo=0;
behav_SW_table=array2table(zeros(0,11),'VariableNames',{'SubID','GroupID','Elec','Block','TrialCond','StimSide','Response','FixBreak','SW','SW_amplitude','SW_downslope'});
behav_SW_table.SubID=categorical(behav_SW_table.SubID);
behav_SW_table.GroupID=categorical(behav_SW_table.GroupID);
behav_SW_table.Elec=categorical(behav_SW_table.Elec);
behav_SW_table.StimSide=categorical(behav_SW_table.StimSide);

load([path_save filesep 'SW_individualThreshold_withICA4_subGroups.mat'],'SW_table');

%     allSW_table=array2table(zeros(0,9),'VariableNames',{'SubID','GroupID','Elec','Block','SW_density','SW_amplitude','SW_frequency','SW_downslope','SW_upslope'});
%     allSW_table.SubID=categorical(allSW_table.SubID);
%     allSW_table.GroupID=categorical(allSW_table.GroupID);
%     allSW_table.Elec=categorical(allSW_table.Elec);

for idx = 1:length(IDs)
    file_names=dir([path_data filesep IDs{idx} filesep IDs{idx} '*.mat']);
    if isempty(file_names)
        continue;
    end
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
        GroupID=cell2mat(SW_table.subGroupID(find(ismember(SW_table.SubID(:,1),SubID),1)));
%         GroupID=FileName(separators(2)+1:separators(end)-1);
        
        load([path_save filesep 'SW_' file_names(nF).name]); %,'slow_Waves','Fs','chan_labels');
        
        my_events=EEG_120Hz.event;
        event_types=[my_events.type];
        event_latencies=[my_events.latency];
        trial_onset_idx=find(event_types==4);
        trial_random_idx=find(event_types==5);
        trial_resp_idx=find(event_types==12);
        trial_break_idx=find(event_types==28);
        trial_coherent_idx=find(ismember(event_types,101:112));
        fprintf('... ... detected %g trials\n',length(trial_onset_idx))
        for nTr=1:length(trial_onset_idx)
            lat_onset=event_latencies(trial_onset_idx(nTr));
            if nTr<length(trial_onset_idx)
                lat_offset=event_latencies(trial_onset_idx(nTr+1));
            else
                lat_offset=event_latencies(end);
            end
            
            temp_slow_Waves=slow_Waves(slow_Waves(:,5)>lat_onset & slow_Waves(:,7)<lat_offset,:);
            
            SW_presence=zeros(size(data,1),1);
            SW_presence(temp_slow_Waves(:,3))=1;
            
            SW_amplitude=nan(size(data,1),1);
            SW_amplitude(temp_slow_Waves(:,3))=temp_slow_Waves(:,4);
            
            SW_downslope=nan(size(data,1),1);
            SW_downslope(temp_slow_Waves(:,3))=temp_slow_Waves(:,12);
            
           
            trial_coherent_cond=event_types(trial_coherent_idx(event_latencies(trial_coherent_idx)>lat_onset & event_latencies(trial_coherent_idx)<lat_offset));
            if ~isempty(trial_coherent_cond)
                if ismember(trial_coherent_cond,Left_Patches)
                    StimSide='Left';
                elseif ismember(trial_coherent_cond,Right_Patches)
                    StimSide='Right';
                end
            else
                StimSide='Undefined';
                trial_coherent_cond=NaN;
            end
            resp_idx=trial_resp_idx(event_latencies(trial_resp_idx)>lat_onset & event_latencies(trial_resp_idx)<lat_offset);
            if isempty(resp_idx)
                Response=0;
            else
                Response=1;
            end
            break_idx=trial_break_idx(event_latencies(trial_break_idx)>lat_onset & event_latencies(trial_break_idx)<lat_offset);
            if isempty(break_idx)
                FixBreak=0;
            else
                FixBreak=1;
            end
            
             table_length=size(behav_SW_table,1);
            %              {'SubID','GroupID','Elec','Block','TrialCond','StimSide','MotionDirection','Response','FixBreak','SW'}
            behav_SW_table.SubID(table_length+(1:length(chan_labels)))=repmat({SubID},length(chan_labels),1);
            behav_SW_table.GroupID(table_length+(1:length(chan_labels)))=repmat({GroupID},length(chan_labels),1);
            behav_SW_table.Block(table_length+(1:length(chan_labels)))=repmat(BlockID,length(chan_labels),1);
            behav_SW_table.Elec(table_length+(1:length(chan_labels)))=chan_labels;
            
            if length(trial_coherent_cond)> 1
                behav_SW_table.TrialCond(table_length+(1:length(chan_labels)))=repmat(trial_coherent_cond(1),length(chan_labels),1);
            else
                behav_SW_table.TrialCond(table_length+(1:length(chan_labels)))=repmat(trial_coherent_cond,length(chan_labels),1);
            end
            behav_SW_table.StimSide(table_length+(1:length(chan_labels)))=repmat({StimSide},length(chan_labels),1);
            behav_SW_table.Response(table_length+(1:length(chan_labels)))=repmat(Response,length(chan_labels),1);
            behav_SW_table.FixBreak(table_length+(1:length(chan_labels)))=repmat(FixBreak,length(chan_labels),1);
            
            behav_SW_table.SW(table_length+(1:length(chan_labels)))=SW_presence;
            behav_SW_table.SW_amplitude(table_length+(1:length(chan_labels)))=SW_amplitude;
            behav_SW_table.SW_downslope(table_length+(1:length(chan_labels)))=SW_downslope;
            
            %             temp=trial_random_idx(trial_random_idx>trial_onset_idx(nTr));
            %             lat_random=event_latencies(temp(1));
            %             temp=trial_coherent_idx(trial_coherent_idx>trial_onset_idx(nTr));
            %             lat_coherent=event_latencies(temp(1));
        end
    end
end
writetable(behav_SW_table,[path_save filesep 'SW_and_Behav_perTrial.csv'])
