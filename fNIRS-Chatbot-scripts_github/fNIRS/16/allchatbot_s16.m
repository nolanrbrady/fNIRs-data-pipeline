clear;
clc;

%% group number
group_number = 1;

%% load data
data_dir = '/Users/xuhan/Desktop/fNIRS_chatbot/realresults/fNIRS/16/2022-04-27_002/';
raw_first = nirs.io.loadDirectory(data_dir, {'group','subject'});

data_dir = '/Users/xuhan/Desktop/fNIRS_chatbot/realresults/fNIRS/16/2022-04-27_003/';
raw_second = nirs.io.loadDirectory(data_dir, {'group','subject'});

data_dir = '/Users/xuhan/Desktop/fNIRS_chatbot/realresults/fNIRS/16/2022-04-27_004/';
raw_third = nirs.io.loadDirectory(data_dir, {'group','subject'});

data_dir = '/Users/xuhan/Desktop/fNIRS_chatbot/realresults/fNIRS/16/2022-04-27_005/';
raw_fourth = nirs.io.loadDirectory(data_dir, {'group','subject'});

root_dir = '/Users/xuhan/Desktop/fNIRS_chatbot/realresults/fNIRS/16/chatbotsection/';

%% Create demographics table

demographics = nirs.createDemographicsTable(raw_first);
disp(demographics)

%% Manually add events info
timestamp_start_chatbot = '15:56:05';
chatbot1 = 'open_neutral';
timestamp_chatbot1_start = '16:15:58';
timestamp_chatbot1_end = '16:16:31';

chatbot2 = 'open_neutral';
timestamp_chatbot2_start = '16:16:56';
timestamp_chatbot2_end = '16:17:43';

chatbot2_survey = 'open_neutral_survey';
timestamp_chatbot2_survey_start = '16:18:06';
timestamp_chatbot2_survey_end = '16:19:55';

chatbot3 = 'choice_agreeable';
timestamp_chatbot3_start = '16:20:18';
timestamp_chatbot3_end = '16:20:55';

chatbot4 = 'choice_agreeable';
timestamp_chatbot4_start = '16:21:19';
timestamp_chatbot4_end = '16:22:21';

timestamp_firstpart_end = '16:22:18';

timestamp_start_chatbot_second = '16:26:48';

%timestamp_start_five2 = '16:27:47';
timestamp_start_ten2 = '16:27:29';

chatbot5 = 'choice_neutral';
timestamp_chatbot5_start = '16:27:52';
timestamp_chatbot5_end = '16:28:32';

chatbot6 = 'choice_neutral';
timestamp_chatbot6_start = '16:28:58';
timestamp_chatbot6_end = '16:29:37';%'16:29:52';

timestamp_secondpart_end = '16:29:37';


timestamp_start_chatbot_third = '16:32:14';

%timestamp_start_five3 = '16:34:05';
timestamp_start_ten3 = '16:33:47';

chatbot7 = 'open_agreeable';
timestamp_chatbot7_start = '16:34:10';
timestamp_chatbot7_end = '16:34:52';

chatbot8 = 'open_agreeable';
timestamp_chatbot8_start = '16:35:17';
timestamp_chatbot8_end = '16:36:23';

timestamp_thirdpart_end = '16:36:46';


timestamp_start_chatbot_fourth = '16:41:57';

%timestamp_start_five4 = '16:43:03';
timestamp_start_ten4 = '16:42:45';

chatbot9 = 'choice_agreeable';
timestamp_chatbot9_start = '16:43:08';
timestamp_chatbot9_end = '16:43:44';

chatbot10 = 'choice_agreeable';
timestamp_chatbot10_start = '16:44:09';
timestamp_chatbot10_end = '16:45:20';

chatbot10_survey = 'choice_agreeable_survey';
timestamp_chatbot10_survey_start = '16:45:43';
timestamp_chatbot10_survey_end = '16:46:53';

chatbot11 = 'choice_neutral';
timestamp_chatbot11_start = '16:47:16';
timestamp_chatbot11_end = '16:47:46';

chatbot12 = 'choice_neutral';
timestamp_chatbot12_start = '16:48:11';
timestamp_chatbot12_end = '16:48:57';

chatbot12_survey = 'choice_neutral_survey';
timestamp_chatbot12_survey_start = '16:49:20';
timestamp_chatbot12_survey_end = '16:50:48';

chatbot13 = 'open_agreeable';
timestamp_chatbot13_start = '16:51:11';
timestamp_chatbot13_end = '16:51:53';

chatbot14 = 'open_agreeable';
timestamp_chatbot14_start = '16:52:19';
timestamp_chatbot14_end = '16:53:19';

chatbot14_survey = 'open_agreeable_survey';
timestamp_chatbot14_survey_start = '16:53:42';
timestamp_chatbot14_survey_end = '16:55:03';

chatbot15 = 'open_neutral';
timestamp_chatbot15_start = '16:55:26';
timestamp_chatbot15_end = '16:56:02';

chatbot16 = 'open_neutral';
timestamp_chatbot16_start = '16:56:27';
timestamp_chatbot16_end = '16:57:25';

timestamp_fourthpart_end = '16:57:48';
%%
%keep the first part raw data until the end of chatbot4 + 23s (clipping)
timestamp_firstpart_relative_end = seconds(datetime(timestamp_firstpart_end,'InputFormat','HH:mm:ss') - datetime(timestamp_start_chatbot,'InputFormat','HH:mm:ss'));
raw_first.data = raw_first.data(1:round(timestamp_firstpart_relative_end*raw_first.Fs),:);
raw_first.time = 0:1/raw_first.Fs:(length(raw_first.data)-1)/raw_first.Fs;

second_relative_start = seconds(datetime(timestamp_start_ten2,'InputFormat','HH:mm:ss') - datetime(timestamp_start_chatbot_second,'InputFormat','HH:mm:ss'));
%timestamp_chatbot6_relative_end = seconds(datetime(timestamp_chatbot6_end,'InputFormat','HH:mm:ss') - datetime(timestamp_start_chatbot_second,'InputFormat','HH:mm:ss'));
timestamp_secondpart_relative_end = seconds(datetime(timestamp_secondpart_end,'InputFormat','HH:mm:ss') - datetime(timestamp_start_chatbot_second,'InputFormat','HH:mm:ss'));
%concat the first part with the second part
raw_second.data = raw_second.data(round(second_relative_start*raw_second.Fs):round(timestamp_secondpart_relative_end*raw_second.Fs),:);
raw_second.time = 0:1/raw_second.Fs:(length(raw_second.data)-1)/raw_second.Fs;

third_relative_start = seconds(datetime(timestamp_start_ten3,'InputFormat','HH:mm:ss') - datetime(timestamp_start_chatbot_third,'InputFormat','HH:mm:ss'));
%timestamp_chatbot8_relative_end = seconds(datetime(timestamp_chatbot8_end,'InputFormat','HH:mm:ss') - datetime(timestamp_start_chatbot_third,'InputFormat','HH:mm:ss'));
timestamp_thirdpart_relative_end = seconds(datetime(timestamp_thirdpart_end,'InputFormat','HH:mm:ss') - datetime(timestamp_start_chatbot_third,'InputFormat','HH:mm:ss'));
%concat the first part with the third part
raw_third.data = raw_third.data(round(third_relative_start*raw_third.Fs):round(timestamp_thirdpart_relative_end*raw_third.Fs),:);
raw_third.time = 0:1/raw_third.Fs:(length(raw_third.data)-1)/raw_third.Fs;

fourth_relative_start = seconds(datetime(timestamp_start_ten4,'InputFormat','HH:mm:ss') - datetime(timestamp_start_chatbot_fourth,'InputFormat','HH:mm:ss'));
%timestamp_chatbot16_relative_end = seconds(datetime(timestamp_chatbot16_end,'InputFormat','HH:mm:ss') - datetime(timestamp_start_chatbot_fourth,'InputFormat','HH:mm:ss'));
timestamp_fourthpart_relative_end = seconds(datetime(timestamp_fourthpart_end,'InputFormat','HH:mm:ss') - datetime(timestamp_start_chatbot_fourth,'InputFormat','HH:mm:ss'));
%concat the first part with the fourth part
raw_fourth.data = raw_fourth.data(round(fourth_relative_start*raw_fourth.Fs):round(timestamp_fourthpart_relative_end*raw_fourth.Fs),:);
raw_fourth.time = 0:1/raw_fourth.Fs:(length(raw_fourth.data)-1)/raw_fourth.Fs;

%% do preprocessing for raw_first and raw_second separately

%% --------------------------raw first------------------------
%% add eventinfo
%first 4 chatbots
temp=Dictionary();

chatbot_keys = string({chatbot1;chatbot2;chatbot3;chatbot4});
chatbot_start = [timestamp_chatbot1_start; timestamp_chatbot2_start; timestamp_chatbot3_start; timestamp_chatbot4_start];
chatbot_end = [timestamp_chatbot1_end; timestamp_chatbot2_end; timestamp_chatbot3_end; timestamp_chatbot4_end];
chatbot_OA_start_relative = [];
chatbot_OA_duration = [];
chatbot_OA_amp = [];
chatbot_ON_start_relative = [];
chatbot_ON_duration = [];
chatbot_ON_amp = [];
chatbot_CA_start_relative = [];
chatbot_CA_duration = [];
chatbot_CA_amp = [];
chatbot_CN_start_relative = [];
chatbot_CN_duration = [];
chatbot_CN_amp = [];

for chatbot_idx = 1:4
    chatbot_start_relative_temp = seconds(datetime(chatbot_start(chatbot_idx,:),'InputFormat','HH:mm:ss') - datetime(timestamp_start_chatbot,'InputFormat','HH:mm:ss'));
    chatbot_duration_temp = seconds(datetime(chatbot_end(chatbot_idx,:),'InputFormat','HH:mm:ss') - datetime(chatbot_start(chatbot_idx,:),'InputFormat','HH:mm:ss'));
    if strcmp(chatbot_keys(chatbot_idx,:), 'open_agreeable')
       chatbot_OA_start_relative = [chatbot_OA_start_relative chatbot_start_relative_temp];
       chatbot_OA_duration = [chatbot_OA_duration chatbot_duration_temp];
       chatbot_OA_amp = [chatbot_OA_amp 1];
       temp_se = nirs.design.StimulusEvents( 'open_agreeable', chatbot_OA_start_relative, chatbot_OA_duration, chatbot_OA_amp);
       temp('open_agreeable') = temp_se;
    end
    
    if strcmp(chatbot_keys(chatbot_idx,:), 'open_neutral')
       chatbot_ON_start_relative = [chatbot_ON_start_relative chatbot_start_relative_temp];
       chatbot_ON_duration = [chatbot_ON_duration chatbot_duration_temp];
       chatbot_ON_amp = [chatbot_ON_amp 1];
       temp_se = nirs.design.StimulusEvents( 'open_neutral', chatbot_ON_start_relative, chatbot_ON_duration, chatbot_ON_amp);
       temp('open_neutral') = temp_se;
    end
    
    if strcmp(chatbot_keys(chatbot_idx,:), 'choice_agreeable')
       chatbot_CA_start_relative = [chatbot_CA_start_relative chatbot_start_relative_temp];
       chatbot_CA_duration = [chatbot_CA_duration chatbot_duration_temp];
       chatbot_CA_amp = [chatbot_CA_amp 1];
       temp_se = nirs.design.StimulusEvents( 'choice_agreeable', chatbot_CA_start_relative, chatbot_CA_duration, chatbot_CA_amp);
       temp('choice_agreeable') = temp_se;
    end
    
    if strcmp(chatbot_keys(chatbot_idx,:), 'choice_neutral')
       chatbot_CN_start_relative = [chatbot_CN_start_relative chatbot_start_relative_temp];
       chatbot_CN_duration = [chatbot_CN_duration chatbot_duration_temp];
       chatbot_CN_amp = [chatbot_CN_amp 1];
       temp_se = nirs.design.StimulusEvents( 'choice_neutral', chatbot_CN_start_relative, chatbot_CN_duration, chatbot_CN_amp);
       temp('choice_neutral') = temp_se;
    end
end

raw_first(1,1).stimulus = temp;
%%
%% Quality checking 
mkdir([strcat(root_dir, 'rawdata_first')])
mkdir([strcat(root_dir, 'qualitycheck_first')])

for s = 1:size(raw_first)
% view each subject's graph
    data = raw_first(s);
    data.draw;
    saveas(gcf,[root_dir,'/rawdata_first/sub',num2str(s),'.png'])
    close;
    
% turncate for quality checking
    pre_Baseline = 23;
    post_Baseline = 23;
    fs = data.Fs;
    [m,n] = size(data.stimulus.keys);
   
    onset_firststi = timestamp_chatbot1_start;
    onset_laststi = timestamp_chatbot8_end;
    if onset_firststi < pre_Baseline
        temp_start = 1;
    else
        temp_start = onset_firststi-pre_Baseline;
    end
        
    [len, chan_num] = size(data.data);
    % if onset_laststi + post_Baseline is larger than the length of the data (the end of
    % the data stream), will cut the data stream before the end of the data
    if (onset_laststi + post_Baseline)*fs > len
        temp_end = len/fs;
    else
        temp_end = onset_laststi + post_Baseline;
    end 
    
    %disp(temp_start)
    %disp(temp_end)
    
    data = raw_first(s).data(round(temp_start*fs):round(temp_end*fs),:);

% quality checking using cv, cv defined as dev/m
    m = mean(data);
    dev = std(data);
    cvThresh = 0.15;  % can be changed
    cvMask = ones(size(m));
    cvOutput = ones(size(m));
    for i = 1:40
        cv = dev(i)/m(i);
        %disp(cv)
        if cv < cvThresh
            cvMask(i) = 1;
        else
            cvMask(i) = 0;
            disp('bad link:');
            disp(raw_first(s).probe.link(i,:));
        end
    end
    T = table(table2array(raw_first(s).probe.link(:,1)), table2array(raw_first(s).probe.link(:,2)),table2array(raw_first(s).probe.link(:,3)), transpose(cvMask));
    xlswrite([root_dir,'qualitycheck_first/sub',num2str(s),'.csv'], table2array(T));

end

%% Prepare the data for preprocessing
%% Remove stimless files
mkdir([strcat(root_dir, 'stimremove_turncate_first')])
j = nirs.modules.RemoveStimless();

j = nirs.modules.KeepStims( j);
j.listOfStims = {'choice_agreeable', 'open_neutral', 'open_agreeable', 'choice_neutral'};


% resample to 4Hz
j = nirs.modules.Resample( j );
j.Fs = 4;  % Sets the new sample rate to 4Hz (was 10Hz).

% Trim pre and post baseline
j = nirs.modules.TrimBaseline( j);
j.preBaseline = 23; %can change these values
j.postBaseline = 23;
raw_first = j.run(raw_first);

for s = 1:size(raw_first)
% view each subject's graph
    raw_first(s).draw;
    saveas(gcf,[root_dir,'stimremove_turncate_first/sub',num2str(s),'.png'])
    close;
end


%%
hb_first = raw_first;
%% %% Signal Preprocessing 
mkdir([strcat(root_dir, 'preprocessing_first')])
% visualize results
for s = 1:size(hb_first)
% view each subject's graph
    hb_first(s).draw;
    saveas(gcf,[root_dir,'preprocessing_first/afterBP_sub_',num2str(s),'.png'])
    close;
end

%% SS label
job_SS = nirs.modules.LabelShortSeperation();
hb_first = job_SS.run(hb_first);
%% Run the conversions 
mkdir([strcat(root_dir, 'conversion_first')])
j = nirs.modules.OpticalDensity();
% the DFP parameters can be modified in nirs.modules.BeerLamberLaw.m

%% %% Convert to hemoglobin
j = nirs.modules.BeerLambertLaw( j);    
hb_first = j.run(hb_first);
% Results visualization
for s = 1:size(hb_first)
% view each subject's graph
    hb_first(s).draw;
    saveas(gcf,[root_dir,'conversion_first/sub',num2str(s),'.png'])
    close;
end

%%  -----------------------------raw second---------------------------
%%
%chatbot 5-6
temp=Dictionary();

chatbot_keys = string({chatbot5;chatbot6});
chatbot_start = [timestamp_chatbot5_start; timestamp_chatbot6_start];
chatbot_end = [timestamp_chatbot5_end; timestamp_chatbot6_end];
chatbot_OA_start_relative = [];
chatbot_OA_duration = [];
chatbot_OA_amp = [];
chatbot_ON_start_relative = [];
chatbot_ON_duration = [];
chatbot_ON_amp = [];
chatbot_CA_start_relative = [];
chatbot_CA_duration = [];
chatbot_CA_amp = [];
chatbot_CN_start_relative = [];
chatbot_CN_duration = [];
chatbot_CN_amp = [];

for chatbot_idx = 1:2   
    chatbot_start_relative_temp = seconds(datetime(chatbot_start(chatbot_idx,:),'InputFormat','HH:mm:ss') - datetime(timestamp_start_ten2,'InputFormat','HH:mm:ss'));
    chatbot_duration_temp = seconds(datetime(chatbot_end(chatbot_idx,:),'InputFormat','HH:mm:ss') - datetime(chatbot_start(chatbot_idx,:),'InputFormat','HH:mm:ss'));
    if strcmp(chatbot_keys(chatbot_idx,:), 'open_agreeable')
       chatbot_OA_start_relative = [chatbot_OA_start_relative chatbot_start_relative_temp];
       chatbot_OA_duration = [chatbot_OA_duration chatbot_duration_temp];
       chatbot_OA_amp = [chatbot_OA_amp 1];
       temp_se = nirs.design.StimulusEvents( 'open_agreeable', chatbot_OA_start_relative, chatbot_OA_duration, chatbot_OA_amp);
       temp('open_agreeable') = temp_se;
    end
    
    if strcmp(chatbot_keys(chatbot_idx,:), 'open_neutral')
       chatbot_ON_start_relative = [chatbot_ON_start_relative chatbot_start_relative_temp];
       chatbot_ON_duration = [chatbot_ON_duration chatbot_duration_temp];
       chatbot_ON_amp = [chatbot_ON_amp 1];
       temp_se = nirs.design.StimulusEvents( 'open_neutral', chatbot_ON_start_relative, chatbot_ON_duration, chatbot_ON_amp);
       temp('open_neutral') = temp_se;
    end
    
    if strcmp(chatbot_keys(chatbot_idx,:), 'choice_agreeable')
       chatbot_CA_start_relative = [chatbot_CA_start_relative chatbot_start_relative_temp];
       chatbot_CA_duration = [chatbot_CA_duration chatbot_duration_temp];
       chatbot_CA_amp = [chatbot_CA_amp 1];
       temp_se = nirs.design.StimulusEvents( 'choice_agreeable', chatbot_CA_start_relative, chatbot_CA_duration, chatbot_CA_amp);
       temp('choice_agreeable') = temp_se;
    end
    
    if strcmp(chatbot_keys(chatbot_idx,:), 'choice_neutral')
       chatbot_CN_start_relative = [chatbot_CN_start_relative chatbot_start_relative_temp];
       chatbot_CN_duration = [chatbot_CN_duration chatbot_duration_temp];
       chatbot_CN_amp = [chatbot_CN_amp 1];
       temp_se = nirs.design.StimulusEvents( 'choice_neutral', chatbot_CN_start_relative, chatbot_CN_duration, chatbot_CN_amp);
       temp('choice_neutral') = temp_se;
    end
end

raw_second(1,1).stimulus = temp;

%% Quality checking 
mkdir([strcat(root_dir, 'rawdata_second')])
mkdir([strcat(root_dir, 'qualitycheck_second')])

for s = 1:size(raw_second)
% view each subject's graph
    data = raw_second(s);
    data.draw;
    saveas(gcf,[root_dir,'/rawdata_second/sub',num2str(s),'.png'])
    close;
    
% turncate for quality checking
    pre_Baseline = 23;
    post_Baseline = 23;
    fs = data.Fs;
    [m,n] = size(data.stimulus.keys);
   
    onset_firststi = timestamp_chatbot8_start;
    onset_laststi = timestamp_chatbot16_end;
    if onset_firststi < pre_Baseline
        temp_start = 1;
    else
        temp_start = onset_firststi-pre_Baseline;
    end
        
    [len, chan_num] = size(data.data);
    % if onset_laststi + post_Baseline is larger than the length of the data (the end of
    % the data stream), will cut the data stream before the end of the data
    if (onset_laststi + post_Baseline)*fs > len
        temp_end = len/fs;
    else
        temp_end = onset_laststi + post_Baseline;
    end 
    
    %disp(temp_start)
    %disp(temp_end)
    
    data = raw_second(s).data(round(temp_start*fs):round(temp_end*fs),:);

% quality checking using cv, cv defined as dev/m
    m = mean(data);
    dev = std(data);
    cvThresh = 0.15;  % can be changed
    cvMask = ones(size(m));
    cvOutput = ones(size(m));
    for i = 1:40
        cv = dev(i)/m(i);
        %disp(cv)
        if cv < cvThresh
            cvMask(i) = 1;
        else
            cvMask(i) = 0;
            disp('bad link:');
            disp(raw_second(s).probe.link(i,:));
        end
    end
    T = table(table2array(raw_second(s).probe.link(:,1)), table2array(raw_second(s).probe.link(:,2)),table2array(raw_second(s).probe.link(:,3)), transpose(cvMask));
    xlswrite([root_dir,'qualitycheck_second/sub',num2str(s),'.csv'], table2array(T));

end

%% Prepare the data for preprocessing
%% Remove stimless files
mkdir([strcat(root_dir, 'stimremove_turncate_second')])
j = nirs.modules.RemoveStimless();

j = nirs.modules.KeepStims( j);
j.listOfStims = {'choice_agreeable', 'open_neutral', 'open_agreeable', 'choice_neutral'};


% resample to 4Hz
j = nirs.modules.Resample( j );
j.Fs = 4;  % Sets the new sample rate to 4Hz (was 10Hz).

% Trim pre and post baseline
j = nirs.modules.TrimBaseline( j);
j.preBaseline = 23; %can change these values
j.postBaseline = 23;
raw_second = j.run(raw_second);

for s = 1:size(raw_second)
% view each subject's graph
    raw_second(s).draw;
    saveas(gcf,[root_dir,'stimremove_turncate_second/sub',num2str(s),'.png'])
    close;
end


%%
hb_second = raw_second;
%% %% Signal Preprocessing 
mkdir([strcat(root_dir, 'preprocessing_second')])
% visualize results
for s = 1:size(hb_second)
% view each subject's graph
    hb_second(s).draw;
    saveas(gcf,[root_dir,'preprocessing_second/afterBP_sub_',num2str(s),'.png'])
    close;
end

%% SS label
job_SS = nirs.modules.LabelShortSeperation();
hb_second = job_SS.run(hb_second);
%% Run the conversions 
mkdir([strcat(root_dir, 'conversion_second')])
j = nirs.modules.OpticalDensity();
% the DFP parameters can be modified in nirs.modules.BeerLamberLaw.m

%% %% Convert to hemoglobin
j = nirs.modules.BeerLambertLaw( j);    
hb_second = j.run(hb_second);
% Results visualization
for s = 1:size(hb_second)
% view each subject's graph
    hb_second(s).draw;
    saveas(gcf,[root_dir,'conversion_second/sub',num2str(s),'.png'])
    close;
end

%% ----------------------------- raw third -----------------------------
%%
%chatbot 7-8
temp=Dictionary();

chatbot_keys = string({chatbot7;chatbot8});
chatbot_start = [timestamp_chatbot7_start; timestamp_chatbot8_start];
chatbot_end = [timestamp_chatbot7_end; timestamp_chatbot8_end];
chatbot_OA_start_relative = [];
chatbot_OA_duration = [];
chatbot_OA_amp = [];
chatbot_ON_start_relative = [];
chatbot_ON_duration = [];
chatbot_ON_amp = [];
chatbot_CA_start_relative = [];
chatbot_CA_duration = [];
chatbot_CA_amp = [];
chatbot_CN_start_relative = [];
chatbot_CN_duration = [];
chatbot_CN_amp = [];

for chatbot_idx = 1:2
    chatbot_start_relative_temp = seconds(datetime(chatbot_start(chatbot_idx,:),'InputFormat','HH:mm:ss') - datetime(timestamp_start_ten3,'InputFormat','HH:mm:ss'));
    chatbot_duration_temp = seconds(datetime(chatbot_end(chatbot_idx,:),'InputFormat','HH:mm:ss') - datetime(chatbot_start(chatbot_idx,:),'InputFormat','HH:mm:ss'));
    if strcmp(chatbot_keys(chatbot_idx,:), 'open_agreeable')
       chatbot_OA_start_relative = [chatbot_OA_start_relative chatbot_start_relative_temp];
       chatbot_OA_duration = [chatbot_OA_duration chatbot_duration_temp];
       chatbot_OA_amp = [chatbot_OA_amp 1];
       temp_se = nirs.design.StimulusEvents( 'open_agreeable', chatbot_OA_start_relative, chatbot_OA_duration, chatbot_OA_amp);
       temp('open_agreeable') = temp_se;
    end
    
    if strcmp(chatbot_keys(chatbot_idx,:), 'open_neutral')
       chatbot_ON_start_relative = [chatbot_ON_start_relative chatbot_start_relative_temp];
       chatbot_ON_duration = [chatbot_ON_duration chatbot_duration_temp];
       chatbot_ON_amp = [chatbot_ON_amp 1];
       temp_se = nirs.design.StimulusEvents( 'open_neutral', chatbot_ON_start_relative, chatbot_ON_duration, chatbot_ON_amp);
       temp('open_neutral') = temp_se;
    end
    
    if strcmp(chatbot_keys(chatbot_idx,:), 'choice_agreeable')
       chatbot_CA_start_relative = [chatbot_CA_start_relative chatbot_start_relative_temp];
       chatbot_CA_duration = [chatbot_CA_duration chatbot_duration_temp];
       chatbot_CA_amp = [chatbot_CA_amp 1];
       temp_se = nirs.design.StimulusEvents( 'choice_agreeable', chatbot_CA_start_relative, chatbot_CA_duration, chatbot_CA_amp);
       temp('choice_agreeable') = temp_se;
    end
    
    if strcmp(chatbot_keys(chatbot_idx,:), 'choice_neutral')
       chatbot_CN_start_relative = [chatbot_CN_start_relative chatbot_start_relative_temp];
       chatbot_CN_duration = [chatbot_CN_duration chatbot_duration_temp];
       chatbot_CN_amp = [chatbot_CN_amp 1];
       temp_se = nirs.design.StimulusEvents( 'choice_neutral', chatbot_CN_start_relative, chatbot_CN_duration, chatbot_CN_amp);
       temp('choice_neutral') = temp_se;
    end
end
raw_third(1,1).stimulus = temp;

%%
%% Quality checking 
mkdir([strcat(root_dir, 'rawdata_third')])
mkdir([strcat(root_dir, 'qualitycheck_third')])

for s = 1:size(raw_third)
% view each subject's graph
    data = raw_third(s);
    data.draw;
    saveas(gcf,[root_dir,'/rawdata_third/sub',num2str(s),'.png'])
    close;
    
% turncate for quality checking
    pre_Baseline = 23;
    post_Baseline = 23;
    fs = data.Fs;
    [m,n] = size(data.stimulus.keys);
   
    onset_firststi = timestamp_chatbot8_start;
    onset_laststi = timestamp_chatbot16_end;
    if onset_firststi < pre_Baseline
        temp_start = 1;
    else
        temp_start = onset_firststi-pre_Baseline;
    end
        
    [len, chan_num] = size(data.data);
    % if onset_laststi + post_Baseline is larger than the length of the data (the end of
    % the data stream), will cut the data stream before the end of the data
    if (onset_laststi + post_Baseline)*fs > len
        temp_end = len/fs;
    else
        temp_end = onset_laststi + post_Baseline;
    end 
    
    %disp(temp_start)
    %disp(temp_end)
    
    data = raw_third(s).data(round(temp_start*fs):round(temp_end*fs),:);

% quality checking using cv, cv defined as dev/m
    m = mean(data);
    dev = std(data);
    cvThresh = 0.15;  % can be changed
    cvMask = ones(size(m));
    cvOutput = ones(size(m));
    for i = 1:40
        cv = dev(i)/m(i);
        %disp(cv)
        if cv < cvThresh
            cvMask(i) = 1;
        else
            cvMask(i) = 0;
            disp('bad link:');
            disp(raw_third(s).probe.link(i,:));
        end
    end
    T = table(table2array(raw_third(s).probe.link(:,1)), table2array(raw_third(s).probe.link(:,2)),table2array(raw_third(s).probe.link(:,3)), transpose(cvMask));
    xlswrite([root_dir,'qualitycheck_third/sub',num2str(s),'.csv'], table2array(T));

end

%% Prepare the data for preprocessing
%% Remove stimless files
mkdir([strcat(root_dir, 'stimremove_turncate_third')])
j = nirs.modules.RemoveStimless();

j = nirs.modules.KeepStims( j);
j.listOfStims = {'choice_agreeable', 'open_neutral', 'open_agreeable', 'choice_neutral'};


% resample to 4Hz
j = nirs.modules.Resample( j );
j.Fs = 4;  % Sets the new sample rate to 4Hz (was 10Hz).

% Trim pre and post baseline
j = nirs.modules.TrimBaseline( j);
j.preBaseline = 23; %can change these values
j.postBaseline = 23;
raw_third = j.run(raw_third);

for s = 1:size(raw_third)
% view each subject's graph
    raw_third(s).draw;
    saveas(gcf,[root_dir,'stimremove_turncate_third/sub',num2str(s),'.png'])
    close;
end


%%
hb_third = raw_third;
%% %% Signal Preprocessing 
mkdir([strcat(root_dir, 'preprocessing_third')])
% visualize results
for s = 1:size(hb_third)
% view each subject's graph
    hb_third(s).draw;
    saveas(gcf,[root_dir,'preprocessing_third/afterBP_sub_',num2str(s),'.png'])
    close;
end

%% SS label
job_SS = nirs.modules.LabelShortSeperation();
hb_third = job_SS.run(hb_third);
%% Run the conversions 
mkdir([strcat(root_dir, 'conversion_third')])
j = nirs.modules.OpticalDensity();
% the DFP parameters can be modified in nirs.modules.BeerLamberLaw.m

%% %% Convert to hemoglobin
j = nirs.modules.BeerLambertLaw( j);    
hb_third = j.run(hb_third);
% Results visualization
for s = 1:size(hb_third)
% view each subject's graph
    hb_third(s).draw;
    saveas(gcf,[root_dir,'conversion_third/sub',num2str(s),'.png'])
    close;
end

%% ----------------------------------- raw fourth ---------------------------------
%%

%last 8 chatbots
temp=Dictionary();
chatbot_keys = string({chatbot9;chatbot10;chatbot11;chatbot12;chatbot13;chatbot14;chatbot15;chatbot16;chatbot10_survey;chatbot12_survey;chatbot14_survey});
chatbot_start = [timestamp_chatbot9_start; timestamp_chatbot10_start; timestamp_chatbot11_start; timestamp_chatbot12_start; timestamp_chatbot13_start; timestamp_chatbot14_start; timestamp_chatbot15_start; timestamp_chatbot16_start; timestamp_chatbot10_survey_start; timestamp_chatbot12_survey_start; timestamp_chatbot14_survey_start];
chatbot_end = [timestamp_chatbot9_end; timestamp_chatbot10_end; timestamp_chatbot11_end; timestamp_chatbot12_end; timestamp_chatbot13_end; timestamp_chatbot14_end; timestamp_chatbot15_end; timestamp_chatbot16_end; timestamp_chatbot10_survey_end; timestamp_chatbot12_survey_end; timestamp_chatbot14_survey_end];
chatbot_OA_start_relative = [];
chatbot_OA_duration = [];
chatbot_OA_amp = [];
chatbot_ON_start_relative = [];
chatbot_ON_duration = [];
chatbot_ON_amp = [];
chatbot_CA_start_relative = [];
chatbot_CA_duration = [];
chatbot_CA_amp = [];
chatbot_CN_start_relative = [];
chatbot_CN_duration = [];
chatbot_CN_amp = [];
chatbot_OA_survey_start_relative = [];
chatbot_OA_survey_duration = [];
chatbot_OA_survey_amp = [];
chatbot_ON_survey_start_relative = [];
chatbot_ON_survey_duration = [];
chatbot_ON_survey_amp = [];
chatbot_CA_survey_start_relative = [];
chatbot_CA_survey_duration = [];
chatbot_CA_survey_amp = [];
chatbot_CN_survey_start_relative = [];
chatbot_CN_survey_duration = [];
chatbot_CN_survey_amp = [];


for chatbot_idx = 1:11
    chatbot_start_relative_temp = seconds(datetime(chatbot_start(chatbot_idx,:),'InputFormat','HH:mm:ss') - datetime(timestamp_start_ten4,'InputFormat','HH:mm:ss'));
    chatbot_duration_temp = seconds(datetime(chatbot_end(chatbot_idx,:),'InputFormat','HH:mm:ss') - datetime(chatbot_start(chatbot_idx,:),'InputFormat','HH:mm:ss'));
    if strcmp(chatbot_keys(chatbot_idx,:), 'open_agreeable')
       chatbot_OA_start_relative = [chatbot_OA_start_relative chatbot_start_relative_temp];
       chatbot_OA_duration = [chatbot_OA_duration chatbot_duration_temp];
       chatbot_OA_amp = [chatbot_OA_amp 1];
       temp_se = nirs.design.StimulusEvents( 'open_agreeable', chatbot_OA_start_relative, chatbot_OA_duration, chatbot_OA_amp);
       temp('open_agreeable') = temp_se;
    end
    
    if strcmp(chatbot_keys(chatbot_idx,:), 'open_neutral')
       chatbot_ON_start_relative = [chatbot_ON_start_relative chatbot_start_relative_temp];
       chatbot_ON_duration = [chatbot_ON_duration chatbot_duration_temp];
       chatbot_ON_amp = [chatbot_ON_amp 1];
       temp_se = nirs.design.StimulusEvents( 'open_neutral', chatbot_ON_start_relative, chatbot_ON_duration, chatbot_ON_amp);
       temp('open_neutral') = temp_se;
    end
    
    if strcmp(chatbot_keys(chatbot_idx,:), 'choice_agreeable')
       chatbot_CA_start_relative = [chatbot_CA_start_relative chatbot_start_relative_temp];
       chatbot_CA_duration = [chatbot_CA_duration chatbot_duration_temp];
       chatbot_CA_amp = [chatbot_CA_amp 1];
       temp_se = nirs.design.StimulusEvents( 'choice_agreeable', chatbot_CA_start_relative, chatbot_CA_duration, chatbot_CA_amp);
       temp('choice_agreeable') = temp_se;
    end
    
    if strcmp(chatbot_keys(chatbot_idx,:), 'choice_neutral')
       chatbot_CN_start_relative = [chatbot_CN_start_relative chatbot_start_relative_temp];
       chatbot_CN_duration = [chatbot_CN_duration chatbot_duration_temp];
       chatbot_CN_amp = [chatbot_CN_amp 1];
       temp_se = nirs.design.StimulusEvents( 'choice_neutral', chatbot_CN_start_relative, chatbot_CN_duration, chatbot_CN_amp);
       temp('choice_neutral') = temp_se;
    end
    
    if strcmp(chatbot_keys(chatbot_idx,:), 'open_agreeable_survey')
       chatbot_OA_survey_start_relative = [chatbot_OA_survey_start_relative chatbot_start_relative_temp];
       chatbot_OA_survey_duration = [chatbot_OA_survey_duration chatbot_duration_temp];
       chatbot_OA_survey_amp = [chatbot_OA_survey_amp 1];
       temp_se = nirs.design.StimulusEvents( 'open_agreeable_survey', chatbot_OA_survey_start_relative, chatbot_OA_survey_duration, chatbot_OA_survey_amp);
       temp('open_agreeable_survey') = temp_se;
    end
    
    if strcmp(chatbot_keys(chatbot_idx,:), 'open_neutral_survey')
       chatbot_ON_survey_start_relative = [chatbot_ON_survey_start_relative chatbot_start_relative_temp];
       chatbot_ON_survey_duration = [chatbot_ON_survey_duration chatbot_duration_temp];
       chatbot_ON_survey_amp = [chatbot_ON_survey_amp 1];
       temp_se = nirs.design.StimulusEvents( 'open_neutral_survey', chatbot_ON_survey_start_relative, chatbot_ON_survey_duration, chatbot_ON_survey_amp);
       temp('open_neutral_survey') = temp_se;
    end
    
    if strcmp(chatbot_keys(chatbot_idx,:), 'choice_agreeable_survey')
       chatbot_CA_survey_start_relative = [chatbot_CA_survey_start_relative chatbot_start_relative_temp];
       chatbot_CA_survey_duration = [chatbot_CA_survey_duration chatbot_duration_temp];
       chatbot_CA_survey_amp = [chatbot_CA_survey_amp 1];
       temp_se = nirs.design.StimulusEvents( 'choice_agreeable_survey', chatbot_CA_survey_start_relative, chatbot_CA_survey_duration, chatbot_CA_survey_amp);
       temp('choice_agreeable_survey') = temp_se;
    end
    
    if strcmp(chatbot_keys(chatbot_idx,:), 'choice_neutral_survey')
       chatbot_CN_survey_start_relative = [chatbot_CN_survey_start_relative chatbot_start_relative_temp];
       chatbot_CN_survey_duration = [chatbot_CN_survey_duration chatbot_duration_temp];
       chatbot_CN_survey_amp = [chatbot_CN_survey_amp 1];
       temp_se = nirs.design.StimulusEvents( 'choice_neutral_survey', chatbot_CN_survey_start_relative, chatbot_CN_survey_duration, chatbot_CN_survey_amp);
       temp('choice_neutral_survey') = temp_se;
    end
end

raw_fourth(1,1).stimulus = temp;
%% Quality checking 
mkdir([strcat(root_dir, 'rawdata_fourth')])
mkdir([strcat(root_dir, 'qualitycheck_fourth')])

for s = 1:size(raw_fourth)
% view each subject's graph
    data = raw_fourth(s);
    data.draw;
    saveas(gcf,[root_dir,'/rawdata_fourth/sub',num2str(s),'.png'])
    close;
    
% turncate for quality checking
    pre_Baseline = 23;
    post_Baseline = 23;
    fs = data.Fs;
    [m,n] = size(data.stimulus.keys);
   
    onset_firststi = timestamp_chatbot8_start;
    onset_laststi = timestamp_chatbot16_end;
    if onset_firststi < pre_Baseline
        temp_start = 1;
    else
        temp_start = onset_firststi-pre_Baseline;
    end
        
    [len, chan_num] = size(data.data);
    % if onset_laststi + post_Baseline is larger than the length of the data (the end of
    % the data stream), will cut the data stream before the end of the data
    if (onset_laststi + post_Baseline)*fs > len
        temp_end = len/fs;
    else
        temp_end = onset_laststi + post_Baseline;
    end 
    
    %disp(temp_start)
    %disp(temp_end)
    
    data = raw_fourth(s).data(round(temp_start*fs):round(temp_end*fs),:);

% quality checking using cv, cv defined as dev/m
    m = mean(data);
    dev = std(data);
    cvThresh = 0.15;  % can be changed
    cvMask = ones(size(m));
    cvOutput = ones(size(m));
    for i = 1:40
        cv = dev(i)/m(i);
        %disp(cv)
        if cv < cvThresh
            cvMask(i) = 1;
        else
            cvMask(i) = 0;
            disp('bad link:');
            disp(raw_fourth(s).probe.link(i,:));
        end
    end
    T = table(table2array(raw_fourth(s).probe.link(:,1)), table2array(raw_fourth(s).probe.link(:,2)),table2array(raw_fourth(s).probe.link(:,3)), transpose(cvMask));
    xlswrite([root_dir,'qualitycheck_fourth/sub',num2str(s),'.csv'], table2array(T));

end

%% Prepare the data for preprocessing
%% Remove stimless files
mkdir([strcat(root_dir, 'stimremove_turncate_fourth')])
j = nirs.modules.RemoveStimless();

j = nirs.modules.KeepStims( j);
j.listOfStims = {'choice_agreeable', 'open_neutral', 'open_agreeable', 'choice_neutral'};


% resample to 4Hz
j = nirs.modules.Resample( j );
j.Fs = 4;  % Sets the new sample rate to 4Hz (was 10Hz).

% Trim pre and post baseline
j = nirs.modules.TrimBaseline( j);
j.preBaseline = 23; %can change these values
j.postBaseline = 23;
raw_fourth = j.run(raw_fourth);

for s = 1:size(raw_fourth)
% view each subject's graph
    raw_fourth(s).draw;
    saveas(gcf,[root_dir,'stimremove_turncate_fourth/sub',num2str(s),'.png'])
    close;
end


%%
hb_fourth = raw_fourth;
%% %% Signal Preprocessing 
% visualize results
mkdir([strcat(root_dir, 'preprocessing_fourth')])
for s = 1:size(hb_fourth)
% view each subject's graph
    hb_fourth(s).draw;
    saveas(gcf,[root_dir,'preprocessing_fourth/afterBP_sub_',num2str(s),'.png'])
    close;
end

%% SS label
job_SS = nirs.modules.LabelShortSeperation();
hb_fourth = job_SS.run(hb_fourth);
%% Run the conversions 
mkdir([strcat(root_dir, 'conversion_fourth')])
j = nirs.modules.OpticalDensity();
% the DFP parameters can be modified in nirs.modules.BeerLamberLaw.m

%% %% Convert to hemoglobin
j = nirs.modules.BeerLambertLaw( j);    
hb_fourth = j.run(hb_fourth);
% Results visualization
for s = 1:size(hb_fourth)
% view each subject's graph
    hb_fourth(s).draw;
    saveas(gcf,[root_dir,'conversion_fourth/sub',num2str(s),'.png'])
    close;
end

%% normalize and concate stream 1, 2, 3, and 4
%normalization
mkdir([strcat(root_dir, 'concatedata_afternormalization')])
hb_first.data = normalize(hb_first.data);
hb_second.data = normalize(hb_second.data);
hb_third.data = normalize(hb_third.data);
hb_fourth.data = normalize(hb_fourth.data);

%%
relative_nback_duration = hb_first.time(1);
hb_first.time = hb_first.time - relative_nback_duration;
relative_firstsession_duration = hb_first.time(length(hb_first.time));
hb_second.time = hb_second.time + relative_firstsession_duration;
relative_firstaddsecond_session_duration = hb_second.time(length(hb_second.time));
hb_third.time = hb_third.time + relative_firstaddsecond_session_duration;
relative_firstsecondandthird_session_duration = hb_third.time(length(hb_third.time));
hb_fourth.time = hb_fourth.time + relative_firstsecondandthird_session_duration;

temp_merge = Dictionary();

chatbot_OA_start_relative = [];
chatbot_OA_duration = [];
chatbot_OA_amp = [];
chatbot_ON_start_relative = [];
chatbot_ON_duration = [];
chatbot_ON_amp = [];
chatbot_CA_start_relative = [];
chatbot_CA_duration = [];
chatbot_CA_amp = [];
chatbot_CN_start_relative = [];
chatbot_CN_duration = [];
chatbot_CN_amp = [];

for idx_se = 1:length(hb_first.stimulus.keys)
    stimulus_key_temp = hb_first.stimulus.values{1,idx_se}.name;
    stimulus_onset_temp = hb_first.stimulus.values{1,idx_se}.onset - relative_nback_duration;
    stimulus_dur_temp = hb_first.stimulus.values{1,idx_se}.dur;
    stimulus_amp_temp = hb_first.stimulus.values{1,idx_se}.amp;
    if strcmp(stimulus_key_temp, 'open_agreeable')
       chatbot_OA_start_relative = [chatbot_OA_start_relative stimulus_onset_temp];
       chatbot_OA_duration = [chatbot_OA_duration stimulus_dur_temp];
       chatbot_OA_amp = [chatbot_OA_amp zeros(1, length(stimulus_amp_temp)) + 1];
    end
    
    if strcmp(stimulus_key_temp, 'open_neutral')
       chatbot_ON_start_relative = [chatbot_ON_start_relative stimulus_onset_temp];
       chatbot_ON_duration = [chatbot_ON_duration stimulus_dur_temp];
       chatbot_ON_amp = [chatbot_ON_amp zeros(1, length(stimulus_amp_temp)) + 1];
    end
    
    if strcmp(stimulus_key_temp, 'choice_agreeable')
       chatbot_CA_start_relative = [chatbot_CA_start_relative stimulus_onset_temp];
       chatbot_CA_duration = [chatbot_CA_duration stimulus_dur_temp];
       chatbot_CA_amp = [chatbot_CA_amp zeros(1, length(stimulus_amp_temp)) + 1];
    end
    
    if strcmp(stimulus_key_temp, 'choice_neutral')
       chatbot_CN_start_relative = [chatbot_CN_start_relative stimulus_onset_temp];
       chatbot_CN_duration = [chatbot_CN_duration stimulus_dur_temp];
       chatbot_CN_amp = [chatbot_CN_amp zeros(1, length(stimulus_amp_temp)) + 1];
    end
end

for idx_se = 1:length(hb_second.stimulus.keys)
    stimulus_key_temp = hb_second.stimulus.values{1,idx_se}.name;
    stimulus_onset_temp = hb_second.stimulus.values{1,idx_se}.onset + relative_firstsession_duration;
    stimulus_dur_temp = hb_second.stimulus.values{1,idx_se}.dur;
    stimulus_amp_temp = hb_second.stimulus.values{1,idx_se}.amp;
    if strcmp(stimulus_key_temp, 'open_agreeable')
       chatbot_OA_start_relative = [chatbot_OA_start_relative stimulus_onset_temp];
       chatbot_OA_duration = [chatbot_OA_duration stimulus_dur_temp];
       chatbot_OA_amp = [chatbot_OA_amp zeros(1, length(stimulus_amp_temp)) + 1];     
    end
    
    if strcmp(stimulus_key_temp, 'open_neutral')
       chatbot_ON_start_relative = [chatbot_ON_start_relative stimulus_onset_temp];
       chatbot_ON_duration = [chatbot_ON_duration stimulus_dur_temp];
       chatbot_ON_amp = [chatbot_ON_amp zeros(1, length(stimulus_amp_temp)) + 1];  
    end
    
    if strcmp(stimulus_key_temp, 'choice_agreeable')
       chatbot_CA_start_relative = [chatbot_CA_start_relative stimulus_onset_temp];
       chatbot_CA_duration = [chatbot_CA_duration stimulus_dur_temp];
       chatbot_CA_amp = [chatbot_CA_amp zeros(1, length(stimulus_amp_temp)) + 1];
    end
    
    if strcmp(stimulus_key_temp, 'choice_neutral')
       chatbot_CN_start_relative = [chatbot_CN_start_relative stimulus_onset_temp];
       chatbot_CN_duration = [chatbot_CN_duration stimulus_dur_temp];
       chatbot_CN_amp = [chatbot_CN_amp zeros(1, length(stimulus_amp_temp)) + 1];
    end
end
for idx_se = 1:length(hb_third.stimulus.keys)
    stimulus_key_temp = hb_third.stimulus.values{1,idx_se}.name;
    stimulus_onset_temp = hb_third.stimulus.values{1,idx_se}.onset + relative_firstaddsecond_session_duration;
    stimulus_dur_temp = hb_third.stimulus.values{1,idx_se}.dur;
    stimulus_amp_temp = hb_third.stimulus.values{1,idx_se}.amp;
    if strcmp(stimulus_key_temp, 'open_agreeable')
       chatbot_OA_start_relative = [chatbot_OA_start_relative stimulus_onset_temp];
       chatbot_OA_duration = [chatbot_OA_duration stimulus_dur_temp];
       chatbot_OA_amp = [chatbot_OA_amp zeros(1, length(stimulus_amp_temp)) + 1];     
    end
    
    if strcmp(stimulus_key_temp, 'open_neutral')
       chatbot_ON_start_relative = [chatbot_ON_start_relative stimulus_onset_temp];
       chatbot_ON_duration = [chatbot_ON_duration stimulus_dur_temp];
       chatbot_ON_amp = [chatbot_ON_amp zeros(1, length(stimulus_amp_temp)) + 1];  
    end
    
    if strcmp(stimulus_key_temp, 'choice_agreeable')
       chatbot_CA_start_relative = [chatbot_CA_start_relative stimulus_onset_temp];
       chatbot_CA_duration = [chatbot_CA_duration stimulus_dur_temp];
       chatbot_CA_amp = [chatbot_CA_amp zeros(1, length(stimulus_amp_temp)) + 1];
    end
    
    if strcmp(stimulus_key_temp, 'choice_neutral')
       chatbot_CN_start_relative = [chatbot_CN_start_relative stimulus_onset_temp];
       chatbot_CN_duration = [chatbot_CN_duration stimulus_dur_temp];
       chatbot_CN_amp = [chatbot_CN_amp zeros(1, length(stimulus_amp_temp)) + 1];
    end
end
for idx_se = 1:length(hb_fourth.stimulus.keys)
    stimulus_key_temp = hb_fourth.stimulus.values{1,idx_se}.name;
    stimulus_onset_temp = hb_fourth.stimulus.values{1,idx_se}.onset + relative_firstsecondandthird_session_duration;
    stimulus_dur_temp = hb_fourth.stimulus.values{1,idx_se}.dur;
    stimulus_amp_temp = hb_fourth.stimulus.values{1,idx_se}.amp;
    if strcmp(stimulus_key_temp, 'open_agreeable')
       chatbot_OA_start_relative = [chatbot_OA_start_relative stimulus_onset_temp];
       chatbot_OA_duration = [chatbot_OA_duration stimulus_dur_temp];
       chatbot_OA_amp = [chatbot_OA_amp zeros(1, length(stimulus_amp_temp)) + 1];     
    end
    
    if strcmp(stimulus_key_temp, 'open_neutral')
       chatbot_ON_start_relative = [chatbot_ON_start_relative stimulus_onset_temp];
       chatbot_ON_duration = [chatbot_ON_duration stimulus_dur_temp];
       chatbot_ON_amp = [chatbot_ON_amp zeros(1, length(stimulus_amp_temp)) + 1];  
    end
    
    if strcmp(stimulus_key_temp, 'choice_agreeable')
       chatbot_CA_start_relative = [chatbot_CA_start_relative stimulus_onset_temp];
       chatbot_CA_duration = [chatbot_CA_duration stimulus_dur_temp];
       chatbot_CA_amp = [chatbot_CA_amp zeros(1, length(stimulus_amp_temp)) + 1];
    end
    
    if strcmp(stimulus_key_temp, 'choice_neutral')
       chatbot_CN_start_relative = [chatbot_CN_start_relative stimulus_onset_temp];
       chatbot_CN_duration = [chatbot_CN_duration stimulus_dur_temp];
       chatbot_CN_amp = [chatbot_CN_amp zeros(1, length(stimulus_amp_temp)) + 1];
    end
end
temp_se = nirs.design.StimulusEvents( 'open_agreeable', chatbot_OA_start_relative, chatbot_OA_duration, chatbot_OA_amp);
temp_merge('open_agreeable') = temp_se;
temp_se = nirs.design.StimulusEvents( 'open_neutral', chatbot_ON_start_relative, chatbot_ON_duration, chatbot_ON_amp);
temp_merge('open_neutral') = temp_se;
temp_se = nirs.design.StimulusEvents( 'choice_agreeable', chatbot_CA_start_relative, chatbot_CA_duration, chatbot_CA_amp);
temp_merge('choice_agreeable') = temp_se;
temp_se = nirs.design.StimulusEvents( 'choice_neutral', chatbot_CN_start_relative, chatbot_CN_duration, chatbot_CN_amp);
temp_merge('choice_neutral') = temp_se;
%%
temp_description = hb_first.description;
temp_data = cat(1, hb_first.data, hb_second.data, hb_third.data, hb_fourth.data);
temp_probe = hb_first.probe;
temp_time = 0:1/hb_first.Fs:(length(temp_data)-1)/hb_first.Fs;
temp_Fm = hb_first.Fm;
temp_auxillary = hb_first.auxillary;
temp_stimulus = temp_merge;
temp_demo = hb_first.demographics;
temp_Fs = hb_first.Fs;
hb = nirs.core.Data(temp_data, temp_time, temp_probe, temp_Fm, temp_stimulus, temp_demo, temp_description);
for s = 1:size(hb)
% view each subject's graph
    hb(s).draw;
    saveas(gcf,[root_dir,'concatedata_afternormalization/sub',num2str(s),'.png'])
    close;
end

%% Block Average for visualizations only
mkdir([strcat(root_dir, 'blockAverage')])
for s = 1:size(hb)
    HRF=BlockAverage2(-5, 40, hb(s),s, root_dir);  %the parameters of start and end here is based on experimental protocols
               
    for i = 1:11
        saveas(gcf,[root_dir,'/blockAverage/sub',num2str(s),'_', num2str(i),'.png'])
        close;
    end
end
%% 
% job1: GLM only
job1 = nirs.modules.GLM();
job1.type = 'AR-IRLS';
%job1.trend_func = @(t) nirs.design.trend.dctmtx(t, 0.005);
job1.basis = Dictionary();
job1.basis('default') = nirs.design.basis.Canonical(); 
job1.basis('open_agreeable') = nirs.design.basis.Canonical();
job1.basis('open_neutral') = nirs.design.basis.Canonical();
job1.basis('choice_agreeable') = nirs.design.basis.Canonical();
job1.basis('choice_neutral') = nirs.design.basis.Canonical();
job1.AddShortSepRegressors = true;
SubjStats = job1.run( hb );
    
%% the parameters can be modified in nirs.design.basis.Canonical
%% save GLM model
mkdir([strcat(root_dir, 'GLM')])
for s = 1:size(transpose(SubjStats))
    writetable(SubjStats(s).table,[root_dir,'/GLM/sub',num2str(s),'.txt'], 'Delimiter', ' '); 
    writetable(SubjStats(s).table,[root_dir,'/GLM/sub',num2str(s),'.xls']);
    save('/Users/xuhan/Desktop/fNIRS_chatbot/realresults/fNIRS/s16_chatbot.mat', 'SubjStats');
end

%% %% data view on individual level
mkdir([strcat(root_dir, 'data_individual')])
for s = 1:size(transpose(SubjStats))
%for s = 1:1
    %SubjStats(s).probe.defaultdrawfcn='3D mesh'/'10-20';  % cannot work, error message 'No public field defaultdrawfcn exists for class nirs.core.Probe.'
    SubjStats(s).probe.defaultdrawfcn='3D-mesh';
    SubjStats(s).draw('tstat', [], 'q < 0.05');
    for i = 1:14
        saveas(gcf,[root_dir,'/data_individual/sub',num2str(s),'_', num2str(i),'.png'])
        close;
    end
    writetable(SubjStats(s).table,[root_dir,'/data_individual/SubjStatssub',num2str(s),'.xls']);
end

%% contrast on individual level and ROI this is only IC and Infant Noise
%% display conditions
disp(SubjStats(1).conditions)
%% Define some contrasts
c = [-1 0 1 0]

mkdir([strcat(root_dir, 'contrast_individual')])
for s = 1:size(transpose(SubjStats))
    ContrastStats = SubjStats(s).ttest(c);
    % Display the contrasts
    ContrastStats.draw('tstat', [], 'q < 0.05');
    
    for i = 1:2
        saveas(gcf,[root_dir,'/contrast_individual/sub',num2str(s),'_', num2str(i),'.png'])
        close;
    end
end

%% export for 2-nd level analysis
mkdir([strcat(root_dir, 'secondlevelAnalysis_exported')])
T_temp = table(repmat(SubjStats(1).demographics.values(1),160,1), SubjStats(1).variables.source, SubjStats(1).variables.detector, SubjStats(1).variables.type, SubjStats(1).variables.cond, SubjStats(1).beta, 'VariableNames', {'subjectID','source','detector', 'type', 'cond', 'beta'});
for sID = 2:size(transpose(SubjStats)) 
    T_temp_loop = table(repmat(SubjStats(sID).demographics.values(1),160,1), SubjStats(sID).variables.source, SubjStats(sID).variables.detector, SubjStats(sID).variables.type, SubjStats(sID).variables.cond, SubjStats(sID).beta, 'VariableNames', {'subjectID', 'source','detector', 'type', 'cond', 'beta'}); 
    T_temp = [T_temp ; T_temp_loop]; 
end
writetable(T_temp,[root_dir, '/secondlevelAnalysis_exported/firstchatbotsection_exportedforSPSS.xls'])