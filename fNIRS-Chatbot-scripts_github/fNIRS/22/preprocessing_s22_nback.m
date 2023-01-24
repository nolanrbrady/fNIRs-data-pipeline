clear;
clc;

%% group number
group_number = 1;

%% load data
data_dir = '/Users/xuhan/Desktop/fNIRS_chatbot/realresults/fNIRS/22/2022-05-02_004/'
raw = nirs.io.loadDirectory(data_dir, {'group','subject'});

root_dir = '/Users/xuhan/Desktop/fNIRS_chatbot/realresults/fNIRS/22/nback/'

%% Create demographics table

demographics = nirs.createDemographicsTable(raw);
disp(demographics)

%% Manually add events info
timestamp_start_nback = '14:16:30'
nback1 = '1_back';
timestamp_nback1_start = '14:20:45';
timestamp_nback1_end = '14:21:27';

nback1_survey = '1_back_survey';
timestamp_nback1_survey_start = '14:21:50';
timestamp_nback1_survey_end = '14:22:47';

nback2 = '1_back';
timestamp_nback2_start = '14:23:10';
timestamp_nback2_end = '14:23:52';

nback2_survey = '1_back_survey';
timestamp_nback2_survey_start = '14:24:15';
timestamp_nback2_survey_end = '14:24:51';

nback3 = '1_back';
timestamp_nback3_start = '14:25:14';
timestamp_nback3_end = '14:25:57';

nback3_survey = '1_back_survey';
timestamp_nback3_survey_start = '14:26:20';
timestamp_nback3_survey_end = '14:27:04';

nback4 = '3_back';
timestamp_nback4_start = '14:27:27';
timestamp_nback4_end = '14:28:09';

nback4_survey = '3_back_survey';
timestamp_nback4_survey_start = '14:28:32';
timestamp_nback4_survey_end = '14:29:13';

nback5 = '3_back';
timestamp_nback5_start = '14:29:36';
timestamp_nback5_end = '14:30:16';

nback5_survey = '3_back_survey';
timestamp_nback5_survey_start = '14:30:39';
timestamp_nback5_survey_end = '14:31:22';

nback6 = '3_back';
timestamp_nback6_start = '14:31:45';
timestamp_nback6_end = '14:32:27';

%% add eventinfo
temp=Dictionary();

nback_keys = string({nback1; nback2; nback3; nback4; nback5; nback6; nback1_survey; nback2_survey; nback3_survey; nback4_survey; nback5_survey});
nback_start = [timestamp_nback1_start; timestamp_nback2_start; timestamp_nback3_start; timestamp_nback4_start; timestamp_nback5_start; timestamp_nback6_start; timestamp_nback1_survey_start; timestamp_nback2_survey_start; timestamp_nback3_survey_start; timestamp_nback4_survey_start; timestamp_nback5_survey_start];
nback_end = [timestamp_nback1_end; timestamp_nback2_end; timestamp_nback3_end; timestamp_nback4_end; timestamp_nback5_end; timestamp_nback6_end; timestamp_nback1_survey_end; timestamp_nback2_survey_end; timestamp_nback3_survey_end; timestamp_nback4_survey_end; timestamp_nback5_survey_end];
oneback_start_relative = [];
oneback_duration = [];
oneback_amp = [];
threeback_start_relative = [];
threeback_duration = [];
threeback_amp = [];
oneback_survey_start_relative = [];
oneback_survey_duration = [];
oneback_survey_amp = [];
threeback_survey_start_relative = [];
threeback_survey_duration = [];
threeback_survey_amp = [];


for nback_idx = 1:11
    nback_start_relative_temp = seconds(datetime(nback_start(nback_idx,:),'InputFormat','HH:mm:ss') - datetime(timestamp_start_nback,'InputFormat','HH:mm:ss'));
    nback_duration_temp = seconds(datetime(nback_end(nback_idx,:),'InputFormat','HH:mm:ss') - datetime(nback_start(nback_idx,:),'InputFormat','HH:mm:ss'));
    if strcmp(nback_keys(nback_idx,:), '1_back')
       oneback_start_relative = [oneback_start_relative nback_start_relative_temp];
       oneback_duration = [oneback_duration nback_duration_temp];
       oneback_amp = [oneback_amp 1];
       temp_se = nirs.design.StimulusEvents( '1_back', oneback_start_relative, oneback_duration, oneback_amp);
       temp('1_back') = temp_se;
    end
    
    if strcmp(nback_keys(nback_idx,:), '3_back')
       threeback_start_relative = [threeback_start_relative nback_start_relative_temp];
       threeback_duration = [threeback_duration nback_duration_temp];
       threeback_amp = [threeback_amp 1];
       temp_se = nirs.design.StimulusEvents( '3_back', threeback_start_relative, threeback_duration, threeback_amp);
       temp('3_back') = temp_se;
    end
    
    if strcmp(nback_keys(nback_idx,:), '1_back_survey')
       oneback_survey_start_relative = [oneback_survey_start_relative nback_start_relative_temp];
       oneback_survey_duration = [oneback_survey_duration nback_duration_temp];
       oneback_survey_amp = [oneback_survey_amp 1];
       temp_se = nirs.design.StimulusEvents( '1_back_survey', oneback_survey_start_relative, oneback_survey_duration, oneback_survey_amp);
       temp('1_back_survey') = temp_se;
    end
    
    if strcmp(nback_keys(nback_idx,:), '3_back_survey')
       threeback_survey_start_relative = [threeback_survey_start_relative nback_start_relative_temp];
       threeback_survey_duration = [threeback_survey_duration nback_duration_temp];
       threeback_survey_amp = [threeback_survey_amp 1];
       temp_se = nirs.design.StimulusEvents( '3_back_survey', threeback_survey_start_relative, threeback_survey_duration, threeback_survey_amp);
       temp('3_back_survey') = temp_se;
    end
end

raw(1,1).stimulus = temp;


%% Quality checking 
mkdir([strcat(root_dir, 'rawdata')])
mkdir([strcat(root_dir, 'qualitycheck')])

for s = 1:size(raw)
% view each subject's graph
    data = raw(s);
    data.draw;
    saveas(gcf,[root_dir,'/rawdata/sub',num2str(s),'.png'])
    close;
    
% turncate for quality checking
    pre_Baseline = 23;
    post_Baseline = 23;
    fs = data.Fs;
    [m,n] = size(data.stimulus.keys);
   
    onset_firststi = timestamp_nback1_start;
    onset_laststi = timestamp_nback6_end;
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
    
    data = raw(s).data(round(temp_start*fs):round(temp_end*fs),:);

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
            disp(raw(s).probe.link(i,:));
        end
    end
    T = table(table2array(raw(s).probe.link(:,1)), table2array(raw(s).probe.link(:,2)),table2array(raw(s).probe.link(:,3)), transpose(cvMask));
    xlswrite([root_dir,'qualitycheck/sub',num2str(s),'.csv'], table2array(T));

end

%% Prepare the data for preprocessing
%% Remove stimless files
mkdir([strcat(root_dir, 'stimremove_turncate')])
j = nirs.modules.RemoveStimless();

j = nirs.modules.KeepStims( j);
j.listOfStims = {'3_back', '1_back'}; 


% resample to 4Hz
j = nirs.modules.Resample( j );
j.Fs = 4;  % Sets the new sample rate to 4Hz (was 10Hz).

% Trim pre and post baseline
j = nirs.modules.TrimBaseline( j);
j.preBaseline = 23; %can change these values
j.postBaseline = 23;
raw = j.run(raw);

for s = 1:size(raw)
% view each subject's graph
    raw(s).draw;
    saveas(gcf,[root_dir,'stimremove_turncate/sub',num2str(s),'.png'])
    close;
end


%% %% Signal Preprocessing
mkdir([strcat(root_dir, 'preprocessing')])
hb = raw;
% visualize results
for s = 1:size(hb)
% view each subject's graph
    hb(s).draw;
    saveas(gcf,[root_dir,'preprocessing/afterBP_sub_',num2str(s),'.png'])
    close;
end

%% SS label
job_SS = nirs.modules.LabelShortSeperation();
hb = job_SS.run(hb);
%% Run the conversions 
mkdir([strcat(root_dir, 'conversion')])
j = nirs.modules.OpticalDensity();
% the DFP parameters can be modified in nirs.modules.BeerLamberLaw.m

%% %% Convert to hemoglobin
j = nirs.modules.BeerLambertLaw( j);    
hb = j.run(hb);
% Results visualization
for s = 1:size(hb)
% view each subject's graph
    hb(s).draw;
    saveas(gcf,[root_dir,'conversion/sub',num2str(s),'.png'])
    close;
end

%% SS filter
%job_SSfilter = advanced.nirs.modules.ShortDistanceFilter();
%hb = job_SSfilter.run(hb);

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
job1.basis('3_back') = nirs.design.basis.Canonical();
job1.basis('1_back') = nirs.design.basis.Canonical();
job1.AddShortSepRegressors = true;
SubjStats = job1.run( hb );
    
%% the parameters can be modified in nirs.design.basis.Canonical
%% save GLM model
mkdir([strcat(root_dir, 'GLM')])
for s = 1:size(transpose(SubjStats))
    writetable(SubjStats(s).table,[root_dir,'/GLM/sub',num2str(s),'.txt'], 'Delimiter', ' '); 
    writetable(SubjStats(s).table,[root_dir,'/GLM/sub',num2str(s),'.xls']);
    save('/Users/xuhan/Desktop/fNIRS_chatbot/realresults/fNIRS/s22_nback.mat', 'SubjStats');
end

%% %% data view on individual level
mkdir([strcat(root_dir, 'data_individual')])
for s = 1:size(transpose(SubjStats))
%for s = 1:1
    %SubjStats(s).probe.defaultdrawfcn='3D mesh'/'10-20';  % cannot work, error message 'No public field defaultdrawfcn exists for class nirs.core.Probe.'
    SubjStats(s).probe.defaultdrawfcn='3D-mesh';
    SubjStats(s).draw('tstat', [], 'q < 0.05');
    for i = 1:12
        saveas(gcf,[root_dir,'/data_individual/sub',num2str(s),'_', num2str(i),'.png'])
        close;
    end
    writetable(SubjStats(s).table,[root_dir,'/data_individual/SubjStatssub',num2str(s),'.xls']);
end

%% contrast on individual level and ROI this is only IC and Infant Noise
%% display conditions
disp(SubjStats(1).conditions)
%% Define some contrasts
c = [-1 1]

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
T_temp = table(repmat(SubjStats(1).demographics.values(1),80,1), SubjStats(1).variables.source, SubjStats(1).variables.detector, SubjStats(1).variables.type, SubjStats(1).variables.cond, SubjStats(1).beta, 'VariableNames', {'subjectID','source','detector', 'type', 'cond', 'beta'});
for sID = 2:size(transpose(SubjStats)) 
    T_temp_loop = table(repmat(SubjStats(sID).demographics.values(1),80,1), SubjStats(sID).variables.source, SubjStats(sID).variables.detector, SubjStats(sID).variables.type, SubjStats(sID).variables.cond, SubjStats(sID).beta, 'VariableNames', {'subjectID', 'source','detector', 'type', 'cond', 'beta'}); 
    T_temp = [T_temp ; T_temp_loop]; 
end
writetable(T_temp,[root_dir, '/secondlevelAnalysis_exported/nback_exportedforSPSS.xls'])