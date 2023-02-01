clc;
clear;

%% load first level stats
load('/Users/xuhan/Desktop/fNIRS_chatbot/realresults/fNIRS/all_nback.mat', 'SubjStats');

%%
job = nirs.modules.MixedEffects();
job.formula = 'beta ~ -1 + cond + (1|subject)';
% if to go with both main effects and interactions, do the following
% formula 
% job.formula = 'beta ~ cond + group + session + cond:group + cond:session +
% group:session + cond:group:session + (1|subject)'
job.dummyCoding = 'full';
% if to go with both main effects and interactions, do the following
% dummyCoding
% "effects"
job.include_diagnostics=true;
GroupStats = job.run(SubjStats);

%%
disp(GroupStats.conditions);
%%
c = [-1 1];
%%
GroupStats.probe.defaultdrawfcn='3D-mesh';
GroupStats.draw('tstat', [], 'p < 0.05')
%%
GroupStats.probe.defaultdrawfcn='3D-mesh';
GroupStats.draw('tstat', [], 'q < 0.05')

%%
ContrastStats = GroupStats.ttest(c);
ContrastStats.draw('tstat', [], 'q < 0.05')

%% ROI Analysis
MeasList=[2 1;...
          2 3;...
          3 3;...
          3 1];
      
Region{1} = table(MeasList(:,1),MeasList(:,2),'VariableNames',{'source','detector'});
ROItable=nirs.util.roiAverage(GroupStats,Region,{'region1'});
disp(ROItable);

job_ROI = nirs.modules.ApplyROI();
job_ROI.listOfROIs = table(MeasList(1,1),MeasList(1,2),{'Edge1'},'VariableNames',{'source','detector','name'});
job_ROI.listOfROIs(1,:) = table(MeasList(1,1),MeasList(1,2),{'Edge1'},'VariableNames',{'source','detector','name'});
job_ROI.listOfROIs(2,:) = table(MeasList(2,1),MeasList(2,2),{'Edge2'},'VariableNames',{'source','detector','name'});
job_ROI.listOfROIs(3,:) = table(MeasList(3,1),MeasList(3,2),{'Edge3'},'VariableNames',{'source','detector','name'});
job_ROI.listOfROIs(4,:) = table(MeasList(4,1),MeasList(4,2),{'Edge4'},'VariableNames',{'source','detector','name'});
job.weighted = false;
dataROI = job_ROI.run( GroupStats );
%%
ContrastStats_ROI = dataROI.ttest(c);
ContrastStats_ROI.draw('tstat', [], 'q < 0.05')


