clc;
clear;

%% load first level stats
load('/Users/xuhan/Desktop/fNIRS_chatbot/realresults/fNIRS/all_chatbot.mat', 'SubjStats');

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
c = [-1 -1 1 1];
%%
GroupStats.probe.defaultdrawfcn='3D-mesh';
GroupStats.draw('tstat', [], 'p < 0.05')
%%
GroupStats.probe.defaultdrawfcn='3D-mesh';
GroupStats.draw('tstat', [], 'q < 0.05')

%%
ContrastStats = GroupStats.ttest(c);
ContrastStats.draw('tstat', [], 'q< 0.05')
%%
%% export for 2-nd level analysis
T_temp = table(GroupStats.variables.source, GroupStats.variables.detector, GroupStats.variables.type, GroupStats.variables.cond, GroupStats.beta, 'VariableNames', {'source','detector', 'type', 'cond', 'beta'});
writetable(T_temp,'/Users/xuhan/Desktop/fNIRS_chatbot/realresults/MixedEffects_chatbot.csv')
