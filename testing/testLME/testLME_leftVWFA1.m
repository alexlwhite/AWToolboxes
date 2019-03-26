%Script to compare fitlme and rm_anova2 on a subset of data from VWFA_Attn5
% Alex White, 8/8/18

clear; close all;

%% load the data table
load('LeftVWFA1_channelResps.mat');
%this table contains "channel responses" in the variable "resp" as
%estimated with an "inverted encoding model" from the left hemisphere VWFA.
%For each subject, each channel's response was estimated in two conditions: attend-left and
%attend-right (coded as 'singTaskL' and 'singTaskR'). 

%The question is: does the focus of attention (left or right) modulate the channel response
%magnitudes, and does it do so differently in the two channels? 

%% fit LME
%Model the responses as an effect of channel (left vs right), cue (attend left
%vs attend right) and their interaction, with subject as a random effect:
eqtn = 'resp ~ channel * cue + (1 | subject)';
lme  = fitlme(T,eqtn,'DummyVarCoding','effects');

fprintf(1,'\n-------------------------------');
fprintf(1,'\nOutput of fitlme for effects of channel*cue:\n');
lme.disp;
display(lme.anova);

%This suggests that there's an effect of channel (R>L) but no effect of cue, and no interaction  

%% try Matlab's anovan 
[pval,anovaTable] = anovan(T.resp,{T.channel, T.cue, T.subject},'model','interaction','varnames',{'channel','cue','subject'},'random',3,'display','off');

%This suggests that there is no main effect of channel but there is a
%(barely significant) interaction. What gives?
%(note: the function rm_anova2 gives the same results)

%display subset of ANOVA table
anovaTable = anovaTable(1:5,:); 
anovaTable = anovaTable(:,1:7);
fprintf(1,'\n-------------------------------');
fprintf(1,'\nOutput of anovan:\n');
disp(anovaTable);

%% bootstrap the effect of cue within each channel, and plot them
nSubjs     = max(T.subject);
rs         = NaN(2,2,nSubjs);
cueEffects = NaN(2,nSubjs);
mrs        = NaN(2,2);
sems       = mrs;
cueCIs     = NaN(2,1,2);

channelNames = unique(T.channel); 
cueNames     = unique(T.cue);

%re-format the data and bootstrap 95% CIs of cue effects
for ki=1:2 %channel
    for ci=1:2 %attn cond
        theseRs     = T.resp(strcmp(T.channel,channelNames{ki}) & strcmp(T.cue, cueNames{ci}));
        rs(ki,ci,:) = theseRs;
        mrs(ki,ci)  = nanmean(theseRs);
        sems(ki,ci) = standardError(theseRs');
    end
    cueEffects(ki,:)   = squeeze(diff(rs(ki,:,:),1,2));
    cueCIs(ki,1,:)     = boyntonBootstrap(@nanmean,cueEffects(ki,:),1000,95);
end

%bar plot
opt = struct;
opt.xTickLabs  = channelNames;
opt.xLab       = 'Channel';
opt.yLab        = 'channel response';
opt.ylims      = [-0.1 1];
opt.doLegend   = true;
opt.legendLabs = cueNames;
figure; subplot(2,1,1); hold on;
barPlot_AW(mrs,sems,opt);
title('Mean channel responses +/- SEM');

opt.doLegend = false;
opt.ylims = [-0.4 0.4];
opt.fillColors(1,1,:) = [0 0.7 0];
opt.fillColors(2,1,:) = [0 0.7 0];
opt.errorBarColors(1,1,:) = [0 0 0];
opt.errorBarColors(2,1,:) = [0 0 0];
opt.yLab = sprintf('%s - %s', cueNames{2},cueNames{1});
subplot(2,1,2); hold on;
barPlot_AW(mean(cueEffects,2), cueCIs, opt);
title('Mean cue effects w/ 95% CIs');

%% Try an LME model on each channel individually 
eqtn = 'resp ~ cue + (1 | subject)';
lmeL = fitlme(T(strcmp(T.channel,'left'),:),eqtn,'DummyVarCoding','effects'); %the result is the same with 'reference' DummyVarCoding

fprintf(1,'\n-------------------------------');
fprintf(1,'\nOutput of fitlme for effects of cue in Left channel only:\n');
lmeL.disp;
display(lmeL.anova);

lmeR = fitlme(T(strcmp(T.channel,'right'),:),eqtn,'DummyVarCoding','effects'); %the result is the same with 'reference' DummyVarCoding
fprintf(1,'\n-------------------------------');
fprintf(1,'\nOutput of fitlme for effects of cue in Right channel only:\n');
lmeR.disp;
display(lmeR.anova);

%So these simpler models on a subset of the data agree with the boostrapping and the ANOVA that there is
%an effect of cue within each channel individually, and they go in opposite
%directions. 

%So shouldn't there then be an interaction in the full LME model? 
%does it have to do with the fact that the LME fit is max likelihood,
%whereas the ANOVA is least-squares? 
