function lc_resilience(whichRoi)

% % %  Description:
% % %  if whichRoi is defined as 3 this script reproduces manuscript figure 1 panel e and figure 3 panels a-f
% % %  if other rois are defined the same plots are produced for the chosen roi data
% % %  please see below for roi options
% % %  please note: all data are physio-corrected, unsmoothed and weighted avegrages %%%
% % % 
% % %  Usage example >> lc_resilience(3)
% % % 
% % %  whichRoi options:
% % % 1 = 'LC-2SD_Binary';
% % % 2 = 'LC (2SD)';
% % % 3 = 'LC (1SD)';
% % % 4 = 'Dorsal Raphe (DR)';
% % % 5 = 'Medial Raphe (MR)';
% % % 6 = 'Amygdala';
% % % 7 = 'Substantia Nigra (SN)';
% % % 8 = 'Ventral Tegmental Area (VTA)';

% % % Author: Marcus Grueschow 
% % % Date: 2019/20

if nargin < 1
    whichRoi = 3; %%% 'LC (1SD)' %%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STATS PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
statsTextFontSize = 18;
axisLabelFontSize = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CONTRAST STATS AND SYMPTOM CORRELATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta1_stai =   [2     3     0     0    -3    -4     4   -13     9    -1     1     2    -4   -12    -3    -4     1    -4     7    -6   -38    -5    -1     3     5    -1    -6     0    -4    -2   -16    -2    -2     6    14    -8    -3     0    -2    -2     2     7    -5    -6   -13     8    -4    10];
delta2_stai =   [3    16   -13     1    -2     1    -5     2    10     0     5     1    -4   -16    -4    -4    -5     1     6    -7   -29   -15    -1    -1    -1     0     4     3     4     4    -5     1    -3     0     8    -2     3     5     5    -6     6     9     4     1    -7   -18    -2    -2];
t0_stai     = [32    39    27    33    28    38    41    28    58    30    24    48    26    41    40    29    62    24    29    28    24    38    33    43    32    23    31    29    38    31    37    23    31    36    51    42    34    51    35    35    38    35    33    26    36    35    24    30];
t1_stai     = t0_stai + delta1_stai;
t2_stai     = t0_stai + delta2_stai;
delta12_stai = mean([delta1_stai; delta2_stai]);

%%% PHQ - Depression
delta1_phq =    [4     4    -2    -2    -3     4     3    -4    -1     0     1     0    -2     1    -4    -3     1    -2     2    -5   -10    -5    -1     0    -1     1    -4    -2    -3     0    -2    -1     1     1     0    -1     1   -11     8     1     0     1     4    -2   -11    -2    -2     2];
delta2_phq =    [3     7    -3    -2    -2     5    -4     0     4     0     3    -2    -3    -3    -1    -2    -2    -2     1    -6    -6     1    -1    -1    -3     1     0    -1    -2     1     0     1    -1     0    -1    -4     0    -7     8     1    -2     2     2    -3    -2   -13    -3     0];
t0_phq     = [15    17     9    10     9    18    13     9    16    10    12    12     9    19    13    10    24    10    12     9     9    12    12    14    10    11    11    11     9    13    12    10    12    16    15    16    13    15    21    13    12    11    14     9    11     9     9    13];
t1_phq     = t0_phq + delta1_phq;
t2_phq     = t0_phq + delta2_phq;
delta12_phq = mean([delta1_phq; delta2_phq]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DEFINE DATA SET (physio-corrected, unsmoothed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
currDataString = 'weightedAvg_PhysioCorrected_UnsmoothedData';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NOW CI>II and Symptom CORRELATION ROI-wise:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Subcortical Control Analyses. Data used --> ', currDataString]);
currRoiNameCell = {'LC_2SD_Binary', 'LC_2SD', 'LC_1SD', 'MR', 'DR', 'AMY', 'SN', 'VTA'};
probCut=0.001; %%% for more stringent probability cut-off increase this variable;

for currRoi = whichRoi
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% GET CURRENT ROI NAME
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if currRoi == 1 %%% original LC-ROI
        currRoiTitle = 'LC-2SD_Binary';
        currRoiColor = [61 245 96]/255;
    elseif currRoi == 2
        currRoiTitle = 'LC (2SD)';
        currRoiColor = [144 219 69]/255;
    elseif currRoi == 3
        currRoiTitle = 'LC (1SD)';
        currRoiColor = [0 176 80]/255;
    elseif currRoi == 4
        currRoiTitle = 'Dorsal Raphe (DR)';
        currRoiColor = [27 83 241]/255;
    elseif currRoi == 5
        currRoiTitle = 'Medial Raphe (MR)';
        currRoiColor = [248 242 0]/255;
    elseif currRoi == 6
        currRoiTitle = 'Amygdala';
        currRoiColor = [184 242 250]/255;
    elseif currRoi == 7
        currRoiTitle = 'Substantia Nigra (SN)';
        currRoiColor = [45 213 231]/255;
    elseif currRoi == 8
        currRoiTitle = 'Ventral Tegmental Area (VTA)';
        currRoiColor = [245 39 230]/255;
    end
    
    currRoiName = currRoiNameCell{currRoi};
    load(['allSubsWeightedAvgVals_', currRoiName, '_', currDataString]); %%% yields allSubsWeightedAvgVals %%%
    load(['allSubsAllRoiVoxCoordsProbValsConVals_', currRoiName, '_', currDataString]) %%% yields allSubsAllRoiVoxCoordsProbValsConVals %%%
    
    allSubsWeightedAvgVals = nan(48,1);
    mySubsStartRow = [1 : size(allSubsAllRoiVoxCoordsProbValsConVals, 1)/48 : size(allSubsAllRoiVoxCoordsProbValsConVals, 1),  size(allSubsAllRoiVoxCoordsProbValsConVals, 1)];
    for currSub = 1 : 48
        currSubCurrRoiVoxRows = allSubsAllRoiVoxCoordsProbValsConVals(mySubsStartRow(currSub) : mySubsStartRow(currSub+1)-1, :);
        currSubCurrRoiVoxRowsCut = currSubCurrRoiVoxRows(currSubCurrRoiVoxRows(:,4) > probCut, :);
        currRoiCurrVoxCount = size(currSubCurrRoiVoxRowsCut, 1);
        %%% Weighted Mean here:
        allSubsWeightedAvgVals(currSub) = mean(currSubCurrRoiVoxRowsCut(:,4) .* currSubCurrRoiVoxRowsCut(:,5));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% ANXIETY Symptom-Changes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% 3 Months
    y= delta1_stai*-1;
    disp([currRoiName ': STAI-1: after 3 months']);
    fh = figure;
    currColor = [0.6 0.6 0.6]; set(gcf,'color','w'); txtColor = [.3 .3 .3];
    yData = y'; xData = allSubsWeightedAvgVals;
    [RHO,PVAL] = corr(xData, yData,'type','Spearman');
    [b,stats] = robustfit(xData, yData); disp(['Robust Fit: p = ', num2str(stats.p(2))]);
    [B,BINT,R,RINT,STATS] = regress( yData, [xData, ones(size(yData)) ]); %%% y = p1*x + p2  p1 = B(1); p2 = B(2)    %%%
    x = linspace(min(xData(:,1))-.1,  max(xData(:,1))+.1, 1000);
    plot(xData(:,1), yData, 'o','markerfacecolor', currColor, 'markeredgecolor', 'k', 'Markersize', 12); hold on; box off; axis square
    plot(x, B(1)*x+B(end), 'linewidth', 3, 'color', currColor/2); hold on
    set(gca, 'fontsize', 15); xlim([-.7 1]); ylim([-20 40])
    xlabel([currRoiTitle ' (CI>II)'], 'fontsize', axisLabelFontSize, 'FontAngle', 'italic');
    ylabel('Anxiety symptoms change', 'fontsize', axisLabelFontSize, 'FontAngle', 'italic');
    myR = num2str(sqrt(STATS(1,1))); myP = num2str(stats.p(2));
    myR = num2str(RHO); myP = num2str(stats.p(2)); th = text(-.5, 42, ['Rho = ' , myR(1:4), ', p = ', myP(1:5)], 'fontsize', statsTextFontSize, 'color', txtColor, 'FontAngle', 'italic'); hold on;
    set(fh, 'position', [30 600 500 400]);

    %%% 6 Months
    y= delta2_stai*-1;
    disp([currRoiName ': STAI-2: after 6 months']);
    fh = figure;
    currColor = [0.6 0.6 0.6]; set(gcf,'color','w'); txtColor = [.3 .3 .3];
    yData = y'; xData = allSubsWeightedAvgVals;
   [RHO,PVAL] = corr(xData, yData,'type','Spearman'); 
    [b,stats] = robustfit(xData, yData); disp(['Robust Fit: p = ', num2str(stats.p(2))]);
    [B,BINT,R,RINT,STATS] = regress( yData, [xData, ones(size(yData)) ]); %%% y = p1*x + p2  p1 = B(1); p2 = B(2)    %%%
    x = linspace(min(xData(:,1))-.1,  max(xData(:,1))+.1, 1000);
    plot(xData(:,1), yData, 'o','markerfacecolor', currColor, 'markeredgecolor', 'k', 'Markersize', 12); hold on; box off; axis square
    plot(x, B(1)*x+B(end), 'linewidth', 3, 'color', currColor/2); hold on
    set(gca, 'fontsize', 15); xlim([-.7 1]); ylim([-20 40])
    xlabel([currRoiTitle ' (CI>II)'], 'fontsize', axisLabelFontSize, 'FontAngle', 'italic');
    ylabel('Anxiety symptoms change', 'fontsize', axisLabelFontSize, 'FontAngle', 'italic');
    myR = num2str(sqrt(STATS(1,1))); myP = num2str(stats.p(2));
    myR = num2str(RHO); myP = num2str(stats.p(2)); th = text(-.5, 42, ['Rho = ' , myR(1:4), ', p = ', myP(1:5)], 'fontsize', statsTextFontSize, 'color', txtColor, 'FontAngle', 'italic'); hold on;
    set(fh, 'position', [530 600 500 400]);
    
    %%% AVERAGE 3 & 6 Months
    y= delta12_stai*-1;
    disp([currRoiName ': STAI-12: Mean Changes']);
    fh = figure;
    currColor = [0.6 0.6 0.6]; set(gcf,'color','w'); txtColor = [.3 .3 .3];
    yData = y'; xData = allSubsWeightedAvgVals;
    [RHO,PVAL] = corr(xData, yData,'type','Spearman'); 
    [b,stats] = robustfit(xData, yData); disp(['Robust Fit: p = ', num2str(stats.p(2))]);
    [B,BINT,R,RINT,STATS] = regress( yData, [xData, ones(size(yData)) ]); %%% y = p1*x + p2  p1 = B(1); p2 = B(2)    %%%
    x = linspace(min(xData(:,1))-.1,  max(xData(:,1))+.1, 1000);
    plot(xData(:,1), yData, 'o','markerfacecolor', currColor, 'markeredgecolor', 'k', 'Markersize', 12); hold on; box off; axis square
    plot(x, B(1)*x+B(end), 'linewidth', 3, 'color', currColor/2); hold on
    set(gca, 'fontsize', 15); xlim([-.7 1]); ylim([-20 40])
    xlabel([currRoiTitle ' (CI>II)'], 'fontsize', axisLabelFontSize, 'FontAngle', 'italic');
    ylabel('Anxiety symptoms change', 'fontsize', axisLabelFontSize, 'FontAngle', 'italic');
    myR = num2str(sqrt(STATS(1,1))); myP = num2str(stats.p(2));
    myR = num2str(RHO); myP = num2str(stats.p(2)); th = text(-.5, 42, ['Rho = ' , myR(1:4), ', p = ', myP(1:5)], 'fontsize', statsTextFontSize, 'color', txtColor, 'FontAngle', 'italic'); hold on;
    set(fh, 'position', [1030 600 500 400]);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% DEPRESSION Symptom-Changes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% 3 Months
    y= delta1_phq*-1;
    disp([currRoiName ': PHQ-1: after 3 months']);
    fh = figure;
    currColor = [0.6 0.6 0.6]; set(gcf,'color','w'); txtColor = [.3 .3 .3];
    yData = y'; xData = allSubsWeightedAvgVals;
    [RHO,PVAL] = corr(xData, yData,'type','Spearman'); 
    [b,stats] = robustfit(xData, yData); disp(['Robust Fit: p = ', num2str(stats.p(2))]);
    [B,BINT,R,RINT,STATS] = regress( yData, [xData, ones(size(yData)) ]); %%% y = p1*x + p2  p1 = B(1); p2 = B(2)    %%%
    x = linspace(min(xData(:,1))-.1,  max(xData(:,1))+.1, 1000);
    plot(xData(:,1), yData, 'o','markerfacecolor', currColor, 'markeredgecolor', 'k', 'Markersize', 12); hold on; box off; axis square
    plot(x, B(1)*x+B(end), 'linewidth', 3, 'color', currColor/2); hold on
    set(gca, 'fontsize', 15); xlim([-.7 1]); ylim([-10 15])
    xlabel([currRoiTitle ' (CI>II)'], 'fontsize', axisLabelFontSize, 'FontAngle', 'italic');
    ylabel('Depression symptoms change', 'fontsize', axisLabelFontSize, 'FontAngle', 'italic');
    myR = num2str(sqrt(STATS(1,1))); myP = num2str(stats.p(2));
    myR = num2str(RHO); myP = num2str(stats.p(2)); th = text(-.5, 15, ['Rho = ' , myR(1:4), ', p = ', myP(1:5)], 'fontsize', statsTextFontSize, 'color', txtColor, 'FontAngle', 'italic'); hold on;
    set(fh, 'position', [30 100 500 400]);
    
    %%% 6 Months
    y= delta2_phq*-1;
    disp([currRoiName ': PHQ-2: after 6 months']);
    fh = figure;
    currColor = [0.6 0.6 0.6]; set(gcf,'color','w'); txtColor = [.3 .3 .3];
    yData = y'; xData = allSubsWeightedAvgVals;
    [RHO,PVAL] = corr(xData, yData,'type','Spearman'); 
    [b,stats] = robustfit(xData, yData); disp(['Robust Fit: p = ', num2str(stats.p(2))]);
    [B,BINT,R,RINT,STATS] = regress( yData, [xData, ones(size(yData)) ]); %%% y = p1*x + p2  p1 = B(1); p2 = B(2)    %%%
    x = linspace(min(xData(:,1))-.1,  max(xData(:,1))+.1, 1000);
    plot(xData(:,1), yData, 'o','markerfacecolor', currColor, 'markeredgecolor', 'k', 'Markersize', 12); hold on; box off; axis square
    plot(x, B(1)*x+B(end), 'linewidth', 3, 'color', currColor/2); hold on
    set(gca, 'fontsize', 15); xlim([-.7 1]); ylim([-10 15])
    xlabel([currRoiTitle ' (CI>II)'], 'fontsize', axisLabelFontSize, 'FontAngle', 'italic');
    ylabel('Depression symptoms change', 'fontsize', axisLabelFontSize, 'FontAngle', 'italic');
    myR = num2str(sqrt(STATS(1,1))); myP = num2str(stats.p(2));
    myR = num2str(RHO); myP = num2str(stats.p(2)); th = text(-.5, 15, ['Rho = ' , myR(1:4), ', p = ', myP(1:5)], 'fontsize', statsTextFontSize, 'color', txtColor, 'FontAngle', 'italic'); hold on;
    set(fh, 'position', [530 100 500 400]);
    
    %%% AVERAGE 3 & 6 Months
    y= delta12_phq*-1;
    disp([currRoiName ': PHQ-12: Mean Changes']);
    fh = figure;
    currColor = [0.6 0.6 0.6]; set(gcf,'color','w'); txtColor = [.3 .3 .3];
    yData = y'; xData = allSubsWeightedAvgVals;
    [RHO,PVAL] = corr(xData, yData,'type','Spearman'); 
    [b,stats] = robustfit(xData, yData); disp(['Robust Fit: p = ', num2str(stats.p(2))]);
    [B,BINT,R,RINT,STATS] = regress( yData, [xData, ones(size(yData)) ]); %%% y = p1*x + p2  p1 = B(1); p2 = B(2)    %%%
    x = linspace(min(xData(:,1))-.1,  max(xData(:,1))+.1, 1000);
    plot(xData(:,1), yData, 'o','markerfacecolor', currColor, 'markeredgecolor', 'k', 'Markersize', 12); hold on; box off; axis square
    plot(x, B(1)*x+B(end), 'linewidth', 3, 'color', currColor/2); hold on
    set(gca, 'fontsize', 15); xlim([-.7 1]); ylim([-10 15])
    xlabel([currRoiTitle ' (CI>II)'], 'fontsize', axisLabelFontSize, 'FontAngle', 'italic');
    ylabel('Depression symptoms change', 'fontsize', axisLabelFontSize, 'FontAngle', 'italic');
    myR = num2str(sqrt(STATS(1,1))); myP = num2str(stats.p(2));
    myR = num2str(RHO); myP = num2str(stats.p(2)); th = text(-.5, 15, ['Rho = ' , myR(1:4), ', p = ', myP(1:5)], 'fontsize', statsTextFontSize, 'color', txtColor, 'FontAngle', 'italic'); hold on;
    set(fh, 'position', [1030 100 500 400]);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% MEDIAN SPLIT ANALYSIS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Anxiety median split
    y= delta12_stai*-1;
    medy = median(y);
    highRiskSubsWeightedAvgVals = allSubsWeightedAvgVals(y > medy); lowRiskSubsWeightedAvgVals = allSubsWeightedAvgVals(y <= medy);
    disp([' Median-Split-Analysis ']);
    [H,P,CI,STATS] = ttest(highRiskSubsWeightedAvgVals);
    disp([currRoiName ': CI>II, High Risk Subs : p=', num2str(round(P,3)), '; T=', num2str(round(STATS.tstat,3)), ';']);
    [H,P,CI,STATS] = ttest(lowRiskSubsWeightedAvgVals);
    disp([currRoiName ': CI>II, Low Risk Subs : p=', num2str(round(P,3)), '; T=', num2str(round(STATS.tstat,3)), ';']);
    [H,P,CI,STATS] = ttest2(highRiskSubsWeightedAvgVals, lowRiskSubsWeightedAvgVals);
    disp([currRoiName ': CI>II, High vs. Low Risk Subs : p=', num2str(round(P,3)), '; T=', num2str(round(STATS.tstat,3)), ';']);
    
    currMeanLowRisk = mean(lowRiskSubsWeightedAvgVals);
    currMeanHighRisk = mean(highRiskSubsWeightedAvgVals);
    currSEMLowRisk = std(lowRiskSubsWeightedAvgVals)/sqrt(length(lowRiskSubsWeightedAvgVals));
    currSEMHighRisk = std(highRiskSubsWeightedAvgVals)/sqrt(length(highRiskSubsWeightedAvgVals));
    
    fh = figure;
    %%% plot it
    bar(1,currMeanLowRisk, 'facecolor', currRoiColor, 'BarWidth', 0.4); hold on
    bar(1.5,currMeanHighRisk, 'facecolor', currRoiColor, 'BarWidth', 0.4); hold on
    errorbar(1,currMeanLowRisk,currSEMLowRisk, 'k');
    errorbar(1.5,currMeanHighRisk,currSEMHighRisk, 'k');
    
    %%% Depression median split
    y= delta12_phq*-1;
    medy = median(y);
    highRiskSubsWeightedAvgVals = allSubsWeightedAvgVals(y > medy); lowRiskSubsWeightedAvgVals = allSubsWeightedAvgVals(y <= medy);
    [H,P,CI,STATS] = ttest(highRiskSubsWeightedAvgVals);
    disp([currRoiName ': CI>II, High Risk Subs : p=', num2str(round(P,3)), '; T=', num2str(round(STATS.tstat,3)), ';']);
    [H,P,CI,STATS] = ttest(lowRiskSubsWeightedAvgVals);
    disp([currRoiName ': CI>II, Low Risk Subs : p=', num2str(round(P,3)), '; T=', num2str(round(STATS.tstat,3)), ';']);
    [H,P,CI,STATS] = ttest2(highRiskSubsWeightedAvgVals, lowRiskSubsWeightedAvgVals);
    disp([currRoiName ': CI>II, High vs. Low Risk Subs : p=', num2str(round(P,3)), '; T=', num2str(round(STATS.tstat,3)), ';']);
    
    currMeanLowRisk = mean(lowRiskSubsWeightedAvgVals);
    currMeanHighRisk = mean(highRiskSubsWeightedAvgVals);
    currSEMLowRisk = std(lowRiskSubsWeightedAvgVals)/sqrt(length(lowRiskSubsWeightedAvgVals));
    currSEMHighRisk = std(highRiskSubsWeightedAvgVals)/sqrt(length(highRiskSubsWeightedAvgVals));
    
    %%% plot it
    bar(2.5,currMeanLowRisk, 'facecolor', currRoiColor, 'BarWidth', 0.4); hold on
    bar(3,currMeanHighRisk, 'facecolor', currRoiColor, 'BarWidth', 0.4); hold on
    errorbar(2.5,currMeanLowRisk,currSEMLowRisk, 'k');
    errorbar(3,currMeanHighRisk,currSEMHighRisk, 'k');
    xlim([.5 3.5]); ylim([-.2 .35]); box off
    set(gca, 'fontsize', 15, 'xticklabel', {})
    axis square
    
    xlabel([currRoiTitle ' (CI>II)'], 'fontsize', axisLabelFontSize, 'FontAngle', 'italic');
    ylabel('Symptoms Changes', 'fontsize', axisLabelFontSize, 'FontAngle', 'italic');
    set(fh, 'position', [1480 300 500 400]);
    
end