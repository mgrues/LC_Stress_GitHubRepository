%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MG: START HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%%% PATHS:
addpath('C:\Users\mgrue\PROJECTS\RARS_SNS_EYE\TOOLS\analysisScriptEyetrackingHelenaMarcB');

%%% OPTIONS:
pairedTTest = 1; %%% 1 = you need to thave two data matrices and take the difference of them means first %%%
%%% Which Filtering
doNoFiltering=0;
doLowPass = 0;
doBandPass = 1;
doMoveAvg = 1;
moveAvgKernel = 500; %%% this is in miliseconds, 500ms is used %%%

%%% some specific Silvia script variables : outsource them later
figures = [228]; %225 216 217 218     %[212 213 216 217 218 224 225]
numRounds = length(figures);

%%% essential variables
create_matrices_CBPT = 1;
create_figures = 1; % switch generates figures

%%% THRESHOLDS:
cluster_forming_threshold = 3; % this denotes a T-stat
cluster_forming_threshold_neg = -cluster_forming_threshold; % to calculate a negative contrast

%%% ITERATIONS:
numIterations = 1000;
%%% Bonferroni correction variables
alpha_level = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GET THE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('C:\Users\mgrue\PROJECTS\LC_Stress')
if pairedTTest
    %%% two data sets
    %%% tic
    disp([' --- Loading data ---'])
    
    if doNoFiltering
        %%% DATASET 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load singlePupilTrialsStimLocked_CI
        data_to_analyze_ons = singlePupilTrialsStimLocked_CI;
        %%% DATASET 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load singlePupilTrialsStimLocked_II
        data_to_analyze_ons2 = singlePupilTrialsStimLocked_II;
    end
    
    if doLowPass
        %%% DATASET 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load singlePupilTrialsStimLocked_CI_LowPass
        data_to_analyze_ons = singlePupilTrialsStimLocked_CI;
        %%% DATASET 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load singlePupilTrialsStimLocked_II_LowPass
        data_to_analyze_ons2 = singlePupilTrialsStimLocked_II;
    end
    
    if doBandPass
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% CI>II
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% DATASET 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load singlePupilTrialsStimLocked_CI_BandPass
        data_to_analyze_ons = singlePupilTrialsStimLocked_CI;
        %%% DATASET 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load singlePupilTrialsStimLocked_II_BandPass
        data_to_analyze_ons2 = singlePupilTrialsStimLocked_II;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% CC>IC
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% DATASET 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load singlePupilTrialsStimLocked_CC_BandPass
        data_to_analyze_ons = singlePupilTrialsStimLocked_CC;
        %%% DATASET 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load singlePupilTrialsStimLocked_IC_BandPass
        data_to_analyze_ons2 = singlePupilTrialsStimLocked_IC;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% I>C
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% DATASET 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load singlePupilTrialsStimLocked_II_BandPass
        load singlePupilTrialsStimLocked_CI_BandPass
        data_to_analyze_ons = cat(1, singlePupilTrialsStimLocked_II, singlePupilTrialsStimLocked_CI);
        %%% DATASET 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load singlePupilTrialsStimLocked_IC_BandPass
        load singlePupilTrialsStimLocked_CC_BandPass
        data_to_analyze_ons2 = cat(1, singlePupilTrialsStimLocked_IC, singlePupilTrialsStimLocked_CC);
        
    end
else
    %%% just one data set
    disp(['Loading data'])
    load data_to_analyze_ons_example
    %%% data_to_analyze_ons = data_to_analyze_ons;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set random number generator seed to computer clock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RandStream.setGlobalStream ...
    (RandStream('mt19937ar','seed',sum(100*clock)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% average matrices first across trial types so you can average across subjects
% in the plots
for j = 1:length(figures)
    % get true within-subject timeseries average
    
    data_to_analyze_sjmean = nan(size(data_to_analyze_ons,3), size(data_to_analyze_ons,2));
    for sub = 1:size(data_to_analyze_ons,3)
        
        % index which rows to average over
        % for a 1-sample t-test: just take the mean per
        % subject over all trials that belong to the condition
        within_sj_cond_rows = (1:size(data_to_analyze_ons,1))';
        
        if pairedTTest
            
            %%% DO SUBSTRACTION HERE:
            tmp1 = nanmean(data_to_analyze_ons(within_sj_cond_rows,:,sub));
            tmp2 = nanmean(data_to_analyze_ons2(within_sj_cond_rows,:,sub));
            if doMoveAvg
                data_to_analyze_sjmean(sub,:) = movmean(tmp1, moveAvgKernel, 2) - movmean(tmp2, moveAvgKernel, 2); %mean for this trial type on the subject level
            else
                data_to_analyze_sjmean(sub,:) = tmp1 - tmp2; %mean for this trial type on the subject level
            end
            %%% keyboard
            %%% figure; plot(tmp1); hold on; plot(tmp2,'r' ); grid on
            %%% figure; plot(nanmean(data_to_analyze_sjmean)); grid on
            
        else
            % average within subject over all trials of one condition -> get one
            % averaged time series for this participant and trial type
            if doMoveAvg
                tmp1 = nanmean(data_to_analyze_ons(within_sj_cond_rows,:,sub));
                data_to_analyze_sjmean(sub,:) = movmean(tmp1, moveAvgKernel, 2); %mean for this trial type on the subject level
            else
                data_to_analyze_sjmean(sub,:) = nanmean(data_to_analyze_ons(within_sj_cond_rows,:,sub)); %mean for this trial type on the subject level
            end
        end
        %%% VISUALIZE THE EFFECT OF SMOOTHING
        % % % figure;
        % % % plot(mean(allSubsPupilStimLocked_II, 1),'b'); hold on
        % % % a = movmean(allSubsPupilStimLocked_II, 500, 2);
        % % % plot(mean(a, 1),'g'); hold on
    end
    
    %%% keyboard
    %%% figure; plot(nanmean(data_to_analyze_sjmean))
    
    % get paired t-stat across time points for each round
    ons_t_stat_CFT_mat = nan(numRounds,size(data_to_analyze_sjmean,2));
    ons_t_stat_CFT_mat_neg = nan(numRounds,size(data_to_analyze_sjmean,2));
    
    for r = 1:numRounds
        
        %one sample t-test
        [h,p_ons{r},ci,stats_ons] = ttest(data_to_analyze_sjmean);
        %%% figure; plot(stats_ons.tstat); grid on;
        %%% figure; plot(p_ons{r}); grid on;
        ons_t_stat_CFT_mat(r,:) = stats_ons.tstat > cluster_forming_threshold; 
        %%% figure; plot(ons_t_stat_CFT_mat(r,:)); ylim([-.1 1.1]) %%%
        ons_t_stat_CFT_mat_neg(r,:) = stats_ons.tstat < cluster_forming_threshold_neg;
        %%% figure; plot(ons_t_stat_CFT_mat_neg(r,:)); ylim([-.1 1.1]) %%%
        
        % find starting and ending indices and cluster sizes of significant
        % clusters IN ORIGINAL DATA
        ons_dsig = diff([0 ons_t_stat_CFT_mat(r,:) 0]); %tells transition point where is change from one time to next -> 0 to 1
        ons_cluster_start_idxs{r} = find(ons_dsig > 0); % check which transitions are positive: start point
        ons_cluster_end_idxs{r} = find(ons_dsig < 0)-1; % where is the end point (diff is negative) -> then go one back bc was "diff" (subtraction)
        ons_cluster_size{r} = ons_cluster_end_idxs{r}-ons_cluster_start_idxs{r}+1; %size of the cluster
        
        %negative direction
        %we are operating on a logical, so it stays the same as above
        ons_dsig_neg = diff([0 ons_t_stat_CFT_mat_neg(r,:) 0]); %tells transition point where is change from one time to next -> 0 to 1
        ons_cluster_start_idxs_neg{r} = find(ons_dsig_neg > 0); % check which transitions are positive: start point
        ons_cluster_end_idxs_neg{r} = find(ons_dsig_neg < 0)-1; % where is the end point (diff is negative) -> then go one back bc was "diff" (subtraction)
        ons_cluster_size_neg{r} = ons_cluster_end_idxs_neg{r}-ons_cluster_start_idxs_neg{r}+1; %size of the cluster
        
        
    end
    
    
    %% PERMUTED DATA
    % cluster-based permutation paired t-test
    % get permuted within-subject timeseries averages
    
    if create_matrices_CBPT
        % pre-allocate
        ons_null_dist = nan(numRounds,numIterations);
        
        for i = 1:numIterations
            
            disp(['Working on permutation: ', num2str(i)])
            
            % pre-allocate
            data_to_analyze_sjmean_CBPT = nan(size(data_to_analyze_sjmean));
            
            for sub = 1:size(data_to_analyze_ons,3) % all subs
                
                % in the case of the 1-sample t-test against zero:
                % randomly assign a vector of 1 and -1 that has the
                % length of the number of trials that we want to
                % permute and shuffle it; then average whole time
                % series over all trials, values go into average
                % as positive in case of a 1 and as negative in
                % case of a -1 index in that trial
                
                %%% SWAPPING SIGNS
                if pairedTTest
                    
                    %%% Swap signs in data from first condition
                    for p = 1:size(data_to_analyze_ons,1) %all trials
                        perm_idx(p) = rand>0.5;
                    end
                    perm_idx = double(perm_idx);
                    perm_idx(perm_idx==0)=-1;
                    % swap around the signs in half of the columns
                    rand_data_to_analyze = bsxfun(@times,data_to_analyze_ons(within_sj_cond_rows,:,sub),perm_idx'); %edit SUM 11.01.18: this was before data_to_analyze, but that does not make sense, because the data we want to compare to are cut and just 7 sec long!
                    
                    %%% Swap signs in data from second condition
                    for p = 1:size(data_to_analyze_ons2,1) %all trials
                        perm_idx(p) = rand>0.5;
                    end
                    perm_idx = double(perm_idx);
                    perm_idx(perm_idx==0)=-1;
                    % swap around the signs in half of the columns
                    rand_data_to_analyze2 = bsxfun(@times,data_to_analyze_ons2(within_sj_cond_rows,:,sub),perm_idx'); %edit SUM 11.01.18: this was before data_to_analyze, but that does not make sense, because the data we want to compare to are cut and just 7 sec long!
                else
                    for p = 1:size(data_to_analyze_ons,1) %all trials
                        perm_idx(p) = rand>0.5;
                    end
                    perm_idx = double(perm_idx);
                    perm_idx(perm_idx==0)=-1;
                    % swap around the signs in half of the columns
                    rand_data_to_analyze = bsxfun(@times,data_to_analyze_ons(within_sj_cond_rows,:,sub),perm_idx'); %edit SUM 11.01.18: this was before data_to_analyze, but that does not make sense, because the data we want to compare to are cut and just 7 sec long!
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% THEN AVERAGE THE PERMUTED DATA
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% then average
                %%% data_to_analyze_sjmean_CBPT(sub,:) = nanmean(rand_data_to_analyze);
                
                if pairedTTest
                    %%% DO SUBSTRACTION HERE:
                    tmp1 = nanmean(rand_data_to_analyze);
                    tmp2 = nanmean(rand_data_to_analyze2);
                    if doMoveAvg
                        data_to_analyze_sjmean_CBPT(sub,:) = movmean(tmp1, moveAvgKernel, 2) - movmean(tmp2, moveAvgKernel, 2); %mean for this trial type on the subject level
                    else
                        data_to_analyze_sjmean_CBPT(sub,:) = tmp1 - tmp2; %mean for this trial type on the subject level
                    end
                    %%% figure; plot(tmp1); hold on; plot(tmp2,'r' )
                    %%% figure; plot(nanmean(data_to_analyze_sjmean))
                    
                else
                    % average within subject over all trials of one condition -> get one
                    % averaged time series for this participant and trial type
                    if doMoveAvg
                        tmp1 = nanmean(rand_data_to_analyze);
                        data_to_analyze_sjmean_CBPT(sub,:) = movmean(tmp1, moveAvgKernel, 2); %mean for this trial type on the subject level
                    else
                        data_to_analyze_sjmean_CBPT(sub,:) = nanmean(rand_data_to_analyze(within_sj_cond_rows,:,sub)); %mean for this trial type on the subject level
                    end
                end
                
                % in case of a paired t-test where we want to
                % compare conditions, eg row 1:20 are condition 1
                % and row 21:40 are condition 2, we calculate a
                % paired t-test of row 1:20 vs row 21:40 to
                % determine the true data outcome, and in order to
                % create the null distribution, we shuffle the row
                % labels for all rows 1:40 and then take the rows
                % 1:20 of the permuted sample to test against rows
                % 21:40 of the permuted sample
                
            end
            
            % get paired t-stat across time points for each round
            % pre-allocate
            ons_t_stat_CFT_mat_CBPT = nan(numRounds,size(data_to_analyze_sjmean_CBPT,2));
            
            for r = 1:numRounds
                
                [h,p,ci,stats_ons] = ttest(data_to_analyze_sjmean_CBPT);
                ons_t_stat_CFT_mat_CBPT(r,:) = stats_ons.tstat > cluster_forming_threshold;
                ons_t_stat_CFT_mat_CBPT_neg(r,:) = stats_ons.tstat < cluster_forming_threshold_neg;
                
                % find starting and ending indices and cluster sizes
                ons_dsig = diff([0 ons_t_stat_CFT_mat_CBPT(r,:) 0]);
                ons_cluster_start_idxs_CBPT = find(ons_dsig > 0);
                ons_cluster_end_idxs_CBPT = find(ons_dsig < 0)-1;
                ons_cluster_size_CBPT = ons_cluster_end_idxs_CBPT-ons_cluster_start_idxs_CBPT+1;
                
                % negative direction
                % find starting and ending indices and cluster sizes
                ons_dsig_neg = diff([0 ons_t_stat_CFT_mat_CBPT_neg(r,:) 0]);
                ons_cluster_start_idxs_CBPT_neg = find(ons_dsig_neg > 0);
                ons_cluster_end_idxs_CBPT_neg = find(ons_dsig_neg < 0)-1;
                ons_cluster_size_CBPT_neg = ons_cluster_end_idxs_CBPT_neg-ons_cluster_start_idxs_CBPT_neg+1;
                
                % only store size of largest significant cluster found in
                % random data:
                if max(ons_cluster_size_CBPT) > 0
                    ons_null_dist(r,i) = max(ons_cluster_size_CBPT);
                else
                    ons_null_dist(r,i) = 0;
                end
                
                if max(ons_cluster_size_CBPT_neg) > 0
                    ons_null_dist_neg(r,i) = max(ons_cluster_size_CBPT_neg);
                else
                    ons_null_dist_neg(r,i) = 0;
                end
                
            end
            
        end % iterations loop - with 5000 iterations and 16000 samples (16 sec) took 11 minutes to compute
        % in Anjali's data: 1.9 seconds per iteration; ok for 1,000 iterations but maybe optimize for 10,000 iterations
        
        % store variables
        O.ons_null_dist = ons_null_dist; % permuted clusters null distribution
        O.ons_null_dist_neg = ons_null_dist_neg; %permuted clusters null distribution - negative direction
        O.ons_cluster_size = ons_cluster_size; % true clusters
        O.ons_cluster_size_neg = ons_cluster_size_neg; % true clusters - negative direction
        O.cluster_start_idxs = ons_cluster_start_idxs;
        O.cluster_start_idxs_neg = ons_cluster_start_idxs_neg;
        O.cluster_end_idxs = ons_cluster_end_idxs;
        O.cluster_end_idxs_neg = ons_cluster_end_idxs_neg;
        O.t_stat_CFT_mat= ons_t_stat_CFT_mat;
        O.t_stat_CFT_mat_neg= ons_t_stat_CFT_mat_neg;
        O.ons_pval = p_ons;
        
        %%% save([num2str(currfig) '_ons_null_dist_cluster_size_CFT_' num2str(cluster_forming_threshold) '_' num2str(numIterations) '_iters.mat'],'-struct','O');
    end
    
    disp([num2str(numIterations), ' permutations done after: '])
    toc
    
    %%% keyboard
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Figures & PLOTTING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Conditions = {'Regulate Pos', 'Regulate Neg', 'View Pos' ,'View Neg', 'View Neut'};
    Colors = {'g','r', 'b', 'y', 'k'};
    axis_fontsize = 18;
    
    if create_figures
        
        %load([num2str(currfig) '_ons_null_dist_cluster_size_CFT_' num2str(cluster_forming_threshold) '_' num2str(numIterations) '_iters.mat']);
        
        % INFERENCE:
        % For each true cluster, ask how many of the maximal permuted clusters in the
        % null distribution were larger than it to get p-value
        % (so if 14 continuous timebins had a t-stat greater than cluster-forming
        % threshold and there were 26 clusters in the null distribution > 14,
        % then p-value = 0.026 for 1000 iterations)
        
        % for each cluster: instead of checking a t-distribution, I want to
        % know how likely it is that I get a cluster of size 14 just by chance;
        % -> if my null distribution has 26 clusters that are greater or equal
        % than 14 -> p = 26/numIterations
        
        % -> only take largest cluster: to make sure that you can't have a
        % value above that ( this accounts for multiple comparison problem: I
        % just make 1 comparison: take the largest cluster; most conservative,
        % not like taking mean cluster size Nichols & Holmes 2002)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% POSITIVE DIRECTION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ons_perm_p = {};
        % get permuted p vals
        for r = 1:numRounds
            for ons_cluster = 1:size(ons_cluster_size{r},2)
                % look for all clusters in the whole trial time
                ons_perm_p{r}(ons_cluster,:) = sum(ons_null_dist(r,:) > ons_cluster_size{r}(ons_cluster))/numIterations;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% NEGATIVE DIRECTION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ons_perm_p_neg = {};
        % get permuted p vals
        for r = 1:numRounds
            for ons_cluster = 1:size(ons_cluster_size_neg{r},2)
                % look for all clusters in the whole trial time
                ons_perm_p_neg{r}(ons_cluster,:) = sum(ons_null_dist_neg(r,:) > ons_cluster_size_neg{r}(ons_cluster))/numIterations;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% DATA FIGURE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FOR CONTRAST IN POSITVE DIRECTION (greater than zero)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(ons_cluster_start_idxs{:})
            
            %%% xes = [0:7000]; %%% original line%%%
            xes = [0:size(data_to_analyze_sjmean,2)-1];
            clear expl;
            % shaded errorbar plot (onset-locked)
            fh1 = figure;
            for r = 1:numRounds
                avg = mean(data_to_analyze_sjmean)';
                sem = squeeze(std(data_to_analyze_sjmean)./sqrt(size(data_to_analyze_ons,3)))';
                expl = shadedErrorBar([xes]/1000,avg,sem, {'-', 'color', [0 0 0], 'lineWidth', 2});
                clear sem;
                hold on;
                
                ylabel([' Pupil Diameter'],'fontsize',axis_fontsize);
                set(gca,'fontsize',axis_fontsize);
                xlabel('Time (seconds)','fontsize',axis_fontsize);
                box("off")
                
                cluster_xes = [];
                midpt_cluster_xes = [];
                
                for ons_cluster = 1:size(ons_cluster_size{r},2)
                    
                    cluster_xes = cat(2,cluster_xes,xes(ons_cluster_start_idxs{r}(ons_cluster):ons_cluster_end_idxs{r}(ons_cluster)));
                    midpt_cluster_xes = cat(2,midpt_cluster_xes,round(mean(xes(ons_cluster_start_idxs{r}(ons_cluster):ons_cluster_end_idxs{r}(ons_cluster)))));
                    
                    if ons_perm_p{r}(ons_cluster,:) < 0.001
                        pval_text{ons_cluster} = '***';
                    elseif ons_perm_p{r}(ons_cluster,:) < 0.01
                        pval_text{ons_cluster} = '**';
                    elseif ons_perm_p{r}(ons_cluster,:) < 0.05
                        pval_text{ons_cluster} = '*';
                    end
                    
                end
                
                %pvals_CBPT = text(cluster_xes/1000,repmat(max(ylim)-0.05,1,size(cluster_xes,2)),'bold'); %,'Color','blue'
                
                %plot line for the significant area:
                hold on
                line([ons_cluster_start_idxs{r}(ons_cluster)/1000, ons_cluster_end_idxs{r}(ons_cluster)/1000], [max(ylim)-0.05,max(ylim)-0.05], 'linewidth', 4, 'Color','black');
                
                for ons_cluster = 1:size(ons_cluster_size{r},2)
                    try
                        stars_CBPT = text((midpt_cluster_xes(ons_cluster)/1000)-0.2,repmat(max(ylim)-0.045,1,1),pval_text{ons_cluster},'Fontsize',32, 'FontWeight', 'bold'); %,'Color','blue'
                    catch
                        disp('Caught !!! Check it!!!')
                    end
                end
                
            end
            set(fh1,'Position',[100 100 700 500]);
            ah1 = get(fh1, 'CurrentAxes');
            set(ah1,'LineWidth',2)
            set(ah1,'TickLength',[0,0])
        end %%% Negative direction ends here %%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FOR CONTRAST IN NEGATIVE DIRECTION (smaller than zero)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(ons_cluster_start_idxs_neg{:})
            
            %%% xes = [0:7000];
            xes = [0:size(data_to_analyze_sjmean,2)-1];
            clear expl;
            % shaded errorbar plot (onset-locked)
            fh1 = figure(1);
            for r = 1:numRounds
                avg = mean(data_to_analyze_sjmean)';
                sem = squeeze(std(data_to_analyze_sjmean)./sqrt(size(data_to_analyze_ons,3)))';
                expl = shadedErrorBar([xes]/1000,avg,sem);
                clear sem;
                hold on;
                
                ylabel([' Pupil Diameter'],'fontsize',axis_fontsize);
                set(gca,'fontsize',axis_fontsize);
                xlabel('Time (s)','fontsize',axis_fontsize);
                xlim([6 15])
                box("off")
                grid on;
                
                cluster_xes = [];
                midpt_cluster_xes = [];
                
                for ons_cluster = 1:size(ons_cluster_size_neg{r},2)
                    
                    cluster_xes = cat(2,cluster_xes,xes(ons_cluster_start_idxs_neg{r}(ons_cluster):ons_cluster_end_idxs_neg{r}(ons_cluster)));
                    midpt_cluster_xes = cat(2,midpt_cluster_xes,round(median(xes(ons_cluster_start_idxs_neg{r}(ons_cluster):ons_cluster_end_idxs_neg{r}(ons_cluster)))));
                    
                    if ons_perm_p_neg{r}(ons_cluster,:) < 0.001
                        pval_text{ons_cluster} = '***';
                    elseif ons_perm_p_neg{r}(ons_cluster,:) < 0.01
                        pval_text{ons_cluster} = '**';
                    elseif ons_perm_p_neg{r}(ons_cluster,:) < 0.05
                        pval_text{ons_cluster} = '*';
                    end
                    
                end
                
                %pvals_CBPT = text(cluster_xes/1000,repmat(max(ylim)-0.05,1,size(cluster_xes,2))); %,'Color','blue'
                line([ons_cluster_start_idxs_neg{r}(ons_cluster)/1000,ons_cluster_end_idxs_neg{r}(ons_cluster)/1000],[max(ylim)-0.05,max(ylim)-0.05], 'linewidth', 4, 'Color','black');
                
                for ons_cluster = 1:size(ons_cluster_size_neg{r},2)
                    try
                        stars_CBPT = text(midpt_cluster_xes(ons_cluster)/1000,repmat(max(ylim)-0.049,1,1),pval_text{ons_cluster}, 'Fontsize',32); %,'Color','blue','Fontsize',32
                    catch
                        disp('Caught !!! Check it!!!')
                    end
                end
                
                %             bonf_xes = xes(find(ons_bonf_p_neg{r}));
                %             stars_bonf = text(bonf_xes/1000,repmat(max(ylim)-0.06,1,size(bonf_xes,2)),'*','fontsize',6);
                %
            end
            
        end %%% Negative direction ends here %%%
        
        axis square
        
    end % switch for base stim figures no group info
end %switch for all contrasts we want to examine
toc

keyboard