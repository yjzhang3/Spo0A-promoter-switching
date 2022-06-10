% get an estimate of how variant the data points are before judging
% if the mean squared error between real and generated are good
% or bad

% every column is a type of mutant and every row is the observation at each
% time point

% for each data points given in the dataset, it is the result of multiple
% observvations/experiments (because each data point has a standard
% deviation and mean)

% there is no point to calculate the average of promoter activity for each
% mutation strain because we know the actvitiy is supposed to be changing
% with time. Each strain is *not* a random variable

% rather, we're interested to see, since each data point is a random variable
% (meaning prmoter activity of a particluar strain at a particular time),
% on average, how much fluctuation is in each data point.

% in other words, we'd like to know the "mean squared error" of all data
% points, where the "error" is exactly the standard deviation (squared in 
% the formula. Then we will divide this standard deviation for each RV by its mean, sum them, and
% take the average. 

% the entire dataset is ONE random variable and each data point is an
% ovsrevton. Solve for the weighted variance of thisshs


load('promoter_activity_single.mat')

real_data_Ps = Ps_promoter_activity_mean(:,2:9); % exclude the real WT, instead WT is Ps only, no mutations
real_data_Ps_std = Ps_promoter_activity_std(:,2:9);
real_data_Pv = Pv_promoter_activity_mean(:,2:9);
real_data_Pv_std = Pv_promoter_activity_std(:,2:9);

sumd = 0;
for ii = 1:length(Ps_promoter_activity_mean(:,1)) % number of rows
    for jj = 1:length(Ps_promoter_activity_mean(1,:)) % number of columns
        m = Ps_promoter_activity_mean(ii,jj);
        s = Ps_promoter_activity_std(ii,jj);
        d = s^2/m^2; % weighted 'mean square error of each obseerevation"
        sumd = sumd + d;
    end
end
Ps_overall_var = sumd/(8*9); % then take the average

clear sumd ii jj d

sumd = 0;
for ii = 1:length(Pv_promoter_activity_mean(:,1)) % number of rows
    for jj = 1:length(Pv_promoter_activity_mean(1,:)) % number of columns
        m = Pv_promoter_activity_mean(ii,jj);
        s = Pv_promoter_activity_std(ii,jj);
        d = s^2/m^2; % weighted 'mean square error of each obseerevation"
        sumd = sumd + d;
    end
end
Pv_overall_var = sumd/(8*9); % then take the average


