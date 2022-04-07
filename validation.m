%% validation file
% This file is just for validation. We first generate "real data" using our
% thermodynamic algorithm, then feed this into the particle swarm
% optimization model to see if we can recover the same energy parameter
% values that are used to generate the "real data".

%% generate enough data points
TF_conc_fun = @(x) x^3+300*x;
tspan = 10;
TF_conc_t = zeros(tspan,1);

for rr = 1:length(TF_conc_t)
    TF_conc_t(rr) = TF_conc_fun(rr-1);
end


%% parameters
nbd = 4;
energyi = [10,10,10,10,-2,-2,10,3,-5,-10];
mut_mat = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[0,0,0]];
RNAp_conc = 500;

real_data = zeros(length(mut_mat),length(TF_conc_t));

lb = zeros(10,1)-10;
% lb(5:6) = -0.0001;
% lb(8) = -0.0001;
ub = zeros(10,1)+10;
% ub(5:6) = 0.0001;
% ub(8) = 0.0001;

%% generate real data
for mm = 1:length(mut_mat)
    subplot(4,2,mm)
    TR_t = time_dep_TR_new(nbd,energyi,TF_conc_t,RNAp_conc,mut_mat(mm,:));
    real_data(mm,:) = TR_t;
    plot(TF_conc_t,TR_t,'LineWidth',4)
    ylim([0 1])
end

%% use particle swarm to recover the energy values
x = fit_data_new(nbd,TF_conc_t,RNAp_conc,mut_mat,real_data,lb,ub);

%% now plot dynamics resulted from each parameter set
figure();
for kk = 1:length(mut_mat)
    subplot(4,2,kk)
    plot(TF_conc_t,real_data(kk,:),'LineWidth',4)
    ylim([0 1])
    xlabel('TF concentration')
    ylabel('transcription rate')
end

figure();
for yy = 1:length(mut_mat)
    subplot(4,2,yy)
    TR = time_dep_TR_new(nbd,x(8,:),TF_conc_t,RNAp_conc,mut_mat(yy,:));
    plot(TF_conc_t,TR,'LineWidth',4)
    ylim([0 1])
    xlabel('TF concentration')
    ylabel('transcription rate')
end

%% preicions, recall, and F1 score
error = zeros(10,1); %  error for each iteration
err_m = mean(x); % mean of each column, which is also each type of energy
err_std = std(x);

for ii = 1:length(x)
    err = sum(((x(ii,:)-energyi)./x(ii,:)).^2)/length(x);
    error(ii) = err;
end

    