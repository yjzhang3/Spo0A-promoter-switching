%% a note on raw data
% Ps: column 2 and onward are one promoter Ps with all mutations. First
% column is true wild type
% Pv: column 2-9 are one promoter Pv with all mutations. First column is
% true wild type. Last column is when both promoters are deleted.
% all data start from T2 to T10 (index 1 to 9)
%% parameters
real_data_Ps = Ps_promoter_activity_mean(2:8,1:8); % exclude the real WT, instead WT is 

nbd = 4;
TF_conc_t = wta(1:5:end);
RNAp_conc = 500;
mut_mat = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[0,0,1];[0,0,0]];
lb = zeros(10,1)-10;
ub = zeros(10,1)+10;

%%
pars = fit_data(nbd,TF_conc_t,RNAp_conc,mut_mat,real_data_Ps,lb,ub);

%% figure();
for kk = 1:length(mut_mat)
    subplot(4,2,kk)
    plot(TF_conc_t,real_data(:,kk),'LineWidth',4)
    ylim([0 1])
    xlabel('TF concentration')
    ylabel('transcription rate')
    title(sprintf('mutation %d',kk))
end

figure();
for yy = 1:length(mut_mat)
    subplot(4,2,yy)
    TR = time_dep_TR_new(nbd,pars(1,:),TF_conc_t,RNAp_conc,mut_mat(yy,:));
    plot(TF_conc_t,TR,'LineWidth',4)
    ylim([0 1])
    xlabel('TF concentration')
    ylabel('transcription rate')
    title(sprintf('mutation %d',kk))
end