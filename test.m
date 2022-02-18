%% test if protein protein interaction helps change transcription rate
energyi_normal = [5,8,10,0,0,0]; % no pp interaction
energyi_actA = [5,8,10,0,-1.7,0]; % TF of site A and mRNA attract
energyi_repB = [5,8,10,0,0,1.5]; % TF of site B and mRNA repel
energyi_AB_actAB = [5,8,10,-1.9,-4,-6.9]; % TF attract, and they each attract mRNA
mRNA_conc = 1000;
TF_conc = 1000;
nbd = 3;

TR_normal = transcription_rate(nbd,energyi_normal,TF_conc,mRNA_conc);
TR_actA = transcription_rate(nbd,energyi_actA,TF_conc,mRNA_conc);
TR_repB = transcription_rate(nbd,energyi_repB,TF_conc,mRNA_conc);
TR_AB_actAB = transcription_rate(nbd,energyi_AB_actAB,TF_conc,mRNA_conc);

%% generate time dependent TR (fake data) to get an intuition about real data
conc_t = @(t) 2*t;
tspan = 10;
TR_overall = zeros(tspan,7);

energyi_normal = [5,8,10,0,0,0]; % no pp interaction
energyi_actA = [5,8,10,0,0,-6.9]; % only site A TF attracts mRNA
energyi_actA_actB = [5,8,10,0,-4,-6.9]; % both TF attract mRNA
energyi_AB_actA_actB = [5,8,10,-1.9,-4,-6.9]; % TF attract each other, and they each attract mRNA

energyi_repB = [5,8,10,0,3,0]; % TF of site B and mRNA repel
energyi_repA_repB = [5,8,10,0,3,4]; % TF of site A and B both repel mRNA
energyi_AB_repA_repB = [5,8,10,1.9,4,6.9]; % TF repel each other, and they both repel mRNA

mRNA_conc = 1000;
nbd = 3;

for tt = 1:tspan
    TR_overall(tt,1) = transcription_rate(nbd,energyi_normal,conc_t(tt),mRNA_conc);
    TR_overall(tt,2) = transcription_rate(nbd,energyi_actA,conc_t(tt),mRNA_conc);
    TR_overall(tt,3) = transcription_rate(nbd,energyi_actA_actB,conc_t(tt),mRNA_conc);
    TR_overall(tt,4) = transcription_rate(nbd,energyi_AB_actA_actB,conc_t(tt),mRNA_conc);
    TR_overall(tt,5) = transcription_rate(nbd,energyi_repB,conc_t(tt),mRNA_conc);
    TR_overall(tt,6) = transcription_rate(nbd,energyi_repA_repB,conc_t(tt),mRNA_conc);
    TR_overall(tt,7) = transcription_rate(nbd,energyi_AB_repA_repB,conc_t(tt),mRNA_conc);
end

figure();
subplot(2,1,1)
plot(1:1:tspan,TR_overall(:,1))
hold on
plot(1:1:tspan,TR_overall(:,2),'-+')
hold on
plot(1:1:tspan,TR_overall(:,3),'-+')
hold on
plot(1:1:tspan,TR_overall(:,4),'-+')
xlabel('time')
ylabel('transcription rate')
set(gca,'FontSize',17)
legend('normal','A-mRNA','A-mRNA and B-mRNA',"A-B, A-mRNA and B-mRNA")
spec = sprintf('Gene regulation with %d binding sites and 1 promoter',nbd-1);
title(spec)

subplot(2,1,2)
plot(1:1:tspan,TR_overall(:,1))
hold on
plot(1:1:tspan,TR_overall(:,5),'-o')
hold on
plot(1:1:tspan,TR_overall(:,6),'-o')
hold on
plot(1:1:tspan,TR_overall(:,7),'-o')
xlabel('time')
ylabel('transcription rate')
set(gca,'FontSize',17)
legend('normal','B|mRNA','A|mRNA and B|mRNA',"A|B, A|mRNA and B|mRNA")
spec = sprintf('Gene regulation with %d binding sites and 1 promoter',nbd-1);
title(spec)
