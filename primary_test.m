%% test if protein protein interaction helps change transcription rate
energyi_normal = [5,5,10,0,0,0]; % no pp interaction
energyi_actA = [5,5,10,0,-1.7,0]; % TF of site A and mRNA attract
energyi_repB = [5,5,10,0,0,1.5]; % TF of site B and mRNA repel
energyi_AB_actAB = [5,5,10,-1.9,-4,-6.9]; % TF attract, and they each attract mRNA
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

energyi_normal = [5,5,10,0,0,0]; % no pp interaction
energyi_actA = [5,5,10,0,0,-6.9]; % only site A TF attracts mRNA
energyi_actA_actB = [5,5,10,0,-4,-6.9]; % both TF attract mRNA
energyi_AB_actA_actB = [5,5,10,-1.9,-4,-6.9]; % TF attract each other, and they each attract mRNA

energyi_repB = [5,5,10,0,0,6.9]; % TF of site B and mRNA repel
energyi_repA_repB = [5,5,10,0,4,6.9]; % TF of site A and B both repel mRNA
energyi_AB_repA_repB = [5,5,10,1.9,4,6.9]; % TF repel each other, and they both repel mRNA

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
plot(1:1:tspan,TR_overall(:,1),'-+','LineWidth',3)
hold on
plot(1:1:tspan,TR_overall(:,2),'-+','LineWidth',3)
hold on
plot(1:1:tspan,TR_overall(:,3),'-+','LineWidth',3)
hold on
plot(1:1:tspan,TR_overall(:,4),'-+','LineWidth',3)
xlabel('time')
ylabel('transcription rate')
set(gca,'FontSize',17)
legend('normal','A-mRNA','A-mRNA and B-mRNA',"A-B, A-mRNA and B-mRNA")
spec = sprintf('Gene regulation with %d binding sites and 1 promoter',nbd-1);
title(spec)

subplot(2,1,2)
plot(1:1:tspan,TR_overall(:,1),'-o','LineWidth',3)
hold on
plot(1:1:tspan,TR_overall(:,5),'-o','LineWidth',3)
hold on
plot(1:1:tspan,TR_overall(:,6),'-o','LineWidth',3)
hold on
plot(1:1:tspan,TR_overall(:,7),'-o','LineWidth',3)
xlabel('time')
ylabel('transcription rate')
set(gca,'FontSize',17)
legend('normal','B|mRNA','A|mRNA and B|mRNA',"A|B, A|mRNA and B|mRNA")
spec = sprintf('Gene regulation with %d binding sites and 1 promoter',nbd-1);
title(spec)

%% explore more complex dynamics

step = 2; % how many steps we want to test up and down
% then we would generate 2*step more energy arrays
% total there is 2*step+1 arrays

pos = 5; % which energy value are you tuning?

TF_conc_t = @(t) 2*t;
tspan = 10;
TR_overall = zeros(tspan,2*step+1);
nbd = 3;
mRNA_conc = 10000000;

% A activates mRNA, B repels mRNA
energyi_normal = [5,5,10,0,-3,3];
energyi_variation = zeros(length(energyi_normal),2*step+1);
energyi_variation(:,end) = energyi_normal;


for ii = 1:step
    currup = [5,5,10,0,energyi_normal(pos)+ii,3];
    curr2 = [5,5,10,0,energyi_normal(pos)-ii,3];
    energyi_variation(:,2*ii-1) = currup;
    energyi_variation(:,2*ii) = curr2;
end

for ss = 1:2*step+1
    for tt = 1:tspan
        TR_overall(tt,ss) = transcription_rate(nbd,energyi_variation(:,ss),TF_conc_t(tt),mRNA_conc);
    end
    plot(1:tspan,TR_overall(:,ss),'LineWidth',3)
    spec = sprintf('A-mRNA = %d',energyi_variation(5,ss)); 
    yline(TR_overall(end-1,ss),'--',spec)
    xlabel('time')
    ylabel('transcription rate')
    title('tuning TF-mRNA attraction')
    hold on
    set(gca,'FontSize',10)
end


%% explore more complex dynamics (2 positions simultaneously)

step = 2; % how many steps we want to test up and down
% then we would generate 2*step more energy arrays
% total there is 2*step+1 arrays
stepsize = 1.5; % how much energy per step are you tuning?

pos_all = [5,6]; % which energy value are you tuning?

TF_conc_t = @(t) 2*t;
tspan = 10;
nbd = 3;
mRNA_conc = 1000;

% A activates mRNA, B repels mRNA
original = [0.5,0.8,20,0,-5,5];
% energyi_variation = energyi_normal,power(2*step,num_pos)+1);
% energyi_variation(:,end) = energyi_normal;


all_combo = nmultichoosek([0-step:0+step],length(pos_all));


for ss = 1:length(all_combo)
    energyi_normal = original;
    for nn = 1:length(pos_all)
        energyi_variation = change_ele(energyi_normal,pos_all(nn),all_combo(ss,nn).*stepsize);
        energyi_normal = energyi_variation;
    end
    
    for tt = 1:tspan
        curr = transcription_rate(nbd,energyi_variation,TF_conc_t(tt),mRNA_conc);
        TR_t_curr(tt) = curr;
    end
    
    plot(1:tspan,TR_t_curr)
    spec1 = sprintf('A-mRNA = %.2f and B-mRNA = %.2f',energyi_variation(pos_all(1)),energyi_variation(pos_all(2))); 
    yline(TR_t_curr(tt-1),'--',spec1)
    xlabel('time')
    ylabel('transcription rate')
    spec2 = sprintf('tuning TF-mRNA attraction (%d sites simultaneously)',length(pos_all));
    title(spec2)
    hold on
    set(gca,'FontSize',15)
end
    



        

    
    

