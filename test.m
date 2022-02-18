energyi_normal = [1,1,10,0,0,0];
energyi_actA = [1,1,10,0,-0.77,0]; % TF of site A and mRNA attract
energyi_repB = [1,1,10,0,0,1]; % TF of site B and mRNA repel
mRNA_conc = 700;
TF_conc = 700;
nbd = 3;

TR_normal = transcription_rate(nbd,energyi_normal,TF_conc,mRNA_conc);
TR_actA = transcription_rate(nbd,energyi_actA,TF_conc,mRNA_conc);
TR_repB = transcription_rate(nbd,energyi_repB,TF_conc,mRNA_conc);
