function all_data = read_data(data)
% input: real data collected from experiment
% output: different time dependent transcription rate (each as an array) at
% different mutation conditions

% every column is one type, from T2 to T10
WT = data(:,1);
mut1 = data(:,2);
mut2 = data(:,3);
mut3 = data(:,4);
mut12 = data(:,5);
mut13 = data(:,6);
mut23 = data(:,7);
mut123 = data(:,8);

all_data = [WT,mut1,mut2,mut3,mut12,mut13,mut23,mut123];




