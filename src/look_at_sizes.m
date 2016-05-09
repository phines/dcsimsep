% look 
data = csvread('../case30x4_results/k2/bo_sizes_1.csv',1);
size_MW = data(:,3);
size_brs = data(:,4);
figure(1); clf;
subplot(2,1,1); cla
hist(size_brs(size_brs>0),10);
xlabel('BO size in # of dep. outages');
ylabel('# of blackouts');
subplot(2,1,2); cla
hist(size_MW(size_MW>0),20);
xlabel('BO size in MW');
ylabel('# of blackouts');
