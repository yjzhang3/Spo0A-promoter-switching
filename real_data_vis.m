%% read
all_data = read_data(promoteractivityA);

%% plot
for ii = 1:length(all_data(1,:))
    plot_data(all_data(:,ii))
    hold on
end

%% plot
plot_data(all_data(:,2))
