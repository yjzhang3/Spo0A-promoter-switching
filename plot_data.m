function success = plot_data(data)

plot(data,'LineWidth',4)
xlabel('time (hours)')
ylabel('promoter activity')
set(gca,'FontSize',12)

success=1;

end