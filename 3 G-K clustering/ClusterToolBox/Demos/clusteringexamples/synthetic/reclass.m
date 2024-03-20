X = xlsread('I:\²©Ò»ÏÂ\01 ½·½­\Clustering\f_value_soilonly.xlsx');
[max_X,position]=max(X,[],2);
% Plot 
for i=1:20
s_plot=position(365*(i-1)+1:365*(i-1)+365,:);
figure
clrs = summer ;
clrs = clrs(end:-1:1,:);
colormap(clrs)
imagesc(s_plot')
xlabel('time')
ylabel('class')
%set(gca,'YTick',1:M,'YTickLabel',xlabel)
%colorbar
title(['Clustering pattern'])
hold on
%plot(M+0.5-flow/max(flow)*M,'k','LineWidth',2)
end