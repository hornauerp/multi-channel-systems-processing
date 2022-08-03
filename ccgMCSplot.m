function ccgMCSplot(ccgR,Bounds,ccg_params)
t = -ccg_params.duration/2:ccg_params.binSize:ccg_params.duration/2;
figure('Color','w');
tiledlayout(size(ccgR,2),size(ccgR,3))
for i = 1:size(ccgR,2)
    for j = 1:size(ccgR,3)
        nexttile
        bar(t,ccgR(:,i,j),'k','BarWidth',1)
        hold on
        plot(t,Bounds(:,i,j,1))
        plot(t,Bounds(:,i,j,2))
        set(gca,'FontSize',7)
        
        
        xlabel('Interval (ms)')
        xline(0,'w');
        % xlabel('Interval [ms]')
%         a1.YAxis.Visible = 'off';
        % xticks([-25])
        box off
        axis tight
        xlims = xlim;
        xticks([xlims(1) 0 xlims(2)])
        x_max = ccg_params.duration/2*1000;
        xticklabels([-x_max 0 x_max])
    end
end