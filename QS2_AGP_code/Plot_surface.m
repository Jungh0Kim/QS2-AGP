% Jungho Kim
% final experimental designs (2-dimensional plot)

figure();
plot(d_opt_sv(:,1),d_opt_sv(:,2),'LineStyle','--','LineWidth',1.5,'Color','k'); hold on
contour(xs1, xs2, -mu_GP_g1,'levellist',0,'LineStyle','-.','LineWidth', 1.5,'Color',[0.09, 0.17, 0.24]);
contour(xs1, xs2, -mu_GP_g2,'levellist',0,'LineStyle','-.','LineWidth', 1.5,'Color','k','handlevisibility','off');
contour(xs1, xs2, reshape(g_p1,size(xs1)),'levellist', 0,'LineStyle','--', 'Color',[0.850, 0.225, 0.098],'LineWidth',1.5)
contour(xs1, xs2, reshape(g_p2,size(xs1)),'levellist', 0,'LineStyle','--', 'Color','r','LineWidth',1.5)
set(gca, 'layer', 'top');
scatter(x_doe(1:N_iniDoE,1),x_doe(1:N_iniDoE,2),20,'k','o')
for k=1:length(comp_id_sv)
    if comp_id_sv(k,:)==1
        scatter(x_doe(N_iniDoE+k,1),x_doe(N_iniDoE+k,2),35,'b','x','LineWidth',1.5)
    elseif comp_id_sv(k,:)==2
        scatter(x_doe(N_iniDoE+k,1),x_doe(N_iniDoE+k,2),35,'x','MarkerEdgeColor',[0.09, 0.17, 0.24],'LineWidth',1.5)
    end
end
scatter(d_opt_sv(:,1),d_opt_sv(:,2),20,'filled','k')
scatter(d_opt_sv(1,1),d_opt_sv(1,2),40,'filled','MarkerFaceColor','b')  % starting point
scatter(d_opt_sv(end,1),d_opt_sv(end,2),40,'filled','MarkerFaceColor',[0.66 0.12 0.74]) % optimum point
xlabel('\theta_{1}(=\mu_{{\itx}1}), {\itx}_1','fontsize',12,'fontname','arial')
ylabel('\theta_{2}(=\mu_{{\itx}2}), {\itx}_2','fontsize',12,'fontname','arial')
axis([0 3.7 0 4])
text(0.97*x0(1),0.97*x0(2),'\theta_{(0)}','HorizontalAlignment','center','Color','b','FontSize',10.5);
set(gcf,'color','w')
grid on; hold off;
