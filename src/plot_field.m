function [] = plot_field (x, v, vname)

load('ocean.mat','ocean');
NPHS = size(v,1);

figure;
set(gcf,'Position',[500,500,400*NPHS,350]);
hAx = tight_subplot(1,NPHS,0.03,[0.05,0.08],0.05);

for iphs = 1:NPHS
    
    axes(hAx(iphs));
    imagesc(x,x,squeeze(v(iphs,:,:)));
    colormap(ocean);
    axis xy equal tight;
    cb = colorbar; set(cb,'TickLabelInterpreter','Latex');
    set(gca,'XTickLabel',[]);
    if iphs>1, set(gca,'YTickLabel',[]); end
    
    if nargin>2, title(['$' vname '^' num2str(iphs) '$']); end
end



end