
% calculate and plot true mms solution

fTrue = MMSsource('Calc_f', 0, X     , Z     , Tmf(:,1), Xmf(:,1), Zmf(:,1), Amf(:,1), dmf(:,1));
pTrue = MMSsource('Calc_p', 0, X     , Z     , Tmf(:,2), Xmf(:,2), Zmf(:,2), Amf(:,2), dmf(:,2));
uTrue = MMSsource('Calc_u', 0, XuGrid, ZuGrid, Tmf(:,3), Xmf(:,3), Zmf(:,3), Amf(:,3), dmf(:,3));
wTrue = MMSsource('Calc_w', 0, XwGrid, ZwGrid, Tmf(:,4), Xmf(:,4), Zmf(:,4), Amf(:,4), dmf(:,4));


axh = 6.00; axw = 7.50;
ahs = 0.90; avs = 0.90;
axb = 0.75; axt = 0.90;
axl = 1.75; axr = 0.90;
fh = axb + NPHS*axh + (NPHS-1)*avs + axt;
fw = axl + 4*axw + 3*ahs + axr;
TX = {'Interpreter','Latex'}; FS = {'FontSize',18};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',14};
UN = {'Units','Centimeters'};
    
% plot
f8 = figure(8); f8.Visible = 'off';
clf; colormap(ocean);
set(f8,UN{:},'Position',[4 4 fw fh]);
set(f8,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f8,'Color','w','InvertHardcopy','off');
set(f8,'Resize','off');

% predefine axis positions
axpos = zeros(NPHS*4, 4);
for nrow = 1:NPHS
    for ncol = 1:4
        axpos((nrow-1)*4 + ncol,:) = [axl+(ncol-1)*axw+(ncol-1)*ahs, axb+(NPHS-nrow)*axh+(NPHS-nrow)*avs, axw, axh];
    end
end

for n=1:NPHS
    axes(UN{:},'position',axpos(4*(n-1)+1,:));
    imagesc(x,z,squeeze(wTrue(n,:,:))); 
    axis xy equal tight; set(gca,'XTickLabel',[]);
    cb = colorbar; set(cb,TL{:},TS{:}); 
    set(gca,TL{:},TS{:});
    title(['$w_{mms}^',num2str(n), '$ [m/s]'],TX{:},FS{:});
    
    axes(UN{:},'position',axpos(4*(n-1)+2,:));
    imagesc(x,z,squeeze(uTrue(n,:,:))); 
    axis xy equal tight; set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    cb = colorbar; set(cb,TL{:},TS{:}); 
    set(gca,TL{:},TS{:}); 
    title(['$u_{mms}^',num2str(n), '$ [m/s]'],TX{:},FS{:});
    
    axes(UN{:},'position',axpos(4*(n-1)+3,:));
    imagesc(x,z,squeeze(pTrue(n,:,:))); 
    axis xy equal tight; set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    cb = colorbar; set(cb,TL{:},TS{:}); 
    set(gca,TL{:},TS{:}); 
    title(['$p_{mms}^',num2str(n), '$ [Pa]'],TX{:},FS{:});
    
    axes(UN{:},'position',axpos(4*(n-1)+4,:));
    imagesc(x,z,squeeze(fTrue(n,:,:))); 
    axis xy equal tight; set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    cb = colorbar; set(cb,TL{:},TS{:}); 
    set(gca,TL{:},TS{:}); 
    title(['$f_{mms}^',num2str(n), '$ [vol]'],TX{:},FS{:});
end

name = [outdir,RunID,'/',RunID,'_mmst'];
print(f8,'-dpdf','-r200','-opengl',name,'-loose');