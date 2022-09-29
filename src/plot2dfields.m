
axh = 6.00; axw = 7.50;
ahs = 0.90; avs = 0.90;
axb = 0.75; axt = 0.90;
axl = 1.75; axr = 0.90;
fh = axb + (NPHS+1)*axh + NPHS*avs + axt;
fw = axl + 4*axw + 3*ahs + axr;

%%
% initialize figure and axes
f1 = figure(1); f1.Visible = figvis;
clf; colormap(ocean);
set(f1,UN{:},'Position',[2 2 fw fh]);
set(f1,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f1,'Color','w','InvertHardcopy','off');
set(f1,'Resize','off');

% predefine axis positions
axpos = zeros((NPHS+1)*4, 4);
for nrow = 1:NPHS+1
    for ncol = 1:4
        axpos((nrow-1)*4 + ncol,:) = [axl + (ncol-1)*axw+(ncol-1)*ahs axb+(NPHS+1-nrow)*axh+(NPHS+1-nrow)*avs axw axh];
    end
end

axes(UN{:},'position',axpos(1,:));
imagesc(x,z,squeeze(wstar(1,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w^*$ [m/s]'],TX{:},FS{:});
set(gca,'XTickLabel',[]);
axes(UN{:},'position',axpos(2,:));
imagesc(x,z,squeeze(ustar(1,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$u^*$ [m/s]'],TX{:},FS{:});
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
text(0,z(1),['time = ',num2str(time,'%.1e'),' [s]'],TX{:},FS{:},'Color','k','VerticalAlignment','bottom','HorizontalAlignment','center');
axes(UN{:},'position',axpos(3,:));
imagesc(x,z,squeeze(pstar(1,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$p^*$ [Pa]' ],TX{:},FS{:});
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
axes(UN{:},'position',axpos(4,:))
imagesc(x,z,squeeze(sum(f.*rho,1))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{\rho}$ [kg/m$^3$]'],TX{:},FS{:});
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
for n=2:NPHS+1
    axes(UN{:},'position',axpos(4*(n-1)+1,:));
    imagesc(x,z,squeeze(w(n-1,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w^',num2str(n-1),'$ [m/s]'],TX{:},FS{:});
    set(gca,'XTickLabel',[]);
    axes(UN{:},'position',axpos(4*(n-1)+2,:));
    imagesc(x,z,squeeze(u(n-1,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$u^',num2str(n-1),'$ [m/s]'],TX{:},FS{:});
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    axes(UN{:},'position',axpos(4*(n-1)+3,:));
    imagesc(x,z,squeeze(p(n-1,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$p^',num2str(n-1),'$ [Pa]' ],TX{:},FS{:});
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    axes(UN{:},'position',axpos(4*(n-1)+4,:));
    imagesc(x,z,squeeze(f(n-1,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$f^',num2str(n-1),'$ [vol]'],TX{:},FS{:});
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
end

%%
f2 = figure(2); f2.Visible = figvis;
clf; colormap(ocean);
fh = axb + NPHS*axh + (NPHS-1)*avs + axt;
fw = axl + 4*axw + 3*ahs + axr;
set(f2,UN{:},'Position',[4 4 fw fh]);
set(f2,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f2,'Color','w','InvertHardcopy','off');
set(f2,'Resize','off');

% predefine axis positions
axpos = zeros(NPHS*4, 4);
for nrow = 1:NPHS
    for ncol = 1:4
        axpos((nrow-1)*4 + ncol,:) = [axl + (ncol-1)*axw+(ncol-1)*ahs axb+(NPHS-nrow)*axh+(NPHS-nrow)*avs axw axh];
    end
end

for n=1:NPHS
    axes(UN{:},'position',axpos(4*(n-1)+1,:));
    imagesc(x,z,squeeze(wsegr(n,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w_\Delta^',num2str(n),'$ [m/s]'],TX{:},FS{:});
    set(gca,'XTickLabel',[]);
    axes(UN{:},'position',axpos(4*(n-1)+2,:));
    imagesc(x,z,squeeze(usegr(n,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$u_\Delta^',num2str(n),'$ [m/s]'],TX{:},FS{:});
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    if n==1, text(0,z(1),['time = ',num2str(time,'%.1e'),' [s]'],TX{:},FS{:},'Color','k','VerticalAlignment','bottom','HorizontalAlignment','center'); end
    axes(UN{:},'position',axpos(4*(n-1)+3,:));
    imagesc(x,z,squeeze(pcmpt(n,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$p_\Delta^',num2str(n),'$ [Pa]' ],TX{:},FS{:});
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    axes(UN{:},'position',axpos(4*(n-1)+4,:));
    imagesc(x,z,squeeze((f(n,:,:)-fo(n,:,:))./dt)); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\partial f^',num2str(n),'/\partial t$ [vol/s]'],TX{:},FS{:});
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
end
%drawnow;

%%
f3 = figure(3); f3.Visible = figvis;
clf; colormap(ocean);
set(f3,UN{:},'Position',[6 6 fw fh]);
set(f3,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f3,'Color','w','InvertHardcopy','off');
set(f3,'Resize','off');

for n=1:NPHS
    axes(UN{:},'position',axpos(4*(n-1)+1,:));
    imagesc(x,z,squeeze(qvzz(n,:,:)-f(n,:,:).*p(n,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$q_{v,zz}^',num2str(n),'$ [Pa]'],TX{:},FS{:});
    set(gca,'XTickLabel',[]);
    axes(UN{:},'position',axpos(4*(n-1)+2,:));
    imagesc(x,z,squeeze(qfz (n,:,:)-(f(n,imz,:)+f(n,ipz,:))./2.*w(n,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$q_{f,z}^',num2str(n),'$ [m/s]'],TX{:},FS{:});
    if n==1, text(0,z(1),['time = ',num2str(time,'%.1e'),' [s]'],TX{:},FS{:},'Color','k','VerticalAlignment','bottom','HorizontalAlignment','center'); end
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    axes(UN{:},'position',axpos(4*(n-1)+3,:));
    imagesc(x,z,squeeze(Gvz (n,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Gamma_{v,z}^',num2str(n),'$ [Pa/m]' ],TX{:},FS{:});
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    axes(UN{:},'position',axpos(4*(n-1)+4,:));
    imagesc(x,z,squeeze(Gf  (n,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Gamma_f^',num2str(n),'$ [1/s]'],TX{:},FS{:});
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
end
drawnow;

%%
f4 = figure(4); f4.Visible = figvis;
clf; colormap(ocean);
set(f4,UN{:},'Position',[8 8 fw fh]);
set(f4,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f4,'Color','w','InvertHardcopy','off');
set(f4,'Resize','off');

for n=1:NPHS
    axes(UN{:},'position',axpos(4*(n-1)+1,:));
    imagesc(x,z,squeeze(log10(Kv(n,:,:)))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$K_v^',num2str(n),'$ [Pas]'],TX{:},FS{:});
    set(gca,'XTickLabel',[]);
    axes(UN{:},'position',axpos(4*(n-1)+2,:));
    imagesc(x,z,squeeze(log10(Kf(n,:,:)))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$K_f^',num2str(n),'$ [m$^2$/Pas]'],TX{:},FS{:});
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    if n==1, text(0,z(1),['time = ',num2str(time,'%.1e'),' [s]'],TX{:},FS{:},'Color','k','VerticalAlignment','bottom','HorizontalAlignment','center'); end
    axes(UN{:},'position',axpos(4*(n-1)+3,:));
    imagesc(x,z,squeeze(log10(Cv(n,:,:)))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$C_v^',num2str(n),'$ [Pas/m$^2$]' ],TX{:},FS{:});
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    axes(UN{:},'position',axpos(4*(n-1)+4,:));
    imagesc(x,z,squeeze(log10(Cf(n,:,:)))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$C_f^',num2str(n),'$ [1/Pas]' ],TX{:},FS{:});
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
end


fh = axb + 2*axh + 1*avs + axt;
fw = axl + 2*axw + 1*ahs + axr;

%%
f5 = figure(5); f5.Visible = figvis;
clf; colormap(ocean);
set(f5,UN{:},'Position',[10 10 fw fh]);
set(f5,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f5,'Color','w','InvertHardcopy','off');
set(f5,'Resize','off');

% predefine axis positions
axpos = zeros(4, 4);
for nrow = 1:2
    for ncol = 1:2
        axpos((nrow-1)*2 + ncol,:) = [axl + (ncol-1)*axw+(ncol-1)*ahs axb+(2-nrow)*axh+(2-nrow)*avs axw axh];
    end
end

axes(UN{:},'position',axpos(1,:));
imagesc(x,z,squeeze(wshr)); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w^{shr}$ [m/s]'],TX{:},FS{:});
set(gca,'XTickLabel',[]);
axes(UN{:},'position',axpos(2,:));
imagesc(x,z,squeeze(ushr)); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$u^{shr}$ [m/s]'],TX{:},FS{:});
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
text(0,z(1),['time = ',num2str(time,'%.1e'),' [s]'],TX{:},FS{:},'Color','k','VerticalAlignment','bottom','HorizontalAlignment','center');
axes(UN{:},'position',axpos(3,:));
imagesc(x,z,squeeze((wstar(1,1:end-1,:)+wstar(1,2:end,:))/2+wshr)); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w^* + w^{shr}$ [m/s]'],TX{:},FS{:});
set(gca,'XTickLabel',[]); 
axes(UN{:},'position',axpos(4,:))
imagesc(x,z,squeeze((ustar(1,:,1:end-1)+ustar(1,:,2:end))/2+ushr)); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$u^* + u^{shr}$ [m/s]'],TX{:},FS{:});
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);

