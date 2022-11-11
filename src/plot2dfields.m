
% axis properties
axh = 6.00; axw = 7.50;
ahs = 0.90; avs = 0.90;
axb = 0.75; axt = 0.90;
axl = 1.75; axr = 0.90;

UN = {'Units','Centimeters'};

%% plot u, w, p, f and reference fields

% initialize figure and axes
f1 = figure(1); f1.Visible = figvis;
clf; colormap(ocean);
fh = axb + (NPHS+1)*axh + NPHS*avs + axt;
fw = axl + 4*axw + 3*ahs + axr;
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
imagesc(x,z,squeeze(wstar(1,:,:))); 
format2dpanels(1, '$w^*$ [m/s]');

axes(UN{:},'position',axpos(2,:));
imagesc(x,z,squeeze(ustar(1,:,:))); 
format2dpanels(2, '$u^*$ [m/s]');

text(0,z(1),['time = ',num2str(time,'%.1e'),' [s]'],TX{:},FS{:},'Color','k','VerticalAlignment','bottom','HorizontalAlignment','center');

axes(UN{:},'position',axpos(3,:));
imagesc(x,z,squeeze(pstar(1,:,:))); 
format2dpanels(3, '$p^*$ [m/s]');

axes(UN{:},'position',axpos(4,:))
imagesc(x,z,squeeze(sum(f.*rho,1))); 
format2dpanels(4, '$\bar{\rho}$ [kg/m$^3$]');

for n=2:NPHS+1
    axes(UN{:},'position',axpos(4*(n-1)+1,:));
    imagesc(x,z,squeeze(w(n-1,:,:))); 
    format2dpanels(1, ['$w^',num2str(n-1),'$ [m/s]']);

    axes(UN{:},'position',axpos(4*(n-1)+2,:));
    imagesc(x,z,squeeze(u(n-1,:,:))); 
    format2dpanels(2, ['$u^',num2str(n-1),'$ [m/s]']);

    axes(UN{:},'position',axpos(4*(n-1)+3,:));
    imagesc(x,z,squeeze(p(n-1,:,:))); 
    format2dpanels(3, ['$p^',num2str(n-1),'$ [Pa]' ]);

    axes(UN{:},'position',axpos(4*(n-1)+4,:));
    imagesc(x,z,squeeze(f(n-1,:,:))); 
    format2dpanels(4, ['$f^',num2str(n-1),'$ [vol]']);
end

%% plot udelta, wdelta, pdelta, dfdt

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
    imagesc(x,z,squeeze(wsegr(n,:,:))); 
    format2dpanels(1, ['$w_\Delta^',num2str(n),'$ [m/s]']);

    axes(UN{:},'position',axpos(4*(n-1)+2,:));
    imagesc(x,z,squeeze(usegr(n,:,:))); 
    format2dpanels(2, ['$u_\Delta^',num2str(n),'$ [m/s]']);

    if n==1, text(0,z(1),['time = ',num2str(time,'%.1e'),' [s]'],TX{:},FS{:},'Color','k','VerticalAlignment','bottom','HorizontalAlignment','center'); end
    
    axes(UN{:},'position',axpos(4*(n-1)+3,:));
    imagesc(x,z,squeeze(pcmpt(n,:,:))); 
    format2dpanels(3, ['$p_\Delta^',num2str(n),'$ [Pa]' ]);
    
    axes(UN{:},'position',axpos(4*(n-1)+4,:));
    imagesc(x,z,squeeze((f(n,:,:)-fo(n,:,:))./dt)); 
    format2dpanels(4, ['$\partial f^',num2str(n),'/\partial t$ [vol/s]']);
end
drawnow;

%% plot qvzz, qfz, Gvz, Gf

f3 = figure(3); f3.Visible = figvis;
clf; colormap(ocean);
set(f3,UN{:},'Position',[6 6 fw fh]);
set(f3,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f3,'Color','w','InvertHardcopy','off');
set(f3,'Resize','off');

for n=1:NPHS
    axes(UN{:},'position',axpos(4*(n-1)+1,:));
    imagesc(x,z,squeeze(qvzz(n,:,:)-f(n,:,:).*p(n,:,:))); 
    format2dpanels(1, ['$q_{v,zz}^',num2str(n),'$ [Pa]']);

    axes(UN{:},'position',axpos(4*(n-1)+2,:));
    imagesc(x,z,squeeze(qfz (n,:,:)-(f(n,imz,:)+f(n,ipz,:))./2.*w(n,:,:))); 
    format2dpanels(2, ['$q_{f,z}^',num2str(n),'$ [m/s]']);
    
    if n==1, text(0,z(1),['time = ',num2str(time,'%.1e'),' [s]'],TX{:},FS{:},'Color','k','VerticalAlignment','bottom','HorizontalAlignment','center'); end
    
    axes(UN{:},'position',axpos(4*(n-1)+3,:));
    imagesc(x,z,squeeze(Gvz (n,:,:))); 
    format2dpanels(3, ['$\Gamma_{v,z}^',num2str(n),'$ [Pa/m]' ]);

    axes(UN{:},'position',axpos(4*(n-1)+4,:));
    imagesc(x,z,squeeze(Gf  (n,:,:))); 
    format2dpanels(4, ['$\Gamma_f^',num2str(n),'$ [1/s]']);
end
drawnow;

%% plot flux and transfer coeffs

f4 = figure(4); f4.Visible = figvis;
clf; colormap(ocean);
set(f4,UN{:},'Position',[8 8 fw fh]);
set(f4,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f4,'Color','w','InvertHardcopy','off');
set(f4,'Resize','off');

for n=1:NPHS
    axes(UN{:},'position',axpos(4*(n-1)+1,:));
    imagesc(x,z,squeeze(log10(Kv(n,:,:)))); 
    format2dpanels(1, ['$K_v^',num2str(n),'$ [Pa s]']);

    axes(UN{:},'position',axpos(4*(n-1)+2,:));
    imagesc(x,z,squeeze(log10(Kf(n,:,:)))); 
    format2dpanels(2, ['$K_f^',num2str(n),'$ [m$^2$/Pa s]']);

    if n==1, text(0,z(1),['time = ',num2str(time,'%.1e'),' [s]'],TX{:},FS{:},'Color','k','VerticalAlignment','bottom','HorizontalAlignment','center'); end
    
    axes(UN{:},'position',axpos(4*(n-1)+3,:));
    imagesc(x,z,squeeze(log10(Cv(n,:,:)))); 
    format2dpanels(3, ['$C_v^',num2str(n),'$ [Pa s/m$^2$]' ]);

    axes(UN{:},'position',axpos(4*(n-1)+4,:));
    imagesc(x,z,squeeze(log10(Cf(n,:,:)))); 
    format2dpanels(4, ['$C_f^',num2str(n),'$ [1/Pa s]' ]);
end



%% plot ushr, wshr velocity fields

f5 = figure(5); f5.Visible = figvis;
clf; colormap(ocean);
fh = axb + 2*axh + 1*avs + axt;
fw = axl + 2*axw + 1*ahs + axr;
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
imagesc(x,z,squeeze(wshr)); 
format2dpanels(1, '$w^{shr}$ [m/s]');

axes(UN{:},'position',axpos(2,:));
imagesc(x,z,squeeze(ushr)); 
format2dpanels(2, '$u^{shr}$ [m/s]');

text(0,z(1),['time = ',num2str(time,'%.1e'),' [s]'],TX{:},FS{:},'Color','k','VerticalAlignment','bottom','HorizontalAlignment','center');

axes(UN{:},'position',axpos(3,:));
imagesc(x,z,squeeze(wstar+wshr)); 
format2dpanels(3, '$w^* + w^{shr}$ [m/s]');

axes(UN{:},'position',axpos(4,:))
imagesc(x,z,squeeze(ustar+ushr)); 
format2dpanels(4, '$u^* + u^{shr}$ [m/s]');

%%

function [] = format2dpanels (plotcol, titletext)

TX = {         'Interpreter','Latex'}; FS = {'FontSize',18};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',14};

axis xy equal tight; 
set(gca,'XTickLabel',[],TL{:},TS{:});
cb = colorbar; set(cb,TL{:},TS{:}); 
title(titletext,TX{:},FS{:});

if plotcol==1
    ylabel('$z$ [m]', TX{:},TS{:}); 
else
    set(gca,'YTickLabel',[]); 
end

end
