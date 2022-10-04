
axh = 20.00; axw = 7.00;
ahs =  0.90; avs = 0.90;
axb =  0.80; axt = 0.90;
axl =  1.75; axr = 0.90;
fh = axb + 1*axh + 0*avs + axt;
fw = axl + 4*axw + 3*ahs + axr;

lcolors = [lines(NPHS);0,0,0];

%%
% initialize figure and axes
f1 = figure(1); f1.Visible = figvis;
clf; colororder(f1,lcolors);
set(f1,UN{:},'Position',[2 2 fw fh]);
set(f1,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f1,'Color','w','InvertHardcopy','off');
set(f1,'Resize','off');

% predefine axis positions
axpos = [axl*ones(4,1) + (0:3)'*axw+(0:3)'*ahs axb*ones(4,1) axw*ones(4,1) axh*ones(4,1)];

axes(UN{:},'position',axpos(1,:));
plot(squeeze(w),[z-h/2,z(end)+h/2],squeeze(wstar),[z-h/2,z(end)+h/2],':'); ylim(0.5*[-D(1),D(1)]);
ylabel('depth [m]'); set(gca,TL{:},TS{:}); 
xlabel('$w$ [m/s]',TX{:},FS{:}); legend([strcat('$w^',num2str((1:NPHS)'),'$');'$w^*$'],'Location','southoutside');

axes(UN{:},'position',axpos(2,:));
plot(squeeze(u(:,:,1)),z,squeeze(ustar(1,:,1)),z,':'); ylim(0.5*[-D(1),D(1)]);
set(gca,TL{:},TS{:}); set(gca,'YTickLabel',[]);
xlabel('$u$ [m/s]',TX{:},FS{:}); legend([strcat('$u^',num2str((1:NPHS)'),'$');'$u^*$'],'Location','southoutside');

axes(UN{:},'position',axpos(3,:));
plot(squeeze(p),z,squeeze(pstar),z,':'); ylim(0.5*[-D(1),D(1)]);
set(gca,TL{:},TS{:}); set(gca,'YTickLabel',[]);
xlabel('$p$ [Pa]',TX{:},FS{:}); legend([strcat('$p^',num2str((1:NPHS)'),'$');'$p^*$'],'Location','southoutside');

axes(UN{:},'position',axpos(4,:));
plot(squeeze(f-f0),z,0,0); ylim(0.5*[-D(1),D(1)]);
set(gca,TL{:},TS{:}); set(gca,'YTickLabel',[]);
xlabel('$f-f_0$',TX{:},FS{:}); legend([strcat('$f^',num2str((1:NPHS)'),'$');'N/A  '],'Location','southoutside');

annotation('textbox','String',['time = ',num2str(time,'%.1e'),' [s]'],'Position',[0.5,0.9,0.1,0.1],'LineStyle','none','HorizontalAlignment','center',TX{:},FS{:});

%%
f2 = figure(2); f2.Visible = figvis;
clf; colororder(f2,lcolors);
set(f2,UN{:},'Position',[4 4 fw fh]);
set(f2,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f2,'Color','w','InvertHardcopy','off');
set(f2,'Resize','off');

axes(UN{:},'position',axpos(1,:));
plot(squeeze(wsegr),[z-h/2,z(end)+h/2]); ylim(0.5*[-D(1),D(1)]);
ylabel('depth [m]'); set(gca,TL{:},TS{:}); 
xlabel('$w_\Delta$ [m/s]',TX{:},FS{:}); legend(strcat('$w_\Delta^',num2str((1:NPHS)'),'$'),'Location','southoutside');

axes(UN{:},'position',axpos(2,:));
plot(squeeze(usegr(:,:,1)),z); ylim(0.5*[-D(1),D(1)]);
set(gca,TL{:},TS{:}); set(gca,'YTickLabel',[]);
xlabel('$u_\Delta$ [m/s]',TX{:},FS{:}); legend(strcat('$u_\Delta^',num2str((1:NPHS)'),'$'),'Location','southoutside');

axes(UN{:},'position',axpos(3,:));
plot(squeeze(pcmpt),z); ylim(0.5*[-D(1),D(1)]);
set(gca,TL{:},TS{:}); set(gca,'YTickLabel',[]);
xlabel('$p_\Delta$ [Pa]',TX{:},FS{:}); legend(strcat('$p_\Delta^',num2str((1:NPHS)'),'$'),'Location','southoutside');

axes(UN{:},'position',axpos(4,:));
plot(squeeze((f-fo)./dt),z); ylim(0.5*[-D(1),D(1)]);
set(gca,TL{:},TS{:}); set(gca,'YTickLabel',[]);
xlabel('$\partial f/\partial t$ [vol/s]',TX{:},FS{:}); legend(strcat('$\partial f^',num2str((1:NPHS)'),'/\partial t$'),'Location','southoutside');

annotation('textbox','String',['time = ',num2str(time,'%.1e'),' [s]'],'Position',[0.5,0.9,0.1,0.1],'LineStyle','none','HorizontalAlignment','center',TX{:},FS{:});

%%
f3 = figure(3); f3.Visible = figvis;
clf; colororder(f3,lcolors);
set(f3,UN{:},'Position',[6 6 fw fh]);
set(f3,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f3,'Color','w','InvertHardcopy','off');
set(f3,'Resize','off');

axes(UN{:},'position',axpos(1,:));
plot(squeeze(qvzz-f.*p),z); ylim(0.5*[-D(1),D(1)]);
ylabel('depth [m]'); set(gca,TL{:},TS{:}); 
xlabel('$q_{v,zz}$ [Pa]',TX{:},FS{:}); legend(strcat('$q_{v,zz}^',num2str((1:NPHS)'),'$'),'Location','southoutside');

axes(UN{:},'position',axpos(2,:));
plot(squeeze(qfz-(f(:,imz,:)+f(:,ipz,:))./2.*w),[z-h/2,z(end)+h/2]); ylim(0.5*[-D(1),D(1)]);
set(gca,TL{:},TS{:}); set(gca,'YTickLabel',[]);
xlabel('$q_{f,z}$ [m/s]',TX{:},FS{:}); legend(strcat('$q_{f,z}^',num2str((1:NPHS)'),'$'),'Location','southoutside');

axes(UN{:},'position',axpos(3,:));
plot(squeeze(Gvz),[z-h/2,z(end)+h/2]); ylim(0.5*[-D(1),D(1)]);
set(gca,TL{:},TS{:}); set(gca,'YTickLabel',[]);
xlabel('$\Gamma_{v,z}$ [Pa/m]',TX{:},FS{:}); legend(strcat('$\Gamma_{v,z}^',num2str((1:NPHS)'),'$'),'Location','southoutside');

axes(UN{:},'position',axpos(4,:));
plot(squeeze(Gf),z); ylim(0.5*[-D(1),D(1)]);
set(gca,TL{:},TS{:}); set(gca,'YTickLabel',[]);
xlabel('$\Gamma_f$ [1/s]',TX{:},FS{:}); legend(strcat('$\Gamma_f^',num2str((1:NPHS)'),'$'),'Location','southoutside');

annotation('textbox','String',['time = ',num2str(time,'%.1e'),' [s]'],'Position',[0.5,0.9,0.1,0.1],'LineStyle','none','HorizontalAlignment','center',TX{:},FS{:});

%%
f4 = figure(4); f4.Visible = figvis;
clf; colororder(f4,lcolors);
set(f4,UN{:},'Position',[8 8 fw fh]);
set(f4,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f4,'Color','w','InvertHardcopy','off');
set(f4,'Resize','off');

axes(UN{:},'position',axpos(1,:));
semilogx(squeeze(Kv),z); ylim(0.5*[-D(1),D(1)]);
ylabel('depth [m]'); set(gca,TL{:},TS{:}); 
xlabel('$K_v$ [Pa s]',TX{:},FS{:}); legend(strcat('$K_v^',num2str((1:NPHS)'),'$'),'Location','southoutside');

axes(UN{:},'position',axpos(2,:));
semilogx(squeeze(Kf),z); ylim(0.5*[-D(1),D(1)]);
set(gca,TL{:},TS{:}); set(gca,'YTickLabel',[]);
xlabel('$K_f$ [m$^2$/Pa s]',TX{:},FS{:}); legend(strcat('$K_f^',num2str((1:NPHS)'),'$'),'Location','southoutside');

axes(UN{:},'position',axpos(3,:));
semilogx(squeeze(Cv),z); ylim(0.5*[-D(1),D(1)]);
set(gca,TL{:},TS{:}); set(gca,'YTickLabel',[]);
xlabel('$C_v$ [Pa s/m$^2$]',TX{:},FS{:}); legend(strcat('$C_v^',num2str((1:NPHS)'),'$'),'Location','southoutside');

axes(UN{:},'position',axpos(4,:));
semilogx(squeeze(Cf),z); ylim(0.5*[-D(1),D(1)]);
set(gca,TL{:},TS{:}); set(gca,'YTickLabel',[]);
xlabel('$C_f$ [1/Pa s]',TX{:},FS{:}); legend(strcat('$C_f^',num2str((1:NPHS)'),'$'),'Location','southoutside');

annotation('textbox','String',['time = ',num2str(time,'%.1e'),' [s]'],'Position',[0.5,0.9,0.1,0.1],'LineStyle','none','HorizontalAlignment','center',TX{:},FS{:});

