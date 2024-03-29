
axh = 20.00; axw = 7.00;
ahs =  0.90; avs = 0.90;
axb =  0.80; axt = 0.90;
axl =  1.75; axr = 0.90;
fh  = axb + 1*axh + 0*avs + axt;
fw  = axl + 4*axw + 3*ahs + axr;
lcolors = [lines(NUM.NPHS);0,0,0];

UN = {'Units','Centimeters'};
TX = {'Interpreter','Latex'}; FS = {'FontSize',18};

%% plot u, w, p, f and reference fields

% initialize figure and axes
f1 = figure(1); f1.Visible = IO.figvis;
clf; colororder(f1,lcolors);
set(f1,UN{:},'Position',[2 2 fw fh]);
set(f1,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f1,'Color','w','InvertHardcopy','off');
set(f1,'Resize','off');

% predefine axis positions
axpos = [axl*ones(4,1) + (0:3)'*axw+(0:3)'*ahs axb*ones(4,1) axw*ones(4,1) axh*ones(4,1)];

axes(UN{:},'position',axpos(1,:));
plot(squeeze(w),[NUM.z-NUM.h/2,NUM.z(end)+NUM.h/2],squeeze(wstar),[NUM.z-NUM.h/2,NUM.z(end)+NUM.h/2],':'); 
format1dpanels(1, NUM.L, '$w$ [m/s]', [strcat('$w^',num2str((1:NUM.NPHS)'),'$');'$w^*$'])

axes(UN{:},'position',axpos(2,:));
plot(squeeze(u(:,:,1)),NUM.z,squeeze(ustar(1,:,1)),NUM.z,':'); 
format1dpanels(2, NUM.L, '$u$ [m/s]', [strcat('$u^',num2str((1:NUM.NPHS)'),'$');'$u^*$']);

axes(UN{:},'position',axpos(3,:));
plot(squeeze(p),NUM.z,squeeze(pstar),NUM.z,':');
format1dpanels(3, NUM.L, '$p$ [Pa]', [strcat('$p^',num2str((1:NUM.NPHS)'),'$');'$p^*$']);

axes(UN{:},'position',axpos(4,:));
plot(squeeze(f-fin),NUM.z,0,0);
format1dpanels(4, NUM.L, '$f-f_0$', [strcat('$f^',num2str((1:NUM.NPHS)'),'$');'N/A  ']);

annotation('textbox','String',['time = ',num2str(time,'%.1e'),' [s]'],'Position',[0.5,0.9,0.1,0.1],'LineStyle','none','HorizontalAlignment','center',TX{:},FS{:});

%% plot udelta, wdelta, pdelta, dfdt

f2 = figure(2); f2.Visible = IO.figvis;
clf; colororder(f2,lcolors);
set(f2,UN{:},'Position',[4 4 fw fh]);
set(f2,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f2,'Color','w','InvertHardcopy','off');
set(f2,'Resize','off');

axes(UN{:},'position',axpos(1,:));
plot(squeeze(fz.*Dw),[NUM.z-NUM.h/2,NUM.z(end)+NUM.h/2]); 
format1dpanels(1, NUM.L, '$w_\Delta$ [m/s]', strcat('$w_\Delta^',num2str((1:NUM.NPHS)'),'$'));

axes(UN{:},'position',axpos(2,:));
plot(squeeze(fx(:,:,1).*u(:,:,1)),NUM.z);
format1dpanels(2, NUM.L, '$u_\Delta$ [m/s]', strcat('$u_\Delta^',num2str((1:NUM.NPHS)'),'$'));

axes(UN{:},'position',axpos(3,:));
plot(squeeze(f.*Dp),NUM.z); 
format1dpanels(3, NUM.L, '$p_\Delta$ [Pa]', strcat('$p_\Delta^',num2str((1:NUM.NPHS)'),'$'));

axes(UN{:},'position',axpos(4,:));
plot(squeeze((f-fo)./dt),NUM.z); 
format1dpanels(4, NUM.L, '$\partial f/\partial t$ [vol/s]', strcat('$\partial f^',num2str((1:NUM.NPHS)'),'/\partial t$'));

annotation('textbox','String',['time = ',num2str(time,'%.1e'),' [s]'],'Position',[0.5,0.9,0.1,0.1],'LineStyle','none','HorizontalAlignment','center',TX{:},FS{:});

%% plot qvzz, qfz, Gvz, Gf

f3 = figure(3); f3.Visible = IO.figvis;
clf; colororder(f3,lcolors);
set(f3,UN{:},'Position',[6 6 fw fh]);
set(f3,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f3,'Color','w','InvertHardcopy','off');
set(f3,'Resize','off');

axes(UN{:},'position',axpos(1,:));
plot(squeeze(qvzz-f.*p),NUM.z); 
format1dpanels(1, NUM.L, '$q_{v,zz}$ [Pa]', strcat('$q_{v,zz}^',num2str((1:NUM.NPHS)'),'$'));

axes(UN{:},'position',axpos(2,:));
plot(squeeze(Div_qf),NUM.z, squeeze(sum(Div_qf,1)),NUM.z, 'k-'); 
format1dpanels(2, NUM.L, '$\nabla \cdot q_{f}$ [m/s]', strcat('$\nabla \cdot q_{f}^',[num2str((1:NUM.NPHS)');'t'],'$'));

axes(UN{:},'position',axpos(3,:));
plot(squeeze(Gvz),[NUM.z-NUM.h/2,NUM.z(end)+NUM.h/2], squeeze(sum(Gvz)),[NUM.z-NUM.h/2,NUM.z(end)+NUM.h/2], 'k-'); 
format1dpanels(3, NUM.L, '$\Gamma_{v,z}$ [Pa/m]', strcat('$\Gamma_{v,z}^',[num2str((1:NUM.NPHS)');'t'],'$'));

axes(UN{:},'position',axpos(4,:));
plot(squeeze(Gf),NUM.z, squeeze(sum(Gf,1)),NUM.z, 'k-'); 
format1dpanels(4, NUM.L, '$\Gamma_f$ [1/s]', strcat('$\Gamma_f^',[num2str((1:NUM.NPHS)');'t'],'$'));

annotation('textbox','String',['time = ',num2str(time,'%.1e'),' [s]'],'Position',[0.5,0.9,0.1,0.1],'LineStyle','none','HorizontalAlignment','center',TX{:},FS{:});

%% plot flux and transfer coeffs

f4 = figure(4); f4.Visible = IO.figvis;
clf; colororder(f4,lcolors);
set(f4,UN{:},'Position',[8 8 fw fh]);
set(f4,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f4,'Color','w','InvertHardcopy','off');
set(f4,'Resize','off');

axes(UN{:},'position',axpos(1,:));
semilogx(squeeze(Kv),NUM.z); 
format1dpanels(1, NUM.L, '$K_v$ [Pa s]', strcat('$K_v^',num2str((1:NUM.NPHS)'),'$'));

axes(UN{:},'position',axpos(2,:));
semilogx(squeeze(Kf),NUM.z);
format1dpanels(2, NUM.L, '$K_f$ [m$^2$/Pa s]', strcat('$K_f^',num2str((1:NUM.NPHS)'),'$'));

axes(UN{:},'position',axpos(3,:));
semilogx(squeeze(Cv),NUM.z); 
format1dpanels(3, NUM.L, '$C_v$ [Pa s/m$^2$]', strcat('$C_v^',num2str((1:NUM.NPHS)'),'$'));

axes(UN{:},'position',axpos(4,:));
semilogx(squeeze(Cf),NUM.z); 
format1dpanels(4, NUM.L, '$C_f$ [1/Pa s]', strcat('$C_f^',num2str((1:NUM.NPHS)'),'$'));

annotation('textbox','String',['time = ',num2str(time,'%.1e'),' [s]'],'Position',[0.5,0.9,0.1,0.1],'LineStyle','none','HorizontalAlignment','center',TX{:},FS{:});

%%

function [] = format1dpanels (plotcol, L, vartext, legtext)

TX = {         'Interpreter','Latex'}; FS = {'FontSize',18};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',14};

% ylim(0.5*[-L(1),L(1)]);
axis tight;
set(gca,TL{:},TS{:});
xlabel(vartext,TX{:},FS{:}); 
legend(legtext,TX{:},FS{:},'Location','southoutside');

if plotcol==1
    ylabel('$z$ [m]', TX{:},TS{:}); 
else
    set(gca,'YTickLabel',[]); 
end
end
