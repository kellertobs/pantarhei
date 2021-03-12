
%% print individual residuals

% generate manufactured solutions
fTrue = MMSsource('Calc_f', time, X     , Z     , Tmf(:,1), Xmf(:,1), Zmf(:,1), Amf(:,1), dmf(:,1));
pTrue = MMSsource('Calc_p', time, X     , Z     , Tmf(:,2), Xmf(:,2), Zmf(:,2), Amf(:,2), dmf(:,2));
uTrue = MMSsource('Calc_u', time, XuGrid, ZuGrid, Tmf(:,3), Xmf(:,3), Zmf(:,3), Amf(:,3), dmf(:,3));
wTrue = MMSsource('Calc_w', time, XwGrid, ZwGrid, Tmf(:,4), Xmf(:,4), Zmf(:,4), Amf(:,4), dmf(:,4));

% 2-norm error
fNormErr = 100*norm(f(:)-fTrue(:),2)./norm(fTrue(:),2);
pNormErr = 100*norm(p(:)-pTrue(:),2)./norm(pTrue(:),2);
uNormErr = 100*norm(u(:)-uTrue(:),2)./norm(uTrue(:),2);
wNormErr = 100*norm(w(:)-wTrue(:),2)./norm(wTrue(:),2);

% maximum error
fMaxErr = 100*norm(f(:)-fTrue(:),inf)./norm(fTrue(:),inf);
pMaxErr = 100*norm(p(:)-pTrue(:),inf)./norm(pTrue(:),inf);
uMaxErr = 100*norm(u(:)-uTrue(:),inf)./norm(uTrue(:),inf);
wMaxErr = 100*norm(w(:)-wTrue(:),inf)./norm(wTrue(:),inf);

fprintf(1, '\n\n');
fprintf(1,'2-norm [f,p,u,w] percent error: \n%.6e, %.6f, %.6f, %.6f.\n\n', ...
    fNormErr, pNormErr, uNormErr, wNormErr); 
fprintf(1,'Max [f,p,u,w] percent error: \n%.6e, %.6f, %.6f, %.6f.', ...
    fMaxErr, pMaxErr, uMaxErr, wMaxErr);
fprintf(1, '\n\n');



%% plot

% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',18};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',14};
UN = {'Units','Centimeters'};

axh = 6.00; axw = 7.50;
ahs = 0.90; avs = 0.90;
axb = 0.75; axt = 0.90;
axl = 1.75; axr = 0.90;
fh = axb + (NPHS+1)*axh + NPHS*avs + axt;
fw = axl + 4*axw + 3*ahs + axr;


f7 = figure(7); clf; colormap(ocean);
fh = axb + NPHS*axh + (NPHS-1)*avs + axt;
fw = axl + 4*axw + 3*ahs + axr;
set(f7,UN{:},'Position',[4 4 fw fh]);
set(f7,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f7,'Color','w','InvertHardcopy','off');
set(f7,'Resize','off');
for n=1:NPHS
    ax((n-1)*4+1) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+(NPHS-n)*axh+(NPHS-n)*avs axw axh]);
    ax((n-1)*4+2) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+(NPHS-n)*axh+(NPHS-n)*avs axw axh]);
    ax((n-1)*4+3) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+(NPHS-n)*axh+(NPHS-n)*avs axw axh]);
    ax((n-1)*4+4) = axes(UN{:},'position',[axl+3*axw+3*ahs axb+(NPHS-n)*axh+(NPHS-n)*avs axw axh]);
end

for n=1:NPHS
    axes(ax((n-1)*4+1));
    imagesc(x,x,squeeze(wTrue(n,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w_\textrm{true}^',num2str(n),'$ [m/s]'],TX{:},FS{:});
    set(gca,'XTickLabel',[]);
    axes(ax((n-1)*4+2));
    imagesc(x,x,squeeze(uTrue(n,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$u_\textrm{true}^',num2str(n),'$ [m/s]'],TX{:},FS{:});
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    axes(ax((n-1)*4+3));
    imagesc(x,x,squeeze(pTrue(n,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$p_\textrm{true}^',num2str(n),'$ [Pa]' ],TX{:},FS{:});
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    axes(ax((n-1)*4+4));
    imagesc(x,x,squeeze(fTrue(n,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$f_\textrm{true}^',num2str(n),'$' ],TX{:},FS{:});
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
end
drawnow;