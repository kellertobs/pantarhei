
%% print individual residuals

% generate true manufactured solutions
fTrue = MMSsource('Calc_f', time, X     , Z     , Tmf(:,1), Xmf(:,1), Zmf(:,1), Amf(:,1), dmf(:,1));
pTrue = MMSsource('Calc_p', time, X     , Z     , Tmf(:,2), Xmf(:,2), Zmf(:,2), Amf(:,2), dmf(:,2));
uTrue = MMSsource('Calc_u', time, XuGrid, ZuGrid, Tmf(:,3), Xmf(:,3), Zmf(:,3), Amf(:,3), dmf(:,3));
wTrue = MMSsource('Calc_w', time, XwGrid, ZwGrid, Tmf(:,4), Xmf(:,4), Zmf(:,4), Amf(:,4), dmf(:,4));

% L2 norm error
fNormErr = calcnormerr(f, fTrue, 2, NPHS);
pNormErr = calcnormerr(p, pTrue, 2, NPHS);
uNormErr = calcnormerr(u, uTrue, 2, NPHS);
wNormErr = calcnormerr(w, wTrue, 2, NPHS);

% maximum error
fMaxErr = calcnormerr(f, fTrue, inf, NPHS);
pMaxErr = calcnormerr(p, pTrue, inf, NPHS);
uMaxErr = calcnormerr(u, uTrue, inf, NPHS);
wMaxErr = calcnormerr(w, wTrue, inf, NPHS);


fprintf(1, '\n\n');
fprintf(1,'2-norm [f,p,u,w] percent error:\n');
for iphs = 1:NPHS
    fprintf(1, '%.6f, %.6f, %.6f, %.6f.\n', ...
        fNormErr(iphs), pNormErr(iphs), uNormErr(iphs), wNormErr(iphs));
end
fprintf(1,'\n Max [f,p,u,w] percent error:\n');
for iphs = 1:NPHS
    fprintf(1, '%.6f, %.6f, %.6f, %.6f.\n', ...
        fMaxErr(iphs), pMaxErr(iphs), uMaxErr(iphs), wMaxErr(iphs));
end
fprintf(1, '\n\n');



%% plot

f7 = figure(7); f7.Visible = figvis;
clf; colormap(ocean);
fh = axb + NPHS*axh + (NPHS-1)*avs + axt;
fw = axl + 4*axw + 3*ahs + axr;
set(f7,UN{:},'Position',[4 4 fw fh]);
set(f7,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f7,'Color','w','InvertHardcopy','off');
set(f7,'Resize','off');


% predefine axis positions
axpos = zeros(NPHS*4, 4);
for nrow = 1:NPHS
    for ncol = 1:4
        axpos((nrow-1)*4 + ncol,:) = [axl + (ncol-1)*axw+(ncol-1)*ahs axb+(NPHS-nrow)*axh+(NPHS-nrow)*avs axw axh];
    end
end

for n=1:NPHS
    axes(UN{:},'position',axpos(4*(n-1)+1,:));
    imagesc(x,z,squeeze(wTrue(n,:,:)-w(n,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w_{mms}^',num2str(n),'-w^',num2str(n) '$ [m/s]'],TX{:},FS{:});
    set(gca,'XTickLabel',[]);
    axes(UN{:},'position',axpos(4*(n-1)+2,:));
    imagesc(x,z,squeeze(uTrue(n,:,:)-u(n,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$u_{mms}^',num2str(n),'-u^',num2str(n) '$ [m/s]'],TX{:},FS{:});
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    if n==1, text(0,z(1),['time = ',num2str(time,'%.1e'),' [s]'],TX{:},FS{:},'Color','k','VerticalAlignment','bottom','HorizontalAlignment','center'); end
    axes(UN{:},'position',axpos(4*(n-1)+3,:));
    imagesc(x,z,squeeze(pTrue(n,:,:)-p(n,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$p_{mms}^',num2str(n),'-p^',num2str(n) '$ [Pa]'],TX{:},FS{:});
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    axes(UN{:},'position',axpos(4*(n-1)+4,:));
    imagesc(x,z,squeeze(fTrue(n,:,:)-f(n,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$f_{mms}^',num2str(n),'-f^',num2str(n) '$ [vol]'],TX{:},FS{:});
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
end

name = [outdir,RunID,'/',RunID,'_mmse_',num2str(step/nop)];
print(f7,'-dpdf','-r200','-opengl',name,'-loose');

%%

function [ne] = calcnormerr (v, vtrue, no, NPHS)
TINY = 1e-16;
ne   = nan(NPHS, 1);

for iphs = 1:NPHS
    ne(iphs) = 100*norm(squeeze(v(iphs,:,:)+TINY - vtrue(iphs,:,:)-TINY),no)./(norm(squeeze(vtrue(iphs,:,:)+TINY),no));
end
end
