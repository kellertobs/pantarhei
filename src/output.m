% print output header
fprintf('\n***** output frame %d\n',step/abs(nop));

% get segregation velocities and compaction pressures
ustar = sum(omvx.*u,1);
wstar = sum(omvz.*w,1);
pstar = sum(omfc.*p,1);
usegr = (f(:,:,im)+f(:,:,ip))./2.*(u-ustar);
wsegr = (f(:,im,:)+f(:,ip,:))./2.*(w-wstar);
pcmpt =  f.*(p-pstar);

% get flux, transfer magnitudes
qvm   = sqrt((((qvxx(:,im,im)+qvxx(:,im,ip)+qvxx(:,ip,im)+qvxx(:,ip,ip))./4).^2 ...
    + ((qvzz(:,im,im)+qvzz(:,im,ip)+qvzz(:,ip,im)+qvzz(:,ip,ip))./4).^2 ...
    +   qvxz.^2 + qvxz.^2)./2);
qfm   = sqrt(((qfx(:,im,:)+qfx(:,ip,:))./2).^2 + ((qfz(:,:,im)+qfz(:,:,ip))./2).^2);
Gvm   = sqrt(((Gvx(:,im,:)+Gvx(:,ip,:))./2).^2 + ((Gvz(:,:,im)+Gvz(:,:,ip))./2).^2);
Gfm   = sqrt(Gf.^2);

% get segregation-compaction lengths
for i=1:NPHS
    for k=1:NPHS
        delta(i,k,:,:) = (f(i,:,:).^2./Cv(i,:,:).*f(k,:,:).^2./Cf(k,:,:)).^0.5;
    end
end

if nop<0 && svop % save
    name = ['../out/',RunID,'/',RunID,'_',num2str(step/abs(nop))];
    save([name,'.mat'],'res','time','x','u','w','p','f',...
        'ustar','wstar','pstar','usegr','wsegr','pcmpt',...
        'qvxx','qvzz','qvxz','qfx','qfz','Gvx','Gvz','Gf',...
        'Kv','Kf','Cv','Cf','delta');
end

if (nop>0) %plot
    
    figvis = 'off'; % toggle for figure visibility on/off
    
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
            iplt = (nrow-1)*4 + ncol;
            axpos(iplt,:) = [axl + (ncol-1)*axw+(ncol-1)*ahs axb+(NPHS+1-nrow)*axh+(NPHS+1-nrow)*avs axw axh];
        end
    end
     
    axes(UN{:},'position',axpos(1,:));
    imagesc(x,x,squeeze(wstar(1,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w^*$ [m/s]'],TX{:},FS{:});
    set(gca,'XTickLabel',[]);
    axes(UN{:},'position',axpos(2,:));
    imagesc(x,x,squeeze(ustar(1,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$u^*$ [m/s]'],TX{:},FS{:});
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    axes(UN{:},'position',axpos(3,:));
    imagesc(x,x,squeeze(pstar(1,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$p^*$ [Pa]' ],TX{:},FS{:});
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    axes(UN{:},'position',axpos(4,:))
    imagesc(x,x,squeeze(sum(f.*rho,1))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{\rho}$ [kg/m$^3$]'],TX{:},FS{:});
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    for n=2:NPHS+1
        axes(UN{:},'position',axpos(4*(n-1)+1,:));
        imagesc(x,x,squeeze(w(n-1,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w^',num2str(n-1),'$ [m/s]'],TX{:},FS{:});
        set(gca,'XTickLabel',[]);
        axes(UN{:},'position',axpos(4*(n-1)+2,:));
        imagesc(x,x,squeeze(u(n-1,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$u^',num2str(n-1),'$ [m/s]'],TX{:},FS{:});
        set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
        axes(UN{:},'position',axpos(4*(n-1)+3,:));
        imagesc(x,x,squeeze(p(n-1,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$p^',num2str(n-1),'$ [Pa]' ],TX{:},FS{:});
        set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
        axes(UN{:},'position',axpos(4*(n-1)+4,:));
        imagesc(x,x,squeeze(f(n-1,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$f^',num2str(n-1),'$ [vol]'],TX{:},FS{:});
        set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    end
    %drawnow;
    
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
            iplt = (nrow-1)*4 + ncol;
            axpos(iplt,:) = [axl + (ncol-1)*axw+(ncol-1)*ahs axb+(NPHS+1-nrow)*axh+(NPHS+1-nrow)*avs axw axh];
        end
    end
    
    for n=1:NPHS
        axes(UN{:},'position',axpos(4*(n-1)+1,:));
        imagesc(x,x,squeeze(wsegr(n,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w_\Delta^',num2str(n),'$ [m/s]'],TX{:},FS{:});
        set(gca,'XTickLabel',[]);
        axes(UN{:},'position',axpos(4*(n-1)+2,:));
        imagesc(x,x,squeeze(usegr(n,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$u_\Delta^',num2str(n),'$ [m/s]'],TX{:},FS{:});
        set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
        axes(UN{:},'position',axpos(4*(n-1)+3,:));
        imagesc(x,x,squeeze(pcmpt(n,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$p_\Delta^',num2str(n),'$ [Pa]' ],TX{:},FS{:});
        set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
        axes(UN{:},'position',axpos(4*(n-1)+4,:));
        imagesc(x,x,squeeze((f(n,:,:)-fo(n,:,:))./dt)); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\partial f^',num2str(n),'/\partial t$ [vol/s]'],TX{:},FS{:});
        set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    end
    %drawnow;
    
    f3 = figure(3); f3.Visible = figvis;
    clf; colormap(ocean);
    set(f3,UN{:},'Position',[6 6 fw fh]);
    set(f3,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(f3,'Color','w','InvertHardcopy','off');
    set(f3,'Resize','off');
    
    for n=1:NPHS
        axes(UN{:},'position',axpos(4*(n-1)+1,:));
        imagesc(x,x,squeeze(qvzz(n,:,:)-f(n,:,:).*p(n,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$q_{v,zz}^',num2str(n),'$ [Pa]'],TX{:},FS{:});
        set(gca,'XTickLabel',[]);
        axes(UN{:},'position',axpos(4*(n-1)+2,:));
        imagesc(x,x,squeeze(qfz (n,:,:)-(f(n,im,:)+f(n,ip,:))./2.*w(n,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$q_{f,z}^',num2str(n),'$ [m/s]'],TX{:},FS{:});
        set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
        axes(UN{:},'position',axpos(4*(n-1)+3,:));
        imagesc(x,x,squeeze(Gvz (n,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Gamma_{v,z}^',num2str(n),'$ [Pa/m]' ],TX{:},FS{:});
        set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
        axes(UN{:},'position',axpos(4*(n-1)+4,:));
        imagesc(x,x,squeeze(Gf  (n,:,:))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Gamma_f^',num2str(n),'$ [1/s]'],TX{:},FS{:});
        set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    end
    drawnow;
    
    f4 = figure(4); f4.Visible = figvis;
    clf; colormap(ocean);
    set(f4,UN{:},'Position',[8 8 fw fh]);
    set(f4,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(f4,'Color','w','InvertHardcopy','off');
    set(f4,'Resize','off');
    
    for n=1:NPHS
        axes(UN{:},'position',axpos(4*(n-1)+1,:));
        imagesc(x,x,squeeze(log10(Kv(n,:,:)))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$K_v^',num2str(n),'$ [Pas]'],TX{:},FS{:});
        set(gca,'XTickLabel',[]);
        axes(UN{:},'position',axpos(4*(n-1)+2,:));
        imagesc(x,x,squeeze(log10(Kf(n,:,:)))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$K_f^',num2str(n),'$ [m$^2$/Pas]'],TX{:},FS{:});
        set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
        axes(UN{:},'position',axpos(4*(n-1)+3,:));
        imagesc(x,x,squeeze(log10(Cv(n,:,:)))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$C_v^',num2str(n),'$ [Pas/m$^2$]' ],TX{:},FS{:});
        set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
        axes(UN{:},'position',axpos(4*(n-1)+4,:));
        imagesc(x,x,squeeze(log10(Cf(n,:,:)))); axis xy equal tight; cb = colorbar; set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$C_f^',num2str(n),'$ [1/Pas]' ],TX{:},FS{:});
        set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    end
    %drawnow;
    
%     figure(5); clf; colormap(ocean);
%     for n=1:NPHS
%         for nn=1:NPHS
%             subplot(NPHS,NPHS,(n-1)*NPHS+nn); imagesc(x,x,squeeze(log10(delta(n,nn,:,:)))); axis xy equal tight; cb = colorbar; title(['log$_{10}$ $\delta_{s/c}^{',[num2str(n),num2str(nn)],'}$ [m]'],TX{:},FS{:}); set(gca,TL{:},TS{:}); set(cb,TL{:},TS{:});
%         end
%     end
    %drawnow;
    
%     figure(6); clf; colormap(ocean);
%     for n=1:NPHS
%         subplot(NPHS,4,(n-1)*4+1); imagesc(x,x,squeeze(-res_w(n,:,:).*dtau_w(n,:,:)./(norm(squeeze(w(n,:,:)),2)./N+1e-32))); axis xy equal tight; cb = colorbar; title(['upd. $w^',num2str(n),'$ [m/s]'],TX{:},FS{:}); set(gca,TL{:},TS{:}); set(cb,TL{:},TS{:});
%         subplot(NPHS,4,(n-1)*4+2); imagesc(x,x,squeeze(-res_u(n,:,:).*dtau_u(n,:,:)./(norm(squeeze(u(n,:,:)),2)./N+1e-32))); axis xy equal tight; cb = colorbar; title(['upd. $u^',num2str(n),'$ [m/s]'],TX{:},FS{:}); set(gca,TL{:},TS{:}); set(cb,TL{:},TS{:});
%         subplot(NPHS,4,(n-1)*4+3); imagesc(x,x,squeeze(-res_p(n,:,:).*dtau_p(n,:,:)./(norm(squeeze(p(n,:,:)),2)./N+1e-32))); axis xy equal tight; cb = colorbar; title(['upd. $p^',num2str(n),'$ [Pa]' ],TX{:},FS{:}); set(gca,TL{:},TS{:}); set(cb,TL{:},TS{:});
%         subplot(NPHS,4,(n-1)*4+4); imagesc(x,x,squeeze(-res_f(n,:,:).*dtau_f(n,:,:)./(norm(squeeze(f(n,:,:)),2)./N+1e-32))); axis xy equal tight; cb = colorbar; title(['upd. $f^',num2str(n),'$ [vol]'],TX{:},FS{:}); set(gca,TL{:},TS{:}); set(cb,TL{:},TS{:});
%     end
    %drawnow;
    
    if svop % save
        name = ['../out/',RunID,'/',RunID,'_sltn_',num2str(step/nop)];
        print(f1,'-dpdf','-r200','-opengl',name,'-loose');
        name = ['../out/',RunID,'/',RunID,'_sgcp_',num2str(step/nop)];
        print(f2,'-dpdf','-r200','-opengl',name,'-loose');
        name = ['../out/',RunID,'/',RunID,'_fltr_',num2str(step/nop)];
        print(f3,'-dpdf','-r200','-opengl',name,'-loose');
        name = ['../out/',RunID,'/',RunID,'_coef_',num2str(step/nop)];
        print(f4,'-dpdf','-r200','-opengl',name,'-loose');
        
        if exist('f10','var')
            figure(f10);
            xlabel('Iterations'); ylabel('residual');
            f10.PaperPositionMode = 'manual';
            f10.PaperPosition = [0 0 f10.PaperPosition(3) f10.PaperPosition(4)];
            f10.PaperSize = [f10.PaperPosition(3) f10.PaperPosition(4)];
            name = ['../out/',RunID,'/',RunID,'_itconv_',num2str(step/nop)];
            print(f10,'-dpdf','-r200','-opengl',name,'-loose');
        end
    end
    close all;
end