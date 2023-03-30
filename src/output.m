% print output header
fprintf('\n***** output frame %d\n',step/abs(nop));

% get segregation velocities and compaction pressures
ustar = sum(omvx.*u,1);
wstar = sum(omvz.*w,1);
pstar = sum(omfc.*p,1);
usegr = (f(:,:,imx)+f(:,:,ipx))./2.*(u-ustar);
wsegr = (f(:,imz,:)+f(:,ipz,:))./2.*(w-wstar);
pcmpt =  f.*(p-pstar);

% get flux, transfer magnitudes
qvm   = sqrt((((qvxx(:,imz,imx)+qvxx(:,imz,ipx)+qvxx(:,ipz,imx)+qvxx(:,ipz,ipx))./4).^2 ...
    + ((qvzz(:,imz,imx)+qvzz(:,imz,ipx)+qvzz(:,ipz,imx)+qvzz(:,ipz,ipx))./4).^2 ...
    +   qvxz.^2 + qvxz.^2)./2);
qfm   = sqrt(((qfx(:,imz,:)+qfx(:,ipz,:))./2).^2 + ((qfz(:,:,imx)+qfz(:,:,ipx))./2).^2);
Gvm   = sqrt(((Gvx(:,imz,:)+Gvx(:,ipz,:))./2).^2 + ((Gvz(:,:,imx)+Gvz(:,:,ipx))./2).^2);
Gfm   = sqrt(Gf.^2);

% get segregation-compaction lengths
for i=1:NPHS
    for k=1:NPHS
        delta(i,k,:,:) = (f(i,:,:).^2./Cv(i,:,:).*f(k,:,:).^2./Cf(k,:,:)).^0.5;
    end
end

if svop % save
    name = [outdir,RunID,'/',RunID,'_frame_',num2str(step/abs(nop), '%03d')];
    save([name,'.mat'],'resflds','res','res0','time','z','x',...
        'u','w','p','f','ustar','wstar','pstar','usegr','wsegr','pcmpt',...
        'qvxx','qvzz','qvxz','qfx','qfz','Gvx','Gvz','Gf',...
        'Kv','Kf','Cv','Cf','delta');

    name = [outdir,RunID,'/',RunID,'_hist'];
    save([name,'.mat'], 'hist');
end


if (nop>0) %plot

    if Nx == 1
        plot1dfields;
    else
        plot2dfields;
    end

    f6 = figure(6); f6.Visible = figvis; clf;
    rc = lines(4);
    for jvar=1:4, semilogy(itvec(:,1),itvec(:,jvar+1),'.','MarkerSize',10,'Color',rc(jvar  ,:),'linewidth',1); hold on; end
    for jvar=5:7, semilogy(itvec(:,1),itvec(:,jvar+1),'*','MarkerSize', 6,'Color',rc(jvar-4,:),'linewidth',0.2); hold on; end
    semilogy(itvec(:,1), sum(itvec(:,2:end),2),'ko','MarkerSize',3,'linewidth',1); hold off;
    axis tight; drawnow;
    xlabel('iteration number'); ylabel('absolute residual');
    legend('u','w','p','f','ustar','wstar','pstar','total','NumColumns',2,'location','northoutside');

    if svop % save
        name = [outdir,RunID,'/',RunID,'_sltn_',num2str(step/nop)];
        print(f1,'-dpdf','-r200','-opengl',name,'-loose');
        name = [outdir,RunID,'/',RunID,'_sgcp_',num2str(step/nop)];
        print(f2,'-dpdf','-r200','-opengl',name,'-loose');
        name = [outdir,RunID,'/',RunID,'_fltr_',num2str(step/nop)];
        print(f3,'-dpdf','-r200','-opengl',name,'-loose');
        name = [outdir,RunID,'/',RunID,'_coef_',num2str(step/nop)];
        print(f4,'-dpdf','-r200','-opengl',name,'-loose');

        if exist('f5','var')
            name = [outdir,RunID,'/',RunID,'_vshr_',num2str(step/nop)];
            print(f5,'-dpdf','-r200','-opengl',name,'-loose');
        end

        fig_pos               = f6.PaperPosition;
        f6.PaperPositionMode = 'manual';
        f6.PaperPosition     = [0 0 fig_pos(3) fig_pos(4)];
        f6.PaperSize         = [fig_pos(3) fig_pos(4)];
        name = [outdir,RunID,'/',RunID,'_itconv_',num2str(step/nop)];
        print(f6,'-dpdf','-r200','-opengl',name,'-loose');

        if pltits
            fig_pos               = f10.PaperPosition;
            f10.PaperPositionMode = 'manual';
            f10.PaperPosition     = [0 0 fig_pos(3) fig_pos(4)];
            f10.PaperSize         = [fig_pos(3) fig_pos(4)];
            name = [outdir,RunID,'/',RunID,'_updp_',num2str(step/nop)];
            print(f10,'-dpdf','-r200','-opengl',name,'-loose');
        end

    end

end

close all;