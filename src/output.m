% print output header
fprintf('\n***** output frame %d\n',step/abs(IO.nop));

% get flux, transfer magnitudes
qvm   = sqrt((((qvxx(:,NUM.imz,NUM.imx)+qvxx(:,NUM.imz,NUM.ipx)+qvxx(:,NUM.ipz,NUM.imx)+qvxx(:,NUM.ipz,NUM.ipx))./4).^2 ...
    + ((qvzz(:,NUM.imz,NUM.imx)+qvzz(:,NUM.imz,NUM.ipx)+qvzz(:,NUM.ipz,NUM.imx)+qvzz(:,NUM.ipz,NUM.ipx))./4).^2 ...
    +   qvxz.^2 + qvxz.^2)./2);
qfm   = sqrt(((qfx(:,NUM.imz,:)+qfx(:,NUM.ipz,:))./2).^2 + ((qfz(:,:,NUM.imx)+qfz(:,:,NUM.ipx))./2).^2);
Gvm   = sqrt(((Gvx(:,NUM.imz,:)+Gvx(:,NUM.ipz,:))./2).^2 + ((Gvz(:,:,NUM.imx)+Gvz(:,:,NUM.ipx))./2).^2);
Gfm   = sqrt(Gf.^2);

% get segregation-compaction lengths
for i=1:NUM.NPHS
    for k=1:NUM.NPHS
        delta(i,k,:,:) = (f(i,:,:).^2./Cv(i,:,:).*f(k,:,:).^2./Cf(k,:,:)).^0.5;
    end
end

if IO.svop % save
    name = [IO.outdir,IO.RunID,'/',IO.RunID,'_frame_',num2str(step/abs(IO.nop), '%03d')];
    save([name,'.mat'],'resflds','res','res0','time', ...
        'u','w','p','f','ustar','wstar','pstar','Du','Dw','Dp',...
        'qvxx','qvzz','qvxz','qfx','qfz','Gvx','Gvz','Gf',...
        'Kv','Kf','Cv','Cf','delta','PHS','NUM');

    name = [IO.outdir,IO.RunID,'/',IO.RunID,'_hist'];
    save([name,'.mat'], 'hist');
end


if (IO.nop>0) %plot

    if NUM.Nx == 1
        plot1dfields;
    else
        plot2dfields;
    end

    f6 = figure(6); f6.Visible = IO.figvis; clf;
    rc = lines(4);
    for jvar=1:3, semilogy(itvec(:,1),itvec(:,jvar+1),'.','MarkerSize',10,'Color',rc(jvar  ,:),'linewidth',1); hold on; end
    semilogy(itvec(:,1), sum(itvec(:,2:end),2),'ko','MarkerSize',3,'linewidth',1); hold off;
    axis tight; drawnow;
    xlabel('iteration number'); ylabel('absolute residual');
    legend('u','w','p','total','NumColumns',2,'location','northoutside');

    if IO.svop % save
        name = [IO.outdir,IO.RunID,'/',IO.RunID,'_sltn_',num2str(step/IO.nop)];
        print(f1,'-dpdf','-r200',name,'-loose');
        name = [IO.outdir,IO.RunID,'/',IO.RunID,'_sgcp_',num2str(step/IO.nop)];
        print(f2,'-dpdf','-r200',name,'-loose');
        name = [IO.outdir,IO.RunID,'/',IO.RunID,'_fltr_',num2str(step/IO.nop)];
        print(f3,'-dpdf','-r200',name,'-loose');
        name = [IO.outdir,IO.RunID,'/',IO.RunID,'_coef_',num2str(step/IO.nop)];
        print(f4,'-dpdf','-r200',name,'-loose');

        if exist('f5','var')
            name = [IO.outdir,IO.RunID,'/',IO.RunID,'_vshr_',num2str(step/IO.nop)];
            print(f5,'-dpdf','-r200',name,'-loose');
        end

        fig_pos               = f6.PaperPosition;
        f6.PaperPositionMode = 'manual';
        f6.PaperPosition     = [0 0 fig_pos(3) fig_pos(4)];
        f6.PaperSize         = [fig_pos(3) fig_pos(4)];
        name = [IO.outdir,IO.RunID,'/',IO.RunID,'_itconv_',num2str(step/IO.nop)];
        print(f6,'-dpdf','-r200',name,'-loose');

        if IO.pltits
            fig_pos               = f10.PaperPosition;
            f10.PaperPositionMode = 'manual';
            f10.PaperPosition     = [0 0 fig_pos(3) fig_pos(4)];
            f10.PaperSize         = [fig_pos(3) fig_pos(4)];
            name = [IO.outdir,IO.RunID,'/',IO.RunID,'_updp_',num2str(step/IO.nop)];
            print(f10,'-dpdf','-r200',name,'-loose');
        end

    end

end

% close all;