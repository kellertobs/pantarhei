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
    name = [outdir,RunID,'/',RunID,'_',num2str(step/abs(nop))];
    save([name,'.mat'],'res','res0','time','z','x','u','w','p','f',...
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
    
    if Nx == 1
        plot1dfields;
    else
        plot2dfields;
    end
    
    if svop % save
        name = [outdir,RunID,'/',RunID,'_sltn_',num2str(step/nop)];
        print(f1,'-dpdf','-r200','-opengl',name,'-loose');
        name = [outdir,RunID,'/',RunID,'_sgcp_',num2str(step/nop)];
        print(f2,'-dpdf','-r200','-opengl',name,'-loose');
        name = [outdir,RunID,'/',RunID,'_fltr_',num2str(step/nop)];
        print(f3,'-dpdf','-r200','-opengl',name,'-loose');
        name = [outdir,RunID,'/',RunID,'_coef_',num2str(step/nop)];
        print(f4,'-dpdf','-r200','-opengl',name,'-loose');
        
        if exist('f10','var')
            fig_pos               = f10.PaperPosition;
            f10.PaperPositionMode = 'manual';
            f10.PaperPosition     = [0 0 fig_pos(3) fig_pos(4)];
            f10.PaperSize         = [fig_pos(3) fig_pos(4)];
            name = [outdir,RunID,'/',RunID,'_itconv_',num2str(step/nop)];
            print(f10,'-dpdf','-r200','-opengl',name,'-loose');
        end

    end
end


close all;


