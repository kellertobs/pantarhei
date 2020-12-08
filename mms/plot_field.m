function [] = plot_field (v,plottitle,X,Z)

phase = {'solid','liquid','gas'};
Nplt  = size(v,1);

figure;
set(gcf,'Position',[600,500,600*Nplt,500]);
for iplt = 1:Nplt
    subplot(1,Nplt,iplt); 
    
    if nargin>2
        surf(squeeze(X(1,:,:)),squeeze(Z(1,:,:)),squeeze(v(iplt,:,:))); 
    else
        surf(squeeze(v(iplt,:,:)));
    end
    
    colorbar;
    vdiff = [min(v(iplt,:,:),[],'all'), max(v(iplt,:,:),[],'all')];
    if vdiff(2)/vdiff(1)>100, set(gca,'ZScale','log'); end
    
    xlabel('x [m]'); ylabel('y [m]'); 
    title(phase{iplt});
    shading interp
end

if nargin>1, sgtitle(plottitle); end

end