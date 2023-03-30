
% checks mass conservation by looking at div fv 

fx = 1/2*(f(:,:,imx) + f(:,:,ipx));   
fz = 1/2*(f(:,imz,:) + f(:,ipz,:));

Div_fv = diff(fx.*u,1,3)./h + diff(fz.*w,1,2)./h;

fprintf('\n norm Div_fv error:\n')
(norm( sum(Div_fv,1), 'fro' ) + TINY-TINY) / ( norm(Div_fv, 'fro') + TINY )

% plot div fv
dqffig = figure;
if ndim==1
    set(gcf,'Position',[10,10,300,600])
    plot(Div_fv, z, sum(Div_fv,1), z, 'k-'); axis tight; grid on;
elseif ndim==2
    set(gcf,'Position',[10,10,800,(NPHS+1)*300])
    for iphs = 1:NPHS
        subplot(NPHS+1,2,2*(iphs-1)+1); 
        imagesc(x, z, squeeze(Div_qf(iphs,:,:))); 
        axis xy equal tight; colormap(ocean); colorbar; 
        title(['$\phi^' num2str(iphs) '$']);
    end
    subplot(NPHS+1,2,2*NPHS+1); 
    imagesc(x, z, squeeze(sum(Div_qf,1))); 
    axis xy equal tight; colormap(ocean); colorbar; title('sum');

    subplot(NPHS+1,2,2*(1:NPHS+1));
    xmid = round(Nx/2);
    plot(Div_fv(:,:,xmid), z, sum(Div_fv(:,:,xmid),1), z, 'k-'); 
    axis tight; grid on;
    title('x midpoint');
end
sgtitle('$\nabla \cdot ( \phi \times v)$')
