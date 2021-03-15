
clear all;
load ../mms_out/olv_plg_bas/olv20_plg20_bas60_mms_NumConvTest_1.mat

N2 = Nvec.^2;
h = D./Nvec;

Color = lines(3);

%% convergence plots using grid spacing

figure;
set(gcf,'Position',[800,400,1200,400], 'defaultlinelinewidth', 1.5, 'defaultlinemarkersize', 15);

subplot(121);
for vi = 2:4
    hCP(vi-1) = loglog(h, NormErr(vi,:), '+:','linewidth',2); hold on;
end
ylim([0.08,20]);
AddConvOrderLine(1); AddConvOrderLine(2); AddConvOrderLine(3); 
legend(hCP, {'pressure, $p$', 'x-velocity, $u$', 'y-velocity, $w$'},'location','northwest');
xlabel('Grid spacing, h'); ylabel('\% error');
title('2-norm error');

subplot(122);
for vi = 2:4
    loglog(h, MaxErr(vi,:), '+:','linewidth',2); hold on;
end
ylim([0.08,20]);
AddConvOrderLine(1); AddConvOrderLine(2); AddConvOrderLine(3); 
xlabel('Grid spacing, h'); ylabel('\% error');
title('Maximum error');

%% convergence plots using number of grid points

figure;
set(gcf,'Position',[800,400,1200,400], 'defaultlinelinewidth', 1.5, 'defaultlinemarkersize', 15);

subplot(121);
for vi = 2:4
    hCP(vi-1) = loglog(Nvec, NormErr(vi,:), '+:','linewidth',2); hold on;
end
ylim([0.08,20]);
AddConvOrderLine(-1); AddConvOrderLine(-2); AddConvOrderLine(-3); 
legend(hCP, {'pressure, $p$', 'x-velocity, $u$', 'y-velocity, $w$'},'location','southwest');
xlabel('Number of grid points, N'); ylabel('\% error');
title('2-norm error');

subplot(122);
for vi = 2:4
    loglog(Nvec, MaxErr(vi,:), '+:','linewidth',2); hold on;
end
ylim([0.08,20]);
AddConvOrderLine(-1); AddConvOrderLine(-2); AddConvOrderLine(-3); 
xlabel('Number of grid points, N'); ylabel('\% error');
title('Maximum error');


%% plot betas used

figure;
semilogx(h(~isnan(NormErr(1,:))), betaOut(~isnan(NormErr(1,:))), '+:');

%% function to add convergence order lines of order k to a plot

function [] = AddConvOrderLine (k)

if k<1
    AxFrac = 0.9;
else
    AxFrac = 0.01;
end

xlimits = xlim; ylimits = ylim;
Intercept = ylimits(1)*10.^(AxFrac*(log10(ylimits(2)./ylimits(1))));
y1 = 10.^(log10(Intercept) - k.*log10(xlimits(1)));
y = y1.*10.^(k.*log10(xlimits));
hold on; axis manual
h = loglog(xlimits, y, '-', 'color', 0.7*ones(1,3), 'linewidth', 1);
hold off;
uistack(h,'bottom');
end