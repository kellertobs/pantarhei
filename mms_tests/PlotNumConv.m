
clear all;
load ../out/plg60_dac20_mvp20_mms/NumConvTest_6.mat

N2 = Nvec.^2;
h = D./Nvec;

Color = lines(3);

%% plot convergence plots

figure;
set(gcf,'Position',[800,400,1200,400], 'defaultlinelinewidth', 1.5, 'defaultlinemarkersize', 10);

subplot(121);
for vi = 2:4
    hCP(vi-1) = loglog(h, NormErr(vi,:), '+:'); hold on;
end
AddConvOrderLine(1); AddConvOrderLine(2); AddConvOrderLine(3); AddConvOrderLine(4);
legend(hCP, {'Pressure', 'Horiz Vel', 'Vertical Vel'},'location','northwest');
xlabel('Grid spacing, h'); ylabel('2-norm \% error');

subplot(122);
for vi = 2:4
    loglog(h, MaxErr(vi,:), '+:'); hold on;
end
AddConvOrderLine(1); AddConvOrderLine(2); AddConvOrderLine(3); AddConvOrderLine(4);
xlabel('Grid spacing, h'); ylabel('Maximum \% error');

%% plot betas used

figure;
semilogx(h(~isnan(NormErr(1,:))), betaOut(~isnan(NormErr(1,:))), '+:');

%% function to add convergence order lines of order k to a plot

function [] = AddConvOrderLine (k)

if k<1
    AxFrac = 0.8;
else
    AxFrac = 0.2;
end

xlimits = xlim; ylimits = ylim;
Intercept = ylimits(1)*10.^(AxFrac*(log10(ylimits(2)./ylimits(1))));
y1 = 10.^(log10(Intercept) - k.*log10(xlimits(1)));
y = y1.*10.^(k.*log10(xlimits));
hold on; axis manual
h = loglog(xlimits, y, '-', 'color', 0.6*ones(1,3), 'linewidth', 1);
hold off;
uistack(h,'bottom');
end