
clear all;

load ../out/plg60_dac20_mvp20_mms/TestCfflim_1.mat

h = D./Nvec;
Color = lines(4);

%%

figure;
set(gcf,'Position',[800,100,1200,1200], 'defaultlinelinewidth', 1.5);

PlotTitles = strcat({'Error in '}, {'pressure'; 'horizontal velocity'; 'vertical velocity'});
LegString = strcat('cfflim=10^', num2str(log10(cffvec')));

for vi = 2:4
    for ci = 1:length(cffvec)
        subplot(3,2,(vi-2)*2+1);
        loglog(h, NormErr(vi,:,ci), 'o:', 'Color', Color(ci,:));
        hold on;
        
        subplot(3,2,(vi-2)*2+2);
        hPlot(ci) = loglog(h, MaxErr(vi,:,ci), 'o:', 'Color', Color(ci,:));
        hold on;
        
    end
    subplot(3,2,(vi-2)*2+1);
    xlabel('Grid spacing, h'); ylabel('2-norm % error');
    title(PlotTitles{vi-1});
    
    subplot(3,2,(vi-2)*2+2);
    xlabel('Grid spacing, h'); ylabel('Maximum % error');
    title(PlotTitles{vi-1});
    if vi==2, legend(hPlot, LegString, 'Location', 'west'); end

end

%%
figure;
UniqueBeta = unique(betaOut(:));
b2 = betaOut(:);
h2 = repmat(h',4,1); 
cff2 = repelem(cffvec',4,1);
colors = parula(length(UniqueBeta)+1);

clear hPlot
for bi = 1:length(UniqueBeta)
    hPlot(bi) = loglog(h2(b2==UniqueBeta(bi)), cff2(b2==UniqueBeta(bi)), 'o', ...
        'MarkerSize', 10, 'MarkerFaceColor',colors(bi,:), 'color', colors(bi,:));
    hold on;
end

LegString = strcat('\beta=', num2str(UniqueBeta));
legend(hPlot, LegString, 'Location', 'east');
    xlabel('Grid spacing, h'); ylabel('cfflim');

