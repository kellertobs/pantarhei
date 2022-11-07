
function [hAx,pos] = default2dpanels (Nrow, Ncol, varargin)
% 
% [hAx,pos] = paperfigdefaults (Nrow, Ncol, varargin)
% 
% function to specify the axis size and margins for figures in paper 

ax.aspectratio = 1.00;

ax.height = 8.00; 
ax.bot    = 1.00; 
ax.top    = 1.50;
ax.gaph   = 1.50; 

ax.left   = 2.20;
ax.right  = 0.90;
ax.gapw   = 0.90;

% allow structure alteration
args = reshape(varargin, 2, []);
for ia = 1:size(args,2)
    ax.(args{1,ia}) = args{2,ia};
end

ax.width  = ax.aspectratio*ax.height + 2.00;

fh = ax.bot  + Nrow*ax.height + (Nrow-1)*ax.gaph + ax.top;
fw = ax.left + Ncol*ax.width  + (Ncol-1)*ax.gapw + ax.right; 

set(gcf,'Units','centimeters','Position',[5 20 fw fh]);
set(gcf,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(gcf,'defaultaxesfontsize',10);
set(gcf,'defaultlinelinewidth',1);

if Nrow>1 || Ncol>1
    [hAx,pos] = tight_subplot(Nrow,Ncol,[ax.gaph/fh,ax.gapw/fw], [ax.bot,ax.top]/fh, [ax.left,ax.right]/fw);
elseif Nrow==1 && Ncol==1
    hAx = axes('Units','centimeters','Position',[ax.left,ax.bot,ax.width,ax.height]);
end


end



