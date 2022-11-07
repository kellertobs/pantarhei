function [] = SaveFigure(figname, fighandle, formattype, varargin)

if nargin < 2, fighandle = gcf; end
if nargin < 3, formattype = '-dpdf'; end

fig_pos = fighandle.PaperPosition;

fighandle.PaperPositionMode = 'manual';
fighandle.PaperPosition     = [0 0 fig_pos(3) fig_pos(4)];
fighandle.PaperSize         = [fig_pos(3) fig_pos(4)];

print(fighandle,  figname, formattype, varargin{:})
end
