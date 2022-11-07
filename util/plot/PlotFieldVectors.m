function [fig, scl] = PlotFieldVectors (folder, RunID, varname, varmat, t, x, z, varargin)
%
% fig = PlotFieldVectors (folder, RunID, varname, varmat, varargin)
% This is very similar to PlotFieldwTime, except that we can have velocity
% vectors overlain
%
% examples
% fig = PlotFieldVectors(folder, RunID, {'f','u','w'});
% fig = PlotFieldVectors(folder, RunID, {'f','u','w'}, [], 'Nstd', 5);
%
%
% INPUTS
% folder    folder name where output folder is stored
% RunID     name of the run so that the total path to a mat file is
%               [folder  RunID / RunID_0.mat]
% varname   cell vector of names of 3 variables to plot:
%               varname{1}: field for the background
%               varname{2}: field for the horizontal vectors
%               varname{3}: field for the vertical vectors
%               NB if you want certain axis limit options, varname has to be precise
% varmat    cell vector of 3 variable matrices to plot
%               if mising or empty, load the variable of varname from file
% t, x      optional arguments (can be empty) if you already have varmat
%           populated so we don't have to load mat files
% varargin  plotting options (see defopts)
%
% 
% OUTPUTS
% fig       figure handle
% varmat    cell of variables plotted {3x1}
% x         x positions
% 
% 
% DEFAULT OPTIONS
% opt.fname  = '';        % specified filename
% opt.save   = 0;         % whether to save
% 
% opt.Nquiv  = 21;        % number of quiver arrows (downsampling)
% opt.bcind  = 0.05;      % exclude arrows from boundary
% opt.uquiv  = 1;         % whether to plot same length quiver for all plts
% 
% opt.zline  = [];        % whether to plot a horizontal zero line
% opt.xdsc   = 0;         % whether to divide by initial max dsc
% opt.rmadv  = 0;         % whether to remove box advection in periodic BC
% opt.boxind = [];        % whether to plot a subset of domain 
% 
% opt.Nplt   = 5;         % number of panels to plot
% opt.iPlt   = [];        % which panels to plot
% opt.uaxes  = 1;         % whether to have uniform axes for all panels
% opt.Nstd   = 3;         % number of stds from mean for axis limits
% 
% YQW, editted 4 May 2021
%

%  check inputs
if nargin<4, varmat = cell(3,1); t = []; x = []; z = []; end
if isempty(varmat), varmat = cell(3,1); t = []; x = []; z = [];  end

% load colormap
load('ocean.mat', 'ocean');

% get output mat files
fp = GetOutputMatFiles(folder, RunID);

% assign variables to background, u vectors, w vectors
[t, x, z, bg, u, w] = LoadPlotVars(folder, RunID, varname, varmat, t, x, z);
Nf = size(bg,4);

% get plotting options
opt = defopts(Nf, varargin{:});

if (opt.xdsc)
    load(fp, 'delta0');
    x = x./max(delta0(:));
    z = z./max(delta0(:));
end

% check sizes and make sure bg, u, w are the same size
NPHS = max(size(bg,1), size(u,1));
if NPHS>1
    if size(bg,1)==1, bg = bg.*ones(NPHS,1); end
    if size( u,1)==1,  u =  u.*ones(NPHS,1); end
    if size( w,1)==1,  w =  w.*ones(NPHS,1); end
end

% check BCs about whether to remove the mean velocity of the box
if opt.rmadv
    load(fp, 'BC');
    if strcmp(BC, 'periodic'), w = w - mean(w,[2,3]); end
end

% check if want to plot a subset of domain
[xplt, zplt, bg, u, w] = domainsubset(opt.boxind, x, z, bg, u, w);

% prepare quiver
scl = 2*prctile(abs(w), 90, [2,3,4]);
[Xquiv, Zquiv, u, w] = quivds(xplt, zplt, opt, scl, u, w);

% get axis limits
if opt.uaxes, [climits, cblimits] = uniformaxislimits(opt.Nstd, varname{1}, bg(:,:,:,opt.iPlt), fp); end


% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',18};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',14};

% initialize figure and axes
fig = figure;
colormap(ocean);
set(fig,'Name',[RunID ', ' varname{1} ' field, ' varname{2:3} ' vectors']);
hAx = default2dpanels(NPHS, opt.Nplt, 'aspectratio', length(xplt)/length(zplt), 'top', 2.00);

% plot panels
for m = 1:opt.Nplt
    
    for iphs=1:NPHS
        axes(hAx((iphs-1)*opt.Nplt+m));
        
        % plot background
        imagesc(xplt,zplt,squeeze(bg(iphs,:,:,opt.iPlt(m))));
        hold on;
        
        % plot zero line
        if ~isempty(opt.zline)        % plot a horizontal line
            plot(xlim, x(opt.zline)*ones(1,2), 'LineWidth',1,'Color',0.5*ones(1,3));
        end
        
        % plot u, w quivers
        hquiv = quiver(Xquiv,Zquiv,squeeze(u(iphs,:,:,opt.iPlt(m))),squeeze(w(iphs,:,:,opt.iPlt(m))),'k','LineWidth',1);
        if opt.uquiv, hquiv.AutoScale = 'off'; end
        hold off;
        
        axis xy equal tight;
        xlim([min(xplt), max(xplt)]); 
        ylim([min(zplt), max(zplt)]); ylabel('depth [m]');
        hAx((iphs-1)*opt.Nplt+m).YAxis.Exponent = 0;
        set(gca,TL{:},TS{:}); set(gca,'XTickLabel',[]);
        if m>1, set(gca,'YTickLabel',[]); end
        
        if (opt.uaxes), caxis(climits(iphs, :)); end
        
        cb = colorbar; set(cb,TL{:},TS{:});
        cb.Limits = cblimits(iphs, :);

        title(['phase ' num2str(iphs) ', t = ' num2str(t(opt.iPlt(m)),'%.1e') ' s'],TX{:},FS{:});
        
    end
end


% plot supertitle to tell us what is plotted
supertitle = ['bg: $' varname{1} '$; vectors: $' varname{2} ', ' varname{3}, '$'];
ha = annotation('textbox','Position',[0.5,0.98,0.05,0.02],'String', supertitle, ...
    'HorizontalAlignment','center','VerticalAlignment','top',...
    'FitBoxToText','on','EdgeColor','none','FontSize',22,TX{:});



if opt.save
    if ~isempty(opt.fname)
        fname = opt.fname; 
    else 
        fname = strjoin(varname,'');
    end
    figname = [folder, RunID '/' RunID '_quivfieldt_' fname];
    SaveFigure(figname, fig);
end


end



function [opt] = defopts (Nf, varargin)

% default opts
opt.fname  = '';        % specified filename
opt.save   = 0;         % whether to save

opt.Nquiv  = 21;        % number of quiver arrows (downsampling)
opt.bcind  = 0.05;      % exclude arrows from boundary
opt.uquiv  = 1;         % whether to plot same length quiver for all plts

opt.zline  = [];        % whether to plot a horizontal zero line
opt.xdsc   = 0;         % whether to divide by initial max dsc
opt.rmadv  = 0;         % whether to remove box advection in periodic BC
opt.boxind = [];        % whether to plot a subset of domain 

opt.Nplt   = 5;         % number of panels to plot
opt.iPlt   = [];        % which panels to plot
opt.uaxes  = 1;         % whether to have uniform axes for all panels
opt.Nstd   = 3;         % number of stds from mean for axis limits


% allow structure alteration
args = reshape(varargin, 2, []);
for ia = 1:size(args,2)
    opt.(args{1,ia}) = args{2,ia};
end

% check form of Nplt
if isempty(opt.iPlt)
    if length(opt.Nplt)==1
        if opt.Nplt==1
            opt.iPlt = 1:Nf;
        else
            opt.iPlt = unique(round(linspace(1,Nf,opt.Nplt)));
        end
    else
        opt.iPlt   = opt.Nplt;
    end
end
opt.Nplt = length(opt.iPlt);



end






