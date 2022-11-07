function [fig,figname,varmat,x,z] = PlotFieldwTime (folder, RunID, varname, varmat, t, x, z, varargin)
%
% fig = PlotFieldwTime (folder, RunID, varname, varmat, varargin)
% use this to plot a field in time. Full 2D field.
%
% examples
% fig = PlotFieldwTime(folder, RunID, 'f');
% fig = PlotFieldwTime(folder, RunID, 'f', [], 'Nstd', 5);
% fig = PlotFieldwTime(folder, RunID, 'rhobar', rhobar, 'Nstd', 5);
%
%
% INPUTS
% folder    folder name where output folder is stored
% RunID     name of the run so that the total path to a mat file is
%               [folder  RunID / RunID_0.mat]
% varname   character vector of the name of variable to plot
%           NB if you want certain axis limit options, varname has to be precise
% varmat    matrix of the variable to plot
%               if mising or empty, load the variable of varname from file
% t, x      optional arguments (can be empty) if you already have varmat
%           populated so we don't have to load mat files
% varargin  plotting options (see defopts)
%
% 
% OUTPUTS
% fig       figure handle
% varmat    matrix of variables plotted [NPHS x Nx x Nx]
% x         x positions
% 
% 
% DEFAULT OPTIONS
% opt.fname  = '';        % specified filename
% opt.save   = 0;         % whether to save
% opt.Nplt   = 5;         % number of panels to plot
% opt.iPlt   = [];        % which panels to plot
% 
% opt.zline  = [];        % whether to plot a horizontal zero line
% opt.xdsc   = 0;         % whether to divide by initial max dsc
% opt.boxind = [];        % whether to plot a subset of domain
% 
% opt.uaxes  = 1;         % whether to have uniform axes for all panels
% opt.Nstd   = 3;         % number of stds from mean for axis limits
% 
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

% load varmat to plot
if isempty(varmat)
    % if varmat undefined and you want to load all the variables from file
    [t, x, z, varmat] = ExtractFieldwTime(folder, RunID, {varname});
end
NPHS = size(varmat, 1);
Nf   = size(varmat, 4);

% get plotting options
opt = defopts(Nf, varargin{:});

if (opt.xdsc)
    load(fp, 'delta0');
    x = x./max(delta0(:));
    z = z./max(delta0(:));
end


% check if want to plot a subset of domain
[xplt, zplt, varmat] = domainsubset(opt.boxind, x, z, varmat);

% get axis limits
if opt.uaxes, climits = uniformaxislimits(opt.Nstd, varname, varmat(:,:,:,opt.iPlt), fp); end


% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',18};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',14};

% initialize figure and axes
fig = figure;
colormap(ocean);
set(fig,'Name',[RunID ', ' varname ' field']);
hAx = default2dpanels(NPHS, opt.Nplt, 'aspectratio', length(xplt)/length(zplt));


% plot panels
for m = 1:opt.Nplt
    
    for iphs=1:NPHS
        axes(hAx((iphs-1)*opt.Nplt+m));
        
        imagesc(xplt, zplt, squeeze(varmat(iphs,:,:,opt.iPlt(m))));
        
        if ~isempty(opt.zline)        % plot a horizontal line
            hold on;
            plot(xlim, x(opt.zline)*ones(1,2), 'LineWidth',1,'Color',0.5*ones(1,3));
            hold off;
        end
        
        axis xy equal tight;
        if (opt.uaxes), caxis(climits(iphs, :)); end
        cb = colorbar; set(cb,TL{:},TS{:});
        
        hAx((iphs-1)*opt.Nplt+m).YAxis.Exponent = 0;
        set(gca,TL{:},TS{:}); set(gca,'XTickLabel',[]);
        
        if m>1
            set(gca,'YTickLabel',[]); 
        else
            ylabel('depth [m]');
        end
        
        if NPHS==1, phssuperscript = '';
        else, phssuperscript = ['^',num2str(iphs)];
        end
        title(['$' varname phssuperscript,'$, t = ' num2str(t(opt.iPlt(m)),'%.1e') ' s'],TX{:},FS{:});
        
    end
end

if opt.save
    if ~isempty(opt.fname)
        fname = opt.fname;
    else
        fname = varname;
    end
    
    figname = [folder, RunID '/' RunID '_fieldt_' fname];
    SaveFigure(figname, fig);
end



end



function [opt] = defopts (Nf, varargin)

% default opts
opt.fname  = '';        % specified filename
opt.save   = 0;         % whether to save
opt.Nplt   = 5;         % number of panels to plot
opt.iPlt   = [];        % which panels to plot

opt.zline  = [];        % whether to plot a horizontal zero line
opt.xdsc   = 0;         % whether to divide by initial max dsc
opt.boxind = [];        % whether to plot a subset of domain

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
















