function [s, e1, e2, qvxzcc] = CalcPrincipalStress (folder, RunID, plot_opt, ti, bgname, bgmat, varargin)
%
% [s, e1, e2] = CalcPrincipalStress (folder, RunID, plot_opt, ti, bgname, bgmat, varargin)
%
% calculates the principal stresses and principal stress directions at each
% grid point (i.e. eigenvalues and eigenvectors of stress matrix).
% 2 principal stresses and 2 directions because 2D.
% easy to calculate analytically.
%
% 
% INPUTS
% folder    folder name where output folder is stored
% RunID     name of the run so that the total path to a mat file is
%               [folder  RunID / RunID_0.mat]
% plot_opt  whether to plot
% ti        which indices to plot. Best to specify else code is slow
% bgname    name of variable plotted in background
% bgmat     matrix of variable plotted in background [NPHS x N x N x Nti]
%
% 
% OUTPUTS
% s         matrix of principal stresses [NPHS x N x N x Nt x 2]
% e1        1st stress direction (corresponding to s1) [NPHS x N x N x Nt x 2]
% e2        2nd stress direction (corresponding to s2) [NPHS x N x N x Nt x 2]
%
% 
% DEFAULT OPTIONS
% opt.fname  = '';        % specified filename
% opt.save   = 0;         % whether to save
% 
% opt.Nquiv  = 21;        % number of quiver arrows (downsampling)
% opt.bcind  = 0.05;      % exclude arrows from boundary
% opt.uquiv  = 1;         % whether to plot same length quiver for all plts
% opt.splt   = 1;         % whether to plot 1 or 2 principal stresses
% 
% opt.zline  = [];        % whether to plot a horizontal zero line
% opt.xdsc   = 0;         % whether to divide by initial max dsc
% opt.boxind = [];        % whether to plot a subset of domain
% 
% opt.Nplt   = 5;         % number of panels to plot
% opt.uaxes  = 1;         % whether to have uniform axes for all panels
% opt.Nstd   = 3;         % number of stds from mean for axis limits
% 
% YQW, 8 June 2021


% check input arguments
if nargin<3, plot_opt = 0;      end
if nargin<4, ti       = [];     end
if nargin<5, bgname   = 'fp';   end

% load components of stress tensor
[t, x, z, varmat] = ExtractFieldwTime(folder, RunID, {'qvxx','qvxz','qvzz','f','p'}, ti);

% take negative of qv to get positive tension stress matrix
% do not remove f x p to show effect of pressure
qvxx = - (varmat{1}); % - varmat{4}.*varmat{5});
qvxz = -  varmat{2};
qvzz = - (varmat{3}); % - varmat{4}.*varmat{5});

% interpolate shear stresses onto pressure grid (cell centers)
qvxzcc = (qvxz(:,1:end-1,1:end-1,:) + qvxz(:,1:end-1,2:end,:) + qvxz(:,2:end,1:end-1,:) + qvxz(:,2:end,2:end,:)) / 4;

% calculate principal stresses
s1 = 0.5*(qvxx+qvzz) + sqrt(0.25*(qvxx-qvzz).^2 + qvxzcc.^2);
s2 = 0.5*(qvxx+qvzz) - sqrt(0.25*(qvxx-qvzz).^2 + qvxzcc.^2);
s  = cat(5, s1, s2);

% calculate first eigenvector (normalised)
e11 = -qvxzcc./(qvxx-s1);
e1  = 1./sqrt(1+e11.^2).*cat(5, e11, ones(size(e11)));

% calculate second eigenvector (normalised)
e21 = -qvxzcc./(qvxx-s2);
e2  = 1./sqrt(1+e21.^2).*cat(5, e21, ones(size(e21)));

if plot_opt    
    switch bgname
        case 's1'   , bgmat = s1;
        case 's2'   , bgmat = s2;
        case 'f'    , bgmat = varmat{4};
        case 'fbas' , bgmat = varmat{4}(2,:,:,:);
        case 'p'    , bgmat = varmat{5};
        case 'fp'   , bgmat = varmat{4}.*varmat{5};
        case 'sdiff', bgmat = s1-s2;
    end
    
    % make sure the correct bgmat indices are used
    if size(bgmat,4)>length(ti), bgmat = bgmat(:,:,:,ti,:); end
    
    plot_pstress(folder, RunID, bgname, bgmat, t, x, z, s, e1, e2, varargin{:});
end

end



function plot_pstress (folder, RunID, bgname, bg, t, x, z, s, e1, e2, varargin)

fp = GetOutputMatFiles(folder, RunID);
Nf = length(t);

% get plotting options
opt = defopts('Nplt', Nf, varargin{:});

% load colormap
load('ocean.mat', 'ocean');

% load plotting variables
[t, x, z, bg] = LoadPlotVars(folder, RunID, {bgname}, bg, t, x, z);

smax = max(abs(s),[],5);
se1  = 1.*e1;
se2  = (1+(s(:,:,:,:,2))./smax)./(1+(s(:,:,:,:,1))./smax).*e2;

if (opt.xdsc)
    load(fp, 'delta0');
    x = x./max(delta0(:));
    z = z./max(delta0(:));
end

% check sizes and make sure bg, s, e are the same size
NPHS = size(se1,1);
if size(bg,1)<NPHS, bg = bg.*ones(NPHS,1); end

% check if want to plot a subset of domain
[xplt, zplt, bg, se1, se2] = domainsubset(opt.boxind, x, z, bg, se1, se2);

% prepare quiver 
[Xquiv, Zquiv, se1, se2] = quivds(xplt, zplt, opt, 2, se1, se2);

% get axis limits
if opt.uaxes, climits = uniformaxislimits(opt.Nstd, bgname, bg, fp); end

% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',18};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',14};

% initialize figure and axes
fig = figure;
colormap(ocean);
set(fig,'Name',[RunID '_stressfield']);

axh = 240; axb =  3; axgh =  3; axt = 40;
axw = 260; axl = 50; axgw = 10; axr = 20; 
fh = axb +     NPHS*axh + (    NPHS-1)*axgh + axt;
fw = axl + opt.Nplt*axw + (opt.Nplt-1)*axgw + axr;
set(fig,'Position',[500,500,fw,fh]);
hAx = tight_subplot(NPHS,opt.Nplt,[axgh/fh,axgw/fw], [axb,axt]/fh, [axl,axr]/fw);

% plot panels
for m = 1:opt.Nplt
    
    for iphs=1:NPHS
        axes(hAx((iphs-1)*opt.Nplt+m));
        
        % plot background
        imagesc(xplt,zplt,squeeze(bg(iphs,:,:,m)));
        hold on;
        
        % plot zero line
        if ~isempty(opt.zline)        % plot a horizontal line
            zloc = x(opt.zline);
            plot(xlim'.*ones(1,length(opt.zline)), zloc.*ones(2,1), 'LineWidth',0.5,'Color',0.5*ones(1,3));
        end
        
        if any(opt.splt == 2)
            hquiv = quiver(Xquiv,Zquiv,squeeze(se2(iphs,:,:,m,1)),  squeeze(se2(iphs,:,:,m,2)),'color', 0.5*ones(1,3),'LineWidth',1);
            hquiv.AutoScale = 'off'; hquiv.ShowArrowHead = 'off';
            hquiv = quiver(Xquiv,Zquiv,-squeeze(se2(iphs,:,:,m,1)),-squeeze(se2(iphs,:,:,m,2)),'color', 0.5*ones(1,3),'LineWidth',1);
            hquiv.AutoScale = 'off'; hquiv.ShowArrowHead = 'off';
        end
        
        if any(opt.splt==1)
            hquiv = quiver(Xquiv,Zquiv,squeeze(se1(iphs,:,:,m,1)),squeeze(se1(iphs,:,:,m,2)),'k','LineWidth',1);
            hquiv.AutoScale = 'off'; hquiv.ShowArrowHead = 'off';
            hquiv = quiver(Xquiv,Zquiv,-squeeze(se1(iphs,:,:,m,1)),-squeeze(se1(iphs,:,:,m,2)),'k','LineWidth',1);
            hquiv.AutoScale = 'off'; hquiv.ShowArrowHead = 'off';
        end

        hold off;
        
        axis xy equal tight;
        xlim([min(xplt), max(xplt)]); ylim([min(zplt), max(zplt)]);
        if (opt.uaxes), caxis(climits(iphs, :)); end
        
        cb = colorbar; set(cb,TL{:},TS{:});
        hAx((iphs-1)*opt.Nplt+m).YAxis.Exponent = 0;
        set(gca,TL{:},TS{:}); set(gca,'XTickLabel',[]);
        if m>1, set(gca,'YTickLabel',[]); end
        
        if iphs==1, title(['t = ' num2str(t(m),'%.1e') ' s'],TX{:},FS{:}); end
        
    end
end



% plot supertitle to tell us what is plotted
supertitle = ['bg: $' bgname '$; vectors: principal stresses'];
ha = annotation('textbox','Position',[0.5,0.98,0.05,0.02],'String', supertitle, ...
    'HorizontalAlignment','center','VerticalAlignment','top',...
    'FitBoxToText','on','EdgeColor','none','FontSize',22,TX{:});



if opt.save
    if ~isempty(opt.fname)
        fname = opt.fname;
    else
        fname = ['s' num2str(1:opt.splt, '%.0f')];
    end
    figname = [folder, RunID '/' RunID '_pstresst_' fname];
    SaveFigure(figname, fig);
end


end







function [opt] = defopts (varargin)

% default opts
opt.fname  = '';        % specified filename
opt.save   = 0;         % whether to save

opt.Nquiv  = 15;        % number of quiver arrows (downsampling)
opt.bcind  = 0.05;      % exclude arrows from boundary
opt.uquiv  = 1;         % whether to plot same length quiver for all plts
opt.splt   = [1,2];         % whether to plot 1 or 2 principal stresses

opt.zline  = [];        % whether to plot a horizontal zero line
opt.xdsc   = 0;         % whether to divide by initial max dsc
opt.boxind = [];        % whether to plot a subset of domain

opt.Nplt   = 1;         % number of panels (depends on input)
opt.uaxes  = 1;         % whether to have uniform axes for all panels
opt.Nstd   = 3;         % number of stds from mean for axis limits


% allow structure alteration
args = reshape(varargin, 2, []);
for ia = 1:size(args,2)
    opt.(args{1,ia}) = args{2,ia};
end


end






