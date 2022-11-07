function [f, t, z] = TranslateProfile (folder, RunID, wref, varargin)
% 
% [f, t, z] = TranslateProfile (folder, RunID, varargin)
% 
% collect profiles of phase fraction and translate by the max reference
% velocity. This assumes periodic boundary conditions. From a simulation
% specified by [folder, RunID]
% 
% INPUTS
% folder  	folder name where the simulation is stored
% RunID 	name of simulation
% varargin  plotting and saving options
% 
% OUTPUTS
% f         translated phase fraction profiles [NPHS x N x Nt]
% t         simulation times [Nt x 1]
% z         vertical positions [1 x N]
% 
% YQW, 23 March 2021

opt = defopts(varargin{:});

% get the phase fraction profiles
[t, zp, v] = GetVertProfiles(folder, RunID, {'f'});

% collect output files
[fp,fn] = GetOutputMatFiles(folder, RunID);
Nf = length(fn);

% get phase names and other params
load(fp,'N','NPHS','PHS','f0');
if ~exist('PHS', 'var'), PHS = strcat({'f'}, num2str((1:NPHS)', '%d')); end

if isempty(wref)
    % collect the maximum value of the reference velocity
    wref = nan(Nf,1);
    for fi = 1:Nf
        load(fn{fi}, 'wstar');
        [~, indstar] = max(abs(wstar), [], [2,3], 'linear');
        wref(fi)   = wstar(indstar);
    end
elseif length(wref)==1
    wref = wref.*ones(Nf,1);
end

% cell offset
z  = nan(1   , N, Nf); z(:,:,1) = zp{1};
f  = nan(NPHS, N, Nf);

% reassign cell values by the cell offset
for fi = 2:Nf
    z(:,:,fi) = z(:,:,fi-1) - wref(fi)*(t(fi)-t(fi-1));
end

f = v{1};



if opt.plot
    % plot outputs
    figure;
    
    ind = [round(opt.wrapind*N):N, 1:(round(opt.wrapind*N)-1)];
    
    set(gcf,'defaultaxescolororder', copper(Nf+1));
    set(gcf,'Position', [500,100,200*NPHS,400]);
    tiledlayout(1, NPHS, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for fi = 1:NPHS
        nexttile;
        plot(squeeze(f(fi,:,:)), squeeze(z(1,:,:)));
        axis manual; 
        hold on; plot(xlim, zeros(1,2), 'k:', 'LineWidth', 1); hold off;
        title([PHS{fi} '(t=0) = ' num2str(100*f0(fi), '%.0f'), '\%'])
        if fi==1,    ylabel(['$z$ [m]']); end
        xlabel(['$\phi^{' PHS{fi} '} (t)$']);
    end
    
    if opt.save
        SaveFigure([folder RunID,'/',RunID,'_fprofiletranslated']);
    end
end

end

function opt = defopts (varargin)

opt.plot = true;    % plot?
opt.save = false;    % save plot?
opt.wrapind = 0.63; % how to plot the periodic wrapping

% allow structure alteration
args = reshape(varargin, 2, []);
for ia = 1:size(args,2)
    opt.(args{1,ia}) = args{2,ia};
end

end