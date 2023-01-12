function [fp, fh, fn, ft, NPHS] = GetOutputMatFiles (folder, RunID)
% 
% [fp, fn, ft, NPHS] = GetOutputMatFiles (folder, RunID)
% 
% example: [fp, fn, ft, NPHS] = GetOutputMatFiles('../out/', 'olv20_bas80')
% 
% extracts the names of simulation and parameter files from a simulation
% specified by [folder, RunID]
% 
% INPUTS
% folder  	folder name where the simulation is stored
% RunID 	name of simulation
% 
% OUTPUTS
% fp        string of parameter file name
% fn        cell vector containing all the names of simulations
% ft        cell vector of other mat files
% NPHS      number of phases
% 
% YQW, 19 January 2021
% 

% get output mat files and sort in order
outdir	= [folder RunID '/'];
f    	= dir([outdir '*.mat']);

% separate files into simulation and parameter files
fname = strcat(outdir, {f.name}');
fp    = fname(contains(fname, '_par'  ));    %  parameter file
fh    = fname(contains(fname, '_hist' ));    %    history file
fn    = fname(contains(fname, '_frame'));    % simulation files
ft    = setdiff(fname, [fp; fh; fn]);        %  other mat files 

% return parameter and history file names as a character vector
fp = fp{1};  fh = fh{1};

% if you only want the parameters file, quit
if nargout<=2, return; end

% collect number at end of simulation file name so we can sort in order
ind = zeros(length(fn),1);
for fi = 1:length(fn)
    tmp     = strsplit(fn{fi},'_');
    tmp2    = strsplit(tmp{end}, '.');
    ind(fi) = str2double(tmp2{1});
end

% sort by simulation index
[~,order] = sort(ind,'ascend');
fn        = fn(order);

% return this last file as a character vector if there is only one
if length(ft)==1, ft = ft{1}; end

% check the number of phases
if nargout>4, load(fp, 'NPHS'); end
end
