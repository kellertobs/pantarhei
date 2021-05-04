function [] = postplot (folder, RunID)

[fn, fp] = GetOutputMatFiles(folder, RunID);

load(fp);
load('../src/ocean.mat');

% initialise indexing for boundary condition stencils
if     strcmp(BC,'periodic'); ic = [N,1:N,1]; im = [N,1:N]; ip = [1:N,1];
elseif strcmp(BC,'open') || strcmp(BC,'closed'); ic = [1,1:N,N]; im = [1,1:N]; ip = [1:N,N]; end

nop = abs(nop);
dt  = 1;

rho = rho0.* ones(NPHS,N,N);

for fi = 1:length(fn)
    step = (fi-1)*nop;
        
    load(fn{fi});
    
    if fi==1
        fo = f; 
    else
        dt = time-to;
    end

    run('../src/closures.m');
    run('../src/output');
    
    to = time;
    fo = f;
    
end


end

function [fn, fp, NPHS] = GetOutputMatFiles (folder, RunID)

% get output mat files and sort in order
outdir	= [folder RunID '/'];
f    	= dir([outdir '*.mat']);

fname = strcat(outdir, {f.name}');
fp    = char(fname( contains(fname, '_par')));
fn    =      fname(~contains(fname, '_par'));

ind = zeros(length(fn),1);
for fi = 1:length(fn)
    tmp     = strsplit(fn{fi},'_');
    tmp2    = strsplit(tmp{end}, '.');
    ind(fi) = str2double(tmp2{1});
end

% sort by time because of numbering issues
[~,order] = sort(ind,'ascend');
fn        = fn(order);

% check the number of phases
if nargout>2, load(fp, 'NPHS'); end
end