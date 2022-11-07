function [t, x, z, varargout] = ExtractFieldwTime (folder, RunID, varnames, ti)
%
% [t, x, varmat] = ExtractFieldwTime(folder, RunID, varnames)
%
% example:
% [t, x, varmat] = ExtractFieldwTime('../out/', 'olv10_plg10_bas80', {'f','wsegr'})
%
% extracts variable fields from a run specified by [folder, RunID]
%
% INPUTS
% folder  	folder name where the simulation is stored
% RunID 	name of simulation
% varnames 	cell vector of variable names [Nvar x 1]
% ti        indices to extract
%
% OUTPUTS
% t         time of simulation [Nf x 1]
% x         x positions of simulation [Nx x 1]
% varmat    cell of variables {Nvar x 1}
%
% YQW, 4 May 2021

% get files from simulations
[~, fn] = GetOutputMatFiles(folder, RunID);

Nf   = length(fn);          % number of files
Nvar = length(varnames);    % number of variables

% which indices to extract? just check some inputs
if nargin<4   , ti = 1:Nf; end
if isempty(ti), ti = 1:Nf; end

% initialise
varargout = cell(Nvar,1);
t         = zeros(length(ti),1);

%  assign variables from file
for i = 1:length(ti)
    fi  = ti(i);
    tmp = load(fn{fi}, 'time','x','z',varnames{:});
    
    t(i) = tmp.time;
    x    = tmp.x;
    z    = tmp.x;
    if isfield(tmp,'z'), z = tmp.z; end
    
    for vi = 1:Nvar
        switch varnames{vi}
            case 'delta'
                % this has a special matrix size NPHS x NPHS x Nz x Nx
                varargout{vi}(:,:,:,:,i) = tmp.(varnames{vi});
            otherwise
                % usual size NPHS x Nz x Nx
                varargout{vi}(:,:,:,i) = tmp.(varnames{vi});
        end
    end
    
end

if nargout==4 && length(varnames)>1, varargout = {varargout}; end


end