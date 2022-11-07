function [t, x, z, varargout] = LoadPlotVars (folder, RunID, varname, varmat, t, x, z)

if ~iscell(varname), varname = {varname}; end
if ~iscell(varmat ), varmat  = {varmat};  end

% load varmat to plot
if isempty(varmat)
    % if varmat undefined and you want to load all the variables from file
    [t, x, z, varmat] = GetVars(folder, RunID, varname);
    
else
    % check which elements of varmat are defined
    vdef = false(length(varmat),1);
    for vi = 1:length(varmat)
        vdef(vi) = isempty(varmat{vi});
    end
    
    if sum(vdef)>0
        % if some but not all of varmat are defined
        [t, x, z, vtmp] = GetVars(folder, RunID, varname(vdef));
        varmat(vdef) = vtmp;
    else
        % if all of varmat are defined, just load t,x
        [t, x, z] = GetVars (folder, RunID, {});
    end
end

if nargout==4 && length(varname)>1
    varargout = {varmat};
else
    varargout = varmat;
end
end

function [t, x, z, varmat] = GetVars (folder, RunID, varname)

% get files from simulations
[~,fn] = GetOutputMatFiles(folder, RunID);

Nf   = length(fn);          % number of files
Nvar = length(varname);    % number of variables

% initialise
varmat = cell(Nvar,1);
t      = zeros(length(fn),1);

%  assign variables from file
for fi = 1:Nf
    tmp = load(fn{fi}, 'time', 'x', 'z', varname{:});
    
    t(fi) = tmp.time;
    x     = tmp.x;
    z     = tmp.x;
    if isfield(tmp,'z'), z = tmp.z;    end
    
    for vi = 1:Nvar
        switch varname{vi}
            case 'delta'
                % this has a special matrix size NPHS x NPHS x Nz x Nx
                varmat{vi}(:,:,:,:,fi) = tmp.(varname{vi});
            otherwise
                % usual size NPHS x Nz x Nx
                varmat{vi}(:,:,:,fi) = tmp.(varname{vi});
        end
    end
    
end

end