function [climits, cblimits] = uniformaxislimits (Nstd, varname, varmat, fp)


% get uniform axes limits?
if isinf(Nstd)
    % limits = range of values
    climits  = [min(varmat, [], [2,3,4]), max(varmat, [], [2,3,4])];
    cblimits = climits;
else
    % limits = mean + some multiple of std
    climits  = mean(varmat, [2,3,4]) + Nstd*std(varmat,0,[2,3,4]).*[-1,1];
    cblimits = climits;
end

% now adjust axes for some specific variables
switch varname
    case {'p','pcmpt','pstar','df','s1','fp'}
        % set c axis to be symmetric about 0
        climits = max(abs(climits),[],2).*[-1,1];
    case {'f'}
        %  set c axis to be symmetric about background phase fraction
        load(fp, 'f0');
        climits  = f0 + max(abs(climits-f0),[],2).*[-1,1];
        cblimits = [max(0,climits(:,1)), min(1,climits(:,2))];
end

end