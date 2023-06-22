function  [upd_u, upd_w, upd_p]  =  direct_solve(u, w, p, f, Qvz, Qf, PHS, NUM)

NPHS  =  NUM.NPHS;
Nz    =  NUM.Nz;
h     =  NUM.h;

ip    =  NUM.ipz;
im    =  NUM.imz;

izc   =  1:Nz+1;           % central z-face node index
izp   =  izc([2:end,end]); % forward z-face node field index
izm   =  izc([1,1:end-1]); % backward z-face node field index

icc   =  1:Nz+2;           % central centre node index (ghosted)
icp   =  icc([2:end,end]); % forward centre node index (ghosted)
icm   =  icc([1,1:end-1]); % backward centre node index (ghosted)

ivc   =  izc + ((1:NPHS).'-1).*(Nz+1);  % central phase velocity global index
ivp   =  izp + ((1:NPHS).'-1).*(Nz+1);  % forward phase velocity global index
ivm   =  izm + ((1:NPHS).'-1).*(Nz+1);  % backward phase velocity global index

ipc   =  icc + ((1:NPHS).'-1).*(Nz+2) + NPHS*(Nz+1);  % central phase pressure global index
ipp   =  icp + ((1:NPHS).'-1).*(Nz+2) + NPHS*(Nz+1);  % forward phase pressure global index
ipm   =  icm + ((1:NPHS).'-1).*(Nz+2) + NPHS*(Nz+1);  % backward phase pressure global index

% initialize operator assembly arrays
AA  =  [];
IA  =  [];
JA  =  [];
RR  =  [];
IR  =  [];


% get coefficient closures and constitutive relations
closures; constitutive;

% SEGREGATION EQUATION

% momentum flux divergence (-Div.Kvi.D(vi))
Kvp = 2/3.*Kv(:,ip);    Kvm = 2/3.*Kv(:,im);

a  =  - Kvp     ./h^2;   i = ivc;   j = ivp;   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
a  =  +(Kvp+Kvm)./h^2;   i = ivc;   j = ivc;   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
a  =  -     Kvm ./h^2;   i = ivc;   j = ivm;   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];

% momentum transfer (Cvi.(vi-v*))
a  =  Cvz;   i = ivc;   j = ivc;   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
for n=1:NPHS
    a  =  -Cvz(n,:).*omvz;   i = repmat(ivc(n,:),NPHS,1);   j = ivc;   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
end

% reference pressure gradient (phii.Grad.P*)
for n=1:NPHS
    a  =  fz(n,:).*omfc(:,ip)./h;  i = repmat(ivc(n,:),NPHS,1);  j = ipc(:,2:end-0);   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
    a  = -fz(n,:).*omfc(:,im)./h;  i = repmat(ivc(n,:),NPHS,1);  j = ipc(:,1:end-1);   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
end

% compaction pressure gradient (Grad.phii.Pi)
a  =  f(:,ip)./h;  i = ivc;  j = ipc(:,2:end-0);   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
a  = -f(:,im)./h;  i = ivc;  j = ipc(:,1:end-1);   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
for n=1:NPHS
    a  = -f(n,ip).*omfc(:,ip)./h;  i = repmat(ivc(n,:),NPHS,1);  j = ipc(:,2:end-0);   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
    a  =  f(n,im).*omfc(:,im)./h;  i = repmat(ivc(n,:),NPHS,1);  j = ipc(:,1:end-1);   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
end

% momentum source term (phii.Drhoi.g)
r  =  Qvz;  i = ivc;  RR = [RR;r(:)]; IR = [IR;i(:)];


% COMPACTION EQUATION

% volume flux divergence (-Div.Kphii.(Grad.Pi-Grad.P*))
Kfp  = Kfz (:,2:end);  Kfm  = Kfz (:,1:end-1);
ompp = ompz(:,2:end);  ompm = ompz(:,1:end-1);
omvp = omvz(:,2:end);  omvm = omvz(:,1:end-1);
fp   = fz  (:,2:end);  fm   = fz  (:,1:end-1);

a  =  - Kfp     ./h^2;   i = ipc(:,2:end-1);   j = ipp(:,2:end-1);   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
a  =  +(Kfp+Kfm)./h^2;   i = ipc(:,2:end-1);   j = ipc(:,2:end-1);   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
a  =  -     Kfm ./h^2;   i = ipc(:,2:end-1);   j = ipm(:,2:end-1);   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
for n=1:NPHS
    a  =  + Kfp(n,:).*ompp                ./h^2;   i = repmat(ipc(n,2:end-1),NPHS,1);   j = ipp(:,2:end-1);   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
    a  =  -(Kfp(n,:).*ompp+Kfm(n,:).*ompm)./h^2;   i = repmat(ipc(n,2:end-1),NPHS,1);   j = ipc(:,2:end-1);   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
    a  =  +                Kfm(n,:).*ompm ./h^2;   i = repmat(ipc(n,2:end-1),NPHS,1);   j = ipm(:,2:end-1);   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
end

% volume transfer (Cphii.(Pi-P*))
a  =  Cf;   i = ipc(:,2:end-1);   j = ipc(:,2:end-1);   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
for n=1:NPHS
    a  =  -Cf(n,:).*omfc;   i = repmat(ipc(n,2:end-1),NPHS,1);   j = ipc(:,2:end-1);   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
end

% reference velocity divergence (phii.Div.v*)
for n=1:NPHS
    a  =  f(n,:).*omvp./h;  i = repmat(ipc(n,2:end-1),NPHS,1);  j = ivc(:,2:end-0);   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
    a  = -f(n,:).*omvm./h;  i = repmat(ipc(n,2:end-1),NPHS,1);  j = ivc(:,1:end-1);   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
end

% segregation velocity divergence (Div.phii.(vi-v*))
a  =  fp./h;  i = ipc(:,2:end-1);  j = ivc(:,2:end-0);   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
a  = -fm./h;  i = ipc(:,2:end-1);  j = ivc(:,1:end-1);   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
for n=1:NPHS
    a  = -fp(n,:).*omvp./h;  i = repmat(ipc(n,2:end-1),NPHS,1);  j = ivc(:,2:end-0);   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
    a  =  fm(n,:).*omvm./h;  i = repmat(ipc(n,2:end-1),NPHS,1);  j = ivc(:,1:end-1);   AA = [AA;a(:)]; IA = [IA;i(:)]; JA = [JA;j(:)];
end

% volume source term (0)
r  =  Qf;  i = ipc(:,2:end-1);  RR = [RR;r(:)]; IR = [IR;i(:)];


% ASSEMBLE OPERATOR
% assemble coefficient matrix
L  =  sparse(IA,JA,AA,NPHS*(Nz+1) + NPHS*(Nz+2),NPHS*(Nz+1) + NPHS*(Nz+2));

% assemble right hand side vector
R  =  sparse(IR,1,RR,NPHS*(Nz+1) + NPHS*(Nz+2),1);


% BOUNDARY CONDITIONS
switch NUM.BC{1}
    
    case 'closed'
        
        % vel boundary conditions
        i  =  ivc(:,1);        L(i(:),:)    = 0;
        j  =  ivc(:,1);        L(i(:),j(:)) = speye(length(i(:)));
        r  =  zeros(NPHS,1);   R(i(:))      = r(:);
        
        i  =  ivc(:,end);      L(i(:),:)    = 0;
        j  =  ivc(:,end);      L(i(:),j(:)) = speye(length(i(:)));
        r  =  zeros(NPHS,1);   R(i(:))      = r(:);
        
        % prs boundary conditions        
        i  =  ipc(:,1);        L(i(:),:)    = 0;
        j  =  ipc(:,1);        L(i(:),j(:)) =  speye(length(i(:)));
        j  =  ipc(:,2);        L(i(:),j(:)) = -speye(length(i(:)));
        r  =  zeros(NPHS,1);   R(i(:))      = r(:);

        i  =  ipc(:,end  );    L(i(:),:)    = 0;
        j  =  ipc(:,end  );    L(i(:),j(:)) =  speye(length(i(:)));
        j  =  ipc(:,end-1);    L(i(:),j(:)) = -speye(length(i(:)));
        r  =  zeros(NPHS,1);   R(i(:))      = r(:);

        i  =  ipc(1,Nz/2+1);   L(i(:),:)    = 0;
        j  =  ipc(1,Nz/2+1);   L(i(:),j(:)) = speye(length(i(:)));
        r  =  0;   R(i(:)) = r(:);

    case 'opentop'
        
        % vel boundary conditions
        i  =  ivc(:,1);        L(i(:),:)    = 0;
        j  =  ivc(:,1);        L(i(:),j(:)) = speye(length(i(:)));
        r  =  zeros(NPHS,1);   R(i(:))      = r(:);
        
        % prs boundary conditions        
        i  =  ipc(:,1);        L(i(:),:)    = 0;
        j  =  ipc(:,1);        L(i(:),j(:)) =  speye(length(i(:)));
        j  =  ipc(:,2);        L(i(:),j(:)) = -speye(length(i(:)));
        r  =  zeros(NPHS,1);   R(i(:))      = r(:);

        i  =  ipc(:,end);      L(i(:),:)    = 0;
        j  =  ipc(:,end);      L(i(:),j(:)) = speye(length(i(:)));
        j  =  ipc(:,end-1);    L(i(:),j(:)) = speye(length(i(:)));
        r  =  zeros(NPHS,1);   R(i(:)) = r(:);

    case 'openbot'
        
        % vel boundary conditions
        i  =  ivc(:,end);      L(i(:),:)    = 0;  
        j  =  ivc(:,end);      L(i(:),j(:)) = speye(length(i(:)));
        r  =  zeros(NPHS,1);   R(i(:))      = r(:);
        
        % closed domain prs boundary conditions        
        i  =  ipc(:,end  );    L(i(:),:)    = 0;
        j  =  ipc(:,end  );    L(i(:),j(:)) =  speye(length(i(:)));
        j  =  ipc(:,end-1);    L(i(:),j(:)) = -speye(length(i(:)));
        r  =  zeros(NPHS,1);   R(i(:))      = r(:);

        i  =  ipc(:,1);        L(i(:),:)    = 0;
        j  =  ipc(:,1);        L(i(:),j(:)) = speye(length(i(:)));
        j  =  ipc(:,2);        L(i(:),j(:)) = speye(length(i(:)));
        r  =  zeros(NPHS,1);   R(i(:)) = r(:);
end

% get current solution vector
S  = zeros(size(R));
S(ivc) = w(:);
S(ipc(:,2:end-1)) = p(:);

% get residual vector
F  = L*S-R;

% scale operator
X  =  sqrt(abs(diag(L)));
X  =  diag(sparse(1./(X+1e-16)));

L  =  X*L*X;
F  =  X*F; 

% SOLVE SYSTEM FOR UPDATE
U    =  X * (L\F);

upd_w  =  reshape(full(U(1:NPHS*(Nz+1))),(Nz+1),NPHS).';
upd_p  =  reshape(full(U(NPHS*(Nz+1)+1:end)),(Nz+2),NPHS).'; upd_p = upd_p(:,2:end-1);
upd_u  =  0.*u;

end