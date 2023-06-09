
% update constitutive relations
constitutive;

% update residual fields
res_u = Div_qvx + Gvx - Qvx;
res_w = Div_qvz + Gvz - Qvz;
res_p = Div_qf  + Gf  - Qf ;

% apply boundary conditions
if strcmp(NUM.BC{1},'closed'), res_w(:,[1,end],:) = 0; end
if strcmp(NUM.BC{2},'closed'), res_u(:,:,[1,end]) = 0; end