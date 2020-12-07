function [] = MMS_coeff_syms ()
% symbolic math to get the partial derivatives of flux and transfer coeffs
% wrt phase fraction
% 
% YQW, 26 November 2020

%----------------- define variables----------------------------------------
syms f [3,1]            % phase fractions
syms eta0 d0 [3,1]      % pure-phase properties
syms A B C [3,3]        % permission weight parameters

VarOrder = {'f1','f2','f3','d01','d02','d03','eta01','eta02','eta03'...
    'A1_1','A1_2','A1_3','A2_1','A2_2','A2_3','A3_1','A3_2','A3_3'...
    'B1_1','B1_2','B1_3','B2_1','B2_2','B2_3','B3_1','B3_2','B3_3',...
    'C1_1','C1_2','C1_3','C2_1','C2_2','C2_3','C3_1','C3_2','C3_3'};

%--------------------------------------------------------------------------
NPHS = length(eta0);

% pure-phase diffusion parameters
kv = eta0;
kf = d0.^2./eta0;

% get diffusivity contrasts
Mv = kv.'./kv;          % momentum diffusivity ratios
Mf = kf.'./kf;          % volume diffusivity ratios

% permission weights
F     = permute(repmat(f,1,1,1,NPHS),[4,1,2,3]);
Sperm = (F./B).^(1./C);
Sperm = Sperm./sum(Sperm,2);
Xperm = sum(A.*Sperm,2).*F + (1-sum(A.*Sperm,2)).*Sperm;

% permissions
thtv = prod(Mv.^Xperm,2);
thtf = prod(Mf.^Xperm,2);

% flux coefficients
Kv = f.*kv.*thtv;
Kf = f.*kf.*thtf;

% % transfer coefficients
% Cv = Kv.*(1-f)./d0.^2;
% Cf = Kf.*(1-f)./d0.^2;

% need derivative for this omega
omKf = Kf./sum(Kf,1);

% get derivatives wrt f
dKvdf1 = diff(Kv,f1);
dKvdf2 = diff(Kv,f2);
dKvdf3 = diff(Kv,f3);

dKfdf1 = diff(Kf,f1);
dKfdf2 = diff(Kf,f2);
dKfdf3 = diff(Kf,f3);

% dCvdf1 = diff(Cv,f1);
% dCvdf2 = diff(Cv,f2);
% dCvdf3 = diff(Cv,f3);
% 
% dCfdf1 = diff(Cf,f1);
% dCfdf2 = diff(Cf,f2);
% dCfdf3 = diff(Cf,f3);

domKfdf1 = diff(omKf, f1);
domKfdf2 = diff(omKf, f2);
domKfdf3 = diff(omKf, f3);

fprintf(1,'Writing derivatives to file...\n');
tic; 
matlabFunction(...
    dKvdf1,dKvdf2,dKvdf3,   dKfdf1,dKfdf2,dKfdf3,...
    domKfdf1, domKfdf2, domKfdf3,...
    'File','SymsFuncs/coeff_derivs','Vars',VarOrder); 
toc;
fprintf(1,'Done writing derivatives to file...\n\n');


end