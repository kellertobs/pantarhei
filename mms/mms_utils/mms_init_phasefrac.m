
% define x and z positions of field variables
XuGrid = cat(3,X-h/2,X(:,:,end)+h/2);   ZuGrid = cat(3,Z,Z(:,:,end));
XwGrid = cat(2,X,X(:,end,:));           ZwGrid = cat(2,Z-h/2,Z(:,end,:)+h/2);

% generate manufactured solutions
f = MMSsource('Calc_f', 0, X, Z, Tmf(:,1), Xmf(:,1), Zmf(:,1), Amf(:,1), dmf(:,1));
fo = f;  fi = f;  res_f = 0*f;  dtau_f = res_f;

% check segregation compaction length is on order of domain
dsc = SegCompLength(f, eta0, d0, A, B, C, thtlim, cfflim);
