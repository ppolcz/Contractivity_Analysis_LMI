%% 
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2021. March 02. (2020b)
%  Revised on 2024. November 20. (2023a)

clear all

%%
% Automatically generated stuff

G_reset(01)
%       │└─ verbosity (0:false, 1:true, 2:do not change)
%       └── scope depth (0:reset, 1:reset if not set)

He = he;
I = @(varargin) eye(varargin{:});
O = @(varargin) zeros(varargin{:});

P_init(12.02)

echo off

%% Select a traffic flow reaction model

[o,r] = TRM_definitions(4,5,'Visualize',2,"Layout","force","View",[35,90], ...
    "SubDIR","TRM_4D_v5");

%%

x = o.x_sym;
x_lfr = o.x_lfr;
xi_lfr = o.xi_lfr;

fi_lfr = o.fi_lfr;

f = o.f_plfr;
A = o.A_plfr;

n = numel(x);

%% Select the structure of PI

PI0 = A.generatePI;
[~,PI,~,~] = P_mingen_for_LFR(PI0,'sym',1);

% Store selection
r.PI = struct('M',PI.M,'blk',PI.blk);


PI_sym = sym(PI);

Mtp = subs(PI_sym,x,ones(n,1));
Mtp(Mtp == 0) = 1;
PI_nice = PI_sym ./ Mtp;

%% Construct further objects 
% Automatically generated from the selected PI

dPI = PI.diff(xi_lfr,fi_lfr);
m = size(PI,1);

PId0 = plfr([ PI.lfr ; PI.lfr*A.lfr + dPI.lfr ]);
J = [ I(m) O(m) ];
Jd = [ O(m) I(m) ];

PId1 = PId0.generatePI;
[Sd,PId,iSd,Kerd] = P_mingen_for_LFR(PId1);

Ia = PI(zeros(n,1))';
Id = PId(zeros(n,1))';

Ed = J*[ PId0.A PId0.B ]*Sd;
Ad = Jd*[ PId0.A PId0.B ]*Sd;

N = P_affine_annihilator_for_LFR(PI,x_lfr,'sym',1);
Nd = P_affine_annihilator_for_LFR(PId,x_lfr);

[s,m] = size(N);
[sd,md] = size(Nd);

if o.Check
    pcz_fhzero_report(N*PI,x,1e-7, 'N(x)*pi(x) = 0');
    pcz_fhzero_report(Nd*PId,x,1e-7, 'Nd(x)*pid(x) = 0');
    pcz_fhzero_report(dPI + PI*A - Ad*PId,x,1e-7, 'Ad * PI(x) = dot PI(x) + PI(x)*A(x)');
    pcz_fhzero_report(Ed*PId - PI,x,1e-7, 'Ed * PId(x) = PI(x)');
    pcz_fhzero_report(Ia*PI - eye(n),x,1e-7, 'Ia * PI(x) = In');
    pcz_fhzero_report(Id*PId - eye(n),x,1e-7, 'Id * PId(x) = In');
end

%%

r = load('/home/ppolcz/Dropbox/Peti/PhD/.Dokumentaciok/37_TRM_Contraction/results/Elso_probalkozasok/trm_dim4_v5_m14_sp0082.mat');

% Sparsity pattern of P is given manually.
Idx0 = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196];

Eps = 0.0;
Clip = [Eps 1-Eps];
x_lim = ones(n,1) * Clip;

X_v = P_ndnorms_of_X(x_lim);

beta = sdpvar;
L = sdpvar(m,s,'full');
Ld = sdpvar(md,sd,'full');
P = sdpvar(m);
P(Idx0) = 0;

CONS = [];

betaM = 0.0001;
for i = 1:size(X_v,1)
    xi = X_v(i,:)';
    
    CONS = [ CONS
        P + He( L*N(xi) ) - 1e-3*(Ia'*Ia) >= 0
        He( Ed'*P*Ad + Ld*Nd(xi) ) + betaM*Ed'*P*Ed <= 0
        ];
end

opts = sdpsettings('solver','sedumi');
sol = optimize(CONS,[],opts);
pcz_feasible(sol,CONS)

r.P = double(P);
r.L = double(L);
r.Ld = double(Ld);
r.sol = sol;
display(r.P,'P')

r.P_zero_in = find(r.P == 0)';

% Save results
if ~startsWith(sol.info,'Infeasible')
    spec = sprintf("_m%d_sp%04g",m,1000-round(numel(r.P_zero_in)/m^2*1000));
    matname = fullfile(o.DIR,o.SubDIR,o.ModelName+spec+".mat");
    save(matname,'-struct','r');
end

%%

PI_sym = sym(PI);

Mtp = subs(PI_sym,x,ones(n,1));
Mtp(Mtp == 0) = 1;


PI_nice = PI_sym ./ Mtp;

P_sym = pcz_sym_symmetric('a',m);
P_sym(abs(r.P) < eps) = 0;
P_vars = symvar(P_sym);

P_str = subs(P_sym,P_vars,sym('p',size(P_vars)));

% Display a parameterized structure of the metric
M = collect(PI_nice.' * P_str * PI_nice,x)

% Display the actual value of the metric (with a given precision)
M = vpa(collect(PI_sym.' * r.P * PI_sym,x),3)


M = expand(PI_sym.' * r.P * PI_sym);

display(vpa(diag(M),3),'Diagonal elements of M')

m_fh = {
    matlabFunction(M(1,1),'Vars',x)
    matlabFunction(M(2,2),'Vars',x)
    matlabFunction(M(3,3),'Vars',x)
    matlabFunction(M(4,4),'Vars',x)
    };

%%

M_fh = matlabFunction(M,'vars',{x});

dM = diff(M,x(1))*o.f_sym(1);
for dim = 2:n
    dM = dM + diff(M,x(dim))*o.f_sym(dim);
end

Ml = PI.lfr' * r.P * PI.lfr;
Q = o.A_sym'*M + M*o.A_sym + dM + betaM*M;
Q2 = sym(PId)' * ( He( Ed'*r.P*Ad ) + betaM*Ed'*r.P*Ed ) * sym(PId);
Q3 = plfr(PId.lfr' * ( He( Ed'*r.P*Ad ) + betaM*Ed'*r.P*Ed ) * PId.lfr);

Q_fh = matlabFunction(Q,'vars',{x});
Q2_fh = matlabFunction(Q2,'vars',{x});


Nr_Samples = 100;
xxr = [ X_v' , rand(n,Nr_Samples).*diff(x_lim,[],2) + x_lim(:,1) ];

printer = @(val) sprintf(strjoin(repmat({'%g'},[1,n]),','),val);

Nem_teljesult = 0;
for x_val = xxr
    lambda = eig(Q_fh(x_val));

    nrm = norm(Q_fh(x_val) - Q3(x_val),'fro');
    % nrm = max(max(abs(Q_fh(x_val) - Q2_fh(x_val))));
    
    if any(real(lambda) >= 0)
        Nem_teljesult = Nem_teljesult + 1;
        fprintf('In (%s), the eigenvalues are: \n',printer(x_val));
        disp(lambda)
    end

    if nrm >= 1e-5
        fprintf('In (%s), the difference is: %g.\n',printer(x_val),nrm);
    end
    
end

if Nem_teljesult == 0
    pcz_OK_FAILED(true,'Lyapunov''s inequality is satisfied in all the %d points including the cornern points within the domain:\n',size(xxr,2))
else
    pcz_OK_FAILED(false,'In %d points out of %d, Lyapunov''s inequality is not satisfied. The domain:\n',Nem_teljesult,size(xxr,2));
end

