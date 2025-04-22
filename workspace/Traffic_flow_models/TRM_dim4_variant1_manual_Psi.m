%% 
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2021. March 02. (2020b)
%  Revised on 2024. November 20. (2023a)

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

%% model

%%%
% A few TRM models are implemented in the following function. We select one of these
% models.
[o,r] = TRM_definitions(3,1,'Visualize',2,"Layout","force","View",[45,90],"SubDIR","TRM_3D_var1");

x = o.x_sym;
x_lfr = o.x_lfr;
xi_lfr = o.xi_lfr;

fi_lfr = o.fi_lfr;

f = o.f_plfr;
A = o.A_plfr;

n = numel(x);

% The studied region of the state-space is the whole admissible region.
x_lim = [
    zeros(n,1) , o.C
    ];

%%

x1 = o.xi_lfr{1};
x2 = o.xi_lfr{2};
x3 = o.xi_lfr{3};

CASE = 1;
switch CASE
    case 1 % ------------------------------------------------
        Pi_lfr = [
            1
            x1
            x2
            x3
            ];
        m1 = size(Pi_lfr,1);
        Z = lfr(zeros(m1,1));
        PI_lfr = [
            Pi_lfr , Z      , Z
            Z      , Pi_lfr , Z
            Z      , Z      , Pi_lfr
            ];

    case 2 % ------------------------------------------------

        Pi = @(x) [
            1
            x
            % x^2
            % x^3
            ];
        m1 = size(Pi(0),1);
        Z = lfr(zeros(m1,1));
        PI_lfr = [
            Pi(x1) , Z      , Z
            Z      , Pi(x2) , Z
            Z      , Z      , Pi(x3)
            ];

    case 3 % ------------------------------------------------

        PI_lfr = [
                 1,          0,          0
                 0,          1,          0
                 0,          0,          1
                x2,         x1,          0
                x3,          0,         x1
                 0,         x3,         x2
              % 2*x1,          0,          0
              %    0,       2*x2,          0
              %    0,          0,       2*x3
           %  3*x1^2,          0,          0
           % 2*x1*x2,       x1^2,          0
           %    x2^2,    2*x1*x2,          0
           %       0,     3*x2^2,          0
           % 2*x1*x3,          0,       x1^2
           %   x2*x3,      x1*x3,      x1*x2
           %       0,    2*x2*x3,       x2^2
           %    x3^2,          0,    2*x1*x3
           %       0,       x3^2,    2*x2*x3
           %       0,          0,     3*x3^2
        %     4*x1^3,          0,          0
        %  3*x1^2*x2,       x1^3,          0
        %  2*x1*x2^2,  2*x1^2*x2,          0
        %       x2^3,  3*x1*x2^2,          0
        %          0,     4*x2^3,          0
        %  3*x1^2*x3,          0,       x1^3
        % 2*x1*x2*x3,    x1^2*x3,    x1^2*x2
        %    x2^2*x3, 2*x1*x2*x3,    x1*x2^2
        %          0,  3*x2^2*x3,       x2^3
        %  2*x1*x3^2,          0,  2*x1^2*x3
        %    x2*x3^2,    x1*x3^2, 2*x1*x2*x3
        %          0,  2*x2*x3^2,  2*x2^2*x3
        %       x3^3,          0,  3*x1*x3^2
        %          0,       x3^3,  3*x2*x3^2
        %          0,          0,     4*x3^3
         ];

    otherwise % ---------------------------------------------

end

PI = plfr(PI_lfr);
dPI = PI.diff(o.xi_lfr,o.fi_lfr);
m = size(PI_lfr,1);

PId0 = plfr([ PI_lfr ; PI_lfr*A.lfr + dPI.lfr ]);
J = [ I(m) O(m) ];
Jd = [ O(m) I(m) ];

PId1 = PId0.generatePI;
[Sd,PId,iSd,Kerd] = P_mingen_for_LFR(PId1);

Ia = PI(zeros(n,1))';
Id = PId(zeros(n,1))';

Ed = J*[ PId0.A PId0.B ]*Sd;
Ad = Jd*[ PId0.A PId0.B ]*Sd;

N = P_affine_annihilator_for_LFR(PI,x_lfr,'sym',1);
Nd = P_affine_annihilator_for_LFR(PId,x_lfr,'sym',1);

[s,m] = size(N);
[sd,md] = size(Nd);

pcz_fhzero_report(N*PI,x,1e-7, 'N(x)*pi(x) = 0');
pcz_fhzero_report(Nd*PId,x,1e-7, 'Nd(x)*pid(x) = 0');
pcz_fhzero_report(dPI + PI*A - Ad*PId,x,1e-7, 'Ad * PI(x) = dot PI(x) + PI(x)*A(x)');
pcz_fhzero_report(Ed*PId - PI,x,1e-7, 'Ed * PId(x) = PI(x)');
pcz_fhzero_report(Ia*PI - eye(3),x,1e-7, 'Ia * PI(x) = In');
pcz_fhzero_report(Id*PId - eye(3),x,1e-7, 'Id * PId(x) = In');

%%

X_v = P_ndnorms_of_X(x_lim);

switch CASE
    case {1,2}

        affineBlk_Impl = @(p) [sdpvar p' ; p zeros(m1-1)];
        affineBlk = @() affineBlk_Impl(sdpvar(m1-1,1));

        multiAffineBlk_Impl = @(c,p,Q) [c p' ; p Q-diag(diag(Q))];
        multiAffineBlk = @() multiAffineBlk_Impl(sdpvar,sdpvar(m1-1,1),sdpvar(m1-1));

        P = diag(sdpvar(m,1));
        % P = blkdiag(sdpvar(m1),sdpvar(m1),sdpvar(m1));
        % P = blkdiag(affineBlk(),affineBlk(),affineBlk());
        % P = blkdiag(multiAffineBlk(),multiAffineBlk(),multiAffineBlk());
    case 3
        P = sdpvar(m);
end
beta = sdpvar;
L = sdpvar(m,s,'full');
Ld = sdpvar(md,sd,'full');

CONS = [];

betaM = 0.5;
for i = 1:size(X_v,1)
    xi = X_v(i,:)';
    
    CONS = [ CONS
        P + He( L*N(xi) ) - 1e-3*(Ia'*Ia) >= 0
        % He( Ed'*P*Ad + Ld*Nd(xi) ) + 1000*(Id'*Id) <= 0
        He( Ed'*P*Ad + Ld*Nd(xi) ) + betaM*Ed'*P*Ed <= 0
        ];
end

sol = optimize(CONS,[],sdpsettings('solver','sedumi')); % or use 'mosek'
pcz_feasible(sol,CONS)

L = double(L);
Ld = double(Ld);
P = double(P)

%%

PI_sym = sym(PI);
M = expand(PI_sym.' * P * PI_sym);

display(vpa(diag(M),3),'Diagonal elements of M')

m_fh = {
    matlabFunction(M(1,1),'Vars',x)
    matlabFunction(M(2,2),'Vars',x)
    matlabFunction(M(3,3),'Vars',x)
    };

fig = figure(132); 
delete(fig.Children)

for dim = 1:3

    ij = setdiff(1:3,dim);

    c = linspace(x_lim(dim,1),x_lim(dim,2),100)';
    [a,b] = ndgrid( ...
        linspace(x_lim(ij(1),1),x_lim(ij(1),2),7), ...
        linspace(x_lim(ij(2),1),x_lim(ij(2),2),7));
    ab = cellfun(@(u) {c*0+u}, num2cell([a(:),b(:)],2));
    abc = cellfun(@(u) { num2cell([u(:,1:dim-1) c u(:,dim:end) ],1) },ab);
    
    nexttile; hold on;
    for Idx = 1:numel(ab)
        plot(c,m_fh{dim}(abc{Idx}{:}));
    end

end

%%

M_fh = matlabFunction(M,'vars',{x});

dM = diff(M,x(1))*o.f_sym(1) + diff(M,x(2))*o.f_sym(2) + diff(M,x(3))*o.f_sym(3);

Q = o.A_sym'*M + M*o.A_sym + dM;
Q_fh = matlabFunction(Q,'vars',{x});


Nr_Samples = 100;
xxr = [ X_v' , rand(3,Nr_Samples).*diff(x_lim,[],2) + x_lim(:,1) ];

for x_val = xxr
    lambda = eig(Q_fh(x_val));

    if any(real(lambda) >= -0.1)
        fprintf('In (%g,%g,%g), the eigenvalues are: \n',x_val);
        disp(lambda)
    end
end
