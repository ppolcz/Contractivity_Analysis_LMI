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

[o,r] = TRM_definitions(4,4,'Visualize',2,"Layout","force","View",[35,90], ...
    "SubDIR","TRM_4D_v4");

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
% Computer algebra manipulations to find out, which elements from P should be eliminated
% to obtain a diagonal metric.

P_sym = pcz_sym_symmetric('a',m);
PI_sym = sym(PI);

M_sym = PI_sym.' * P_sym * PI_sym;

[i,j] = ndgrid(1:n);
abdiagM = M_sym(i~=j);

cc = sym([]);
for k = 1:numel(abdiagM)
    [c,t] = coeffs(abdiagM(k),x);
    cc = [cc c];
end

Sparsity_Pattern = setdiff(symvar(P_sym),symvar(cc));
P_sym0 = subs(P_sym,Sparsity_Pattern,Sparsity_Pattern*0);
Idx0 = find(P_sym0); % <---- These elements should be zero

%%

% These elements can also be eliminated to obtain a sparser solution:
Idx1 = [8, 10, 13, 19, 25, 34, 37, 42, 49, 54, 58, 67, 73, 79, 84, 88, 96, 99, 108, 111, 115, 118, 126, 127, 134, 139, 142, 145, 151, 158, 161, 166, 169, 176, 178, 185, 188, 191, 196];

% Sparsity pattern of P is given manually.
% Idx0 = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196];

Eps = 0.0;
Clip = [Eps 1-Eps];
x_lim = ones(n,1) * Clip;

X_v = P_ndnorms_of_X(x_lim);

beta = sdpvar;
L = sdpvar(m,s,'full');
Ld = sdpvar(md,sd,'full');
P = sdpvar(m);
P(Idx0) = 0;

% ---------------------------------------------------------------------
% Further elements are eliminated manually. (You can comment this out.)
P(Idx1) = 0; % -------------------------------------------------------- 

CONS = [];

epsilon1 = 0.001;
betaM = 0.0001;
for i = 1:size(X_v,1)
    xi = X_v(i,:)';
    
    CONS = [ CONS
        P + He( L*N(xi) ) - epsilon1*(Ia'*Ia) >= 0
        He( Ed'*P*Ad + Ld*Nd(xi) ) + betaM*Ed'*P*Ed <= 0
        ];
end

opts = sdpsettings('solver','sedumi');
sol = optimize(CONS,[],opts);
pcz_feasible(sol,CONS)

%{
SeDuMi 1.3 by AdvOL, 2005-2008 and Jos F. Sturm, 1998-2003.
Alg = 2: xz-corrector, theta = 0.250, beta = 0.500
eqs m = 6363, order n = 1025, dim = 43137, blocks = 33
nnz(A) = 4880220 + 0, nnz(ADA) = 37413369, nnz(L) = 18709866
 it :     b*y       gap    delta  rate   t/tP*  t/tD*   feas cg cg  prec
  0 :            8.25E+01 0.000
  1 :   0.00E+00 2.88E+01 0.000 0.3492 0.9000 0.9000   0.95  1  1  1.0E+02
  2 :   0.00E+00 9.78E+00 0.000 0.3395 0.9000 0.9000   0.95  1  1  3.6E+01
  3 :   0.00E+00 3.69E+00 0.000 0.3768 0.9000 0.9000   1.00  1  1  1.3E+01
  4 :   0.00E+00 7.74E-01 0.000 0.2099 0.9000 0.9000   1.01  1  1  2.8E+00
  5 :   0.00E+00 1.55E-02 0.000 0.0200 0.9900 0.9900   1.00  1  1  5.6E-02
  6 :   0.00E+00 1.13E-07 0.000 0.0000 1.0000 1.0000   1.00  1  1  4.2E-07
Run into numerical problems.

iter seconds digits       c*x               b*y
  6    662.4   1.3 -9.8883720065e-11  0.0000000000e+00
|Ax-b| =   5.5e-07, [Ay-c]_+ =   0.0E+00, |x|=  1.7e-06, |y|=  2.6e+04

Detailed timing (sec)
   Pre          IPM          Post
3.181E+00    4.928E+02    6.676E-02    
Max-norms: ||b||=0, ||c|| = 1.000000e-03,
Cholesky |add|=0, |skip| = 890, ||L.L|| = 153.72.

sol = 
    yalmipversion: '20230622'
    matlabversion: '24.2.0.2773142 (R2024b) Update 2'
       yalmiptime: 1.4879
       solvertime: 496.2733
             info: 'Successfully solved (SeDuMi)'
          problem: 0

[   OK   ] Successfully solved (SeDuMi). Solver time: 496.273
[   OK   ] The solution is feasible. Tolerance: 1e-10.
%}

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

%%

m_fh = matlabFunction(diag(M).','Vars',{o.x_sym.'});
f_ode = matlabFunction(o.f_sym,'Vars',{'t',o.x_sym});

x0s = [
    0.7084 0.9667 , 0.9621 0.8078 , 0.9680 0.6470 , 0.9610 0.7809 , 0.8976 0.9894 , 0.7508 0.8989 , 0.8408 0.9879
    0.8673 0.9766 , 0.6357 0.5580 , 0.9912 0.8220 , 0.8829 0.8459 , 0.9626 0.9898 , 0.5030 0.5625 , 0.7565 0.8259
    0.0267 0.1284 , 0.1109 0.0242 , 0.1911 0.0022 , 0.0940 0.0371 , 0.0108 0.0361 , 0.0053 0.0823 , 0.0929 0.1825
    0.2561 0.8484 , 0.7439 0.3593 , 0.7379 0.1603 , 0.3903 0.1410 , 0.2559 0.4442 , 0.1015 0.4649 , 0.0046 0.3553
    ];

x0s = x0s(:,9+[0,1]);

x0s = round(x0s,2,"significant");

fprintf('%g, ',x0s(:,1))
disp ' '
fprintf('%g, ',x0s(:,2))
disp ' '

N = size(x0s,2);

T = 10;
dt = 0.0001;

tt = linspace(0,T,T/dt);

xx = zeros(numel(tt),4,size(x0s,2));

for i = 1:N
    [~,xx(:,:,i)] = ode45(f_ode,tt,x0s(:,i));
end

fig = figure(12);
fig.Position(3:4) = [549 293];
delete(fig.Children);
ax = axes(fig);
hold on;

dtt = tt(1:end-1);
for i = 1:N-1
    for j = i+1:N
        deltax = xx(:,:,i) - xx(:,:,j);
        ndeltax = diff(vecnorm((deltax)')');
        if any(ndeltax > 0)

            mdeltax = diff(sqrt(sum( m_fh(xx(:,:,i)) .* deltax .* deltax , 2 )));
            mdeltax2 = diff(sqrt(sum( m_fh(xx(:,:,j)) .* deltax .* deltax , 2 )));

            Pl1 = plot(dtt,ndeltax,'k','LineWidth',2);
            Pl2 = plot(dtt,mdeltax,'--','Color',Pl1.Color,'LineWidth',2);
            Pl3 = plot(dtt,mdeltax2,':','Color',Pl1.Color,'LineWidth',2);
            % plot(dtt,mdeltax,':');

            [ xx(1,:,i)' , xx(1,:,j)' ]
        end
    end
end
grid on, box on,
yline(0,'-','LineWidth',1,'Color',[0.4 0.4 0.4])
ax.XAxis.Scale = 'log';
ax.XLim = [1e-4 10];

xlabel('time in logarithmic scale','Interpreter','latex');
ylabel({'~~~~~~~~~~time derivative of the norms' '~~~~~~~~of the virtual displacement'},'Interpreter','latex');

Leg = legend('Rate of $\|{\delta x(t)}\|$, where $\delta x(t) = x(t) - x''(t)$', ...
    'Rate of $\big(\delta x(t)\big)^\top M\big(x(t)\big) \, \delta x(t)$', ...
    'Rate of $\big(\delta x(t)\big)^\top M\big(x''(t)\big) \, \delta x(t)$');
Leg.Interpreter = 'latex';
Leg.FontSize = 13;
Leg.Location = "northoutside";
Leg.Box = 'off';

ax.FontSize = 14;
ax.TickLabelInterpreter = 'latex';

exportgraphics(fig,'/home/ppolcz/Dropbox/Peti/PhD/.Dokumentaciok/37_TRM_Contraction/NOLCOS/fig/contractive.pdf','ContentType','vector')
