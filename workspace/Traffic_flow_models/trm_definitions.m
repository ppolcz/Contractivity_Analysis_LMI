%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2024. November 22. (2023a)
function [s,r] = trm_definitions(dim,variant,s)
arguments
    dim
    variant = 1
    s.Check = 0
    s.Eq = 0
    s.Visualize = 0
    s.Layout = "auto"
    s.View = [0 90];
    s.DIR = "/home/ppolcz/Dropbox/Peti/PhD/.Dokumentaciok/37_TRM_Contraction/results";
    s.SubDIR = string(datetime('today','Format','yyyy-MM-dd'))
end
%% Select [K k_on k_off C], which unique determines the traffic flow reaction model
% Every other objects and variables are generated automatically based on these matrices.

switch dim

    case 3
        K_ = dim3(variant);

    case 4
        K_ = dim4(variant);

    case 5
        K_ = dim5(variant);

    case 6
        K_ = dim6(variant);

end

% Store this matrix in r (i.e., struct for the results, which will be backed up later).
r.K_on_off_C = K_;

%% Construct the symbolic model

mkdir(fullfile(s.DIR,s.SubDIR));

% Detect the dimension of the model
n = height(K_);

s.ModelName = sprintf('trm_dim%d_v%d',n,variant);
s.dim = dim;
s.variant = variant;

% Generate SMT symbolic objects
x = sym('x',[n,1]);
assumeAlso(in(x,'real'));

% Split up matrix K_, whic contains the transition graph's matrix (K), k_on, k_off, and C
K = K_(:,1:n);
k_on = K_(:,n+1);
k_off = K_(:,n+2);
C = K_(:,n+3);

Transitions = diag(x) * K * diag(C-x);

f = sum(Transitions,1).'-sum(Transitions,2) - k_off.*x + k_on.*(C-x);
A = jacobian(f,x);

s.x_sym = x;
s.f_sym = f;
s.A_sym = A;

%% Generate LFR Toolbox objects

Eps = 0.05;
technical_limits_for_x = C .* [Eps 1-Eps];
[s.x_lfr,s.xi_lfr] = pcz_generateLFRStateVector('x',technical_limits_for_x);

% Convert f(x) to LFR
f_fh_cell = cellfun(@(f) {matlabFunction(f,'vars',x)}, num2cell(f));
s.fi_lfr = cellfun(@(f) {f(s.xi_lfr{:})},f_fh_cell);
s.f_plfr = plfr(vertcat(s.fi_lfr{:}));

% Convert A(x) to LFR
A_fh_cell = cellfun(@(f) {matlabFunction(f,'vars',x)}, num2cell(A));
s.Aij_lfr = cellfun(@(f) {f(s.xi_lfr{:})},A_fh_cell);

s.Ai_lfr = cellfun(@(c) {[c{:}]}, num2cell(s.Aij_lfr,2));
s.A_plfr = plfr(vertcat(s.Ai_lfr{:}));

if s.Check
    pcz_symzero_report(sym(s.A_plfr) - A,'A_lfr = A_sym')
end

%% Compute the equilibrium point(s)

if s.Eq
    
    import casadi.*
    
    Eps = 0.01;
    x_lim = C .* [Eps 1-Eps];
    
    helper = Pcz_CasADi_Helper('SX');
    
    x = helper.new_var('x',[n,1],'lb',C*0,'ub',C);
    
    Transitions = diag(x) * K * diag(C-x);
    
    f = sum(Transitions,1)'-sum(Transitions,2) - k_off.*x + k_on.*(C-x);
    helper.add_eq_con(f);
    
    Nr_Samples = 30;
    IC = rand(n,Nr_Samples).*diff(x_lim,[],2) + x_lim(:,1);
    solver = helper.get_nl_solver;
    sol = solver.solve([],IC);
    
    Eq = full(sol.x');
    
    [~,Idx] = unique(round(Eq,3),'rows');
    
    % return values:
    % ----------------
    s.Eq = Eq(Idx,:)';
    % ----------------
end

%% Visualize

if s.Visualize
    On = zeros(n);
    In = eye(n);
    
    ion = find(k_on > 0);
    ioff = find(k_off > 0);
    
    K__ = [
        K On(:,ion) In(:,ioff)
        In(ion,:) zeros(numel(ion),numel(ioff)+numel(ion))
        zeros(numel(ioff),n+numel(ioff)+numel(ion))
        ];
    
    Names = cellfun(@(i) {num2str(i)},num2cell(1:n));
    On = cellfun(@(i){['on' num2str(i)]},num2cell(1:numel(ion)));
    Off = cellfun(@(i){['off' num2str(i)]},num2cell(1:numel(ioff)));
    Names = [Names, On, Off];
    
    G = digraph(K__,Names);
    
    fig = figure(673);
    Tl = tiledlayout(1,1,'Padding','none');
    ax = nexttile;
    plot(G,"ArrowSize",15,"MarkerSize",8,"NodeFontSize",15,"Layout",s.Layout)
    ax.Visible = 'off';
    view(s.View);
   
    if s.Visualize > 1
        fname = fullfile(s.DIR,s.SubDIR,s.ModelName + ".pdf");
        exportgraphics(gcf,fname,'ContentType','vector');
    end
end

end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function K_on_off_C = dim4(variant)
%%
    switch variant
        case 1
            K_on_off_C = [
                0   2   0   0 , 1 , 0 , 1
                0   0   3   0 , 0 , 0 , 1
                0   0   0   2 , 0 , 1 , 1
                1.5 0   0   0 , 0 , 0 , 1
                ];

        case 2
            K_on_off_C = [
                0   2   0   0 , 0 , 0 , 1
                0   0   3   0 , 0 , 1 , 1
                0   0   0   2 , 0 , 0 , 1
                1.5 0   0   0 , 1 , 0 , 1
                ];

        case 3
            K_on_off_C = [
                0   2   0   0 , 0 , 0 , 1
                0   0   3   0 , 0 , 1 , 1
                0   0   0   2 , 1 , 0 , 1
                1.5 0   0   0 , 0 , 0 , 1
                ];

        case 4
            K_on_off_C = [
                0   2   0.5 0 , 0 , 0 , 1
                0   0   3   0 , 0 , 1 , 1
                0   0   0   2 , 1 , 0 , 1
                1.5 0   0   0 , 0 , 0 , 1
                ];
    end

end

function K_on_off_C = dim5(variant)
%%
    
    switch variant
        case 1
        
            K_on_off_C = [
                0   1   0   0   0   , 1 , 0 , 1
                0   0   3   0   1   , 0 , 0 , 1
                0   0   0   2   0   , 0 , 1 , 1
                0   0   0   0   1.5 , 0 , 0 , 1
                1.5 0   0   0   0   , 0 , 0 , 1
                ];
    
    end
end
function K_on_off_C = dim6(variant)
%%
    
    switch variant
        case 1
        
            K_on_off_C = [
                0   2   0   0   0   0 , 1 , 0 , 1
                0   0   1   0   1   0 , 0 , 0 , 1
                0   0   0   1   0   0 , 0 , 1 , 1
                0   0   0   0   1   0 , 0 , 0 , 1
                0   0   0   0   0   1 , 0 , 0 , 1
                1   0   0   0   0   0 , 0 , 1 , 1
                ];
    
    end
end

