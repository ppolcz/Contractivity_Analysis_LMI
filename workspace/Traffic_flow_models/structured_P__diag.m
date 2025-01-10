function [P_sdp,P_str,M_str] = structured_P__diag(x,PI)
%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2024. November 22. (2023a)
%


PI_sym = sym(PI);
[m,n] = size(PI_sym);

Mtp = subs(PI_sym,x,ones(n,1));
Mtp(Mtp == 0) = 1;
PI_nice = PI_sym ./ Mtp;

[i,j] = meshgrid(1:m);
P = sym('p',m);
P(i<j) = 0;
P = P + P.' - diag(diag(P));

M = collect(PI_nice.' * P * PI_nice,x);

% diagM = diag(M);

[i,j] = meshgrid(1:n);
abdiagM = M(i~=j);

cc = {};
tt = {};
for k = 1:numel(abdiagM)
    [c,t] = coeffs(abdiagM(k),x);
    cc = [cc {c}];
    tt = [tt {t}];
end
Eq = [cc{:}].';

ZeroVars = symvar(Eq);
P_str = subs(P,ZeroVars,ZeroVars*0);
M_str = collect(PI_nice.' * P_str * PI_nice,x);

Vars = symvar(P_str);
P_fh = matlabFunction(P_str,"Vars",{Vars});
P_sdp = P_fh(sdpvar(1,numel(Vars)));

end