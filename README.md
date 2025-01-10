# Contractivity analysis of nonlinear rational autonomous models

This code uses the LFR Toolbox [1] and the GSS Library [2] of the SMAC Toolbox [7]. For modeling semidefinite program, we use YALMIP [4] with SeDuMi [3] or Mosek [6] solvers. To formulate sufficient polytopic linear matrix inequalities to solve the contraction problem, we used the numerical methods presented in [5].

[1] J-F. Magni, "Linear Fractional Representation toolbox for use with Matlab", February 2006, available with the SMAC Toolbox at http://w3.onera.fr/smac/lfrt.

[2] J-M. Biannic and C. Roos, "Generalized State Space: a new Matlab class to model uncertain and nonlinear systems as Linear Fractional Representations", 2012-2020, available with the SMAC Toolbox at http://w3.onera.fr/smac/gss.

[3] Sturm, J. F. (1999) ‘Using SeDuMi 1.02, A Matlab toolbox for optimization over symmetric cones’, Optimization Methods and Software, 11(1–4), pp. 625–653. doi: 10.1080/10556789908805766.

[4] J. Lofberg, "YALMIP : a toolbox for modeling and optimization in MATLAB," 2004 IEEE International Conference on Robotics and Automation (IEEE Cat. No.04CH37508), Taipei, Taiwan, 2004, pp. 284-289, doi: 10.1109/CACSD.2004.1393890

[5] Polcz, P., Péni, T., Kulcsár, B., and Szederkényi, G.
(2020). Induced L2-gain computation for rational LPV
systems using Finsler’s lemma and minimal generators.
Systems & Control Letters, 142, 104738. doi:10.1016/j.
sysconle.2020.104738.

[6] MOSEK ApS (2015). The MOSEK optimization toolbox
for MATLAB manual. version 7.1 (revision 28). URL
http://docs.mosek.com.

[7] Biannic, J.M., Burlion, L., Demourant, F., Ferreres, G.,
Hardier, G., Loquen, T., and Roos, C. (2016). The
SMAC toolbox: A collection of libraries for systems
modeling, analysis and control. online available at
http://w3.onera.fr/smac/