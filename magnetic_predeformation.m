% --- Magnetic field-induced predeformation ---

function [theta_st,ex_st,ey_st] = magnetic_predeformation(B)

% read parameters from global varibales
global k_theta k_l m a theta0 Vol

syms theta_pre ex_pre ey_pre

equ(1) = -8*k_theta*theta_pre-k_l*a^2/cos(theta0)*sin(theta_pre+theta0)*(2+ex_pre+ey_pre)...
    +2*k_l*a^2/cos(theta0)/cos(theta0)*sin(theta_pre+theta0)*cos(theta_pre+theta0)...
    -m*B*Vol*cos(pi/4+theta_pre+theta0); % rotate motion

equ(2) = k_l*a*(1+ex_pre-cos(theta_pre+theta0)/cos(theta0)); % motion along x direction

equ(3) = k_l*a*(1+ey_pre-cos(theta_pre+theta0)/cos(theta0)); % motion along y direction

S = vpasolve([equ(1)==0,equ(2)==0,equ(3)==0],[theta_pre,ex_pre,ey_pre]);

ex_st=S.ex_pre;
ey_st=S.ey_pre;
theta_st=S.theta_pre;

end