% --- theoretical solution generator ---

function [T_out,theta_out,Ux_out,Uy_out] = input_generator_Norm(c,B)
% read parameters from global varibales
global a theta0 Vol M J k_theta k_s k_l m fai theta_st ex_st phase_shift
a_st = a*(1+ex_st);
Theta0 = theta0+theta_st;

Ex = cos(fai)^2/(cos(fai)^2+k_s/k_l*sin(fai)^2-M*c^2/k_l/a_st^2);
Ey = sin(fai)^2/(sin(fai)^2+k_s/k_l*cos(fai)^2-M*c^2/k_l/a_st^2);
Eyy = sin(fai)/(sin(fai)^2+k_s/k_l*cos(fai)^2-M*c^2/k_l/a_st^2);
F = -J*c^2/k_l-sin(Theta0)^2/4/cos(theta0)^2*a^2*a_st^2 ...
    +cos(Theta0)^2/4/cos(theta0)^2*k_s/k_l*a^2*a_st^2 ...
    -a_st^2*k_theta/k_l;

Ex = double(Ex);
Ey = double(Ey);
Eyy = double(Eyy);
F = double(F);
Theta0 = double(Theta0);

C1 = (-a^2*sin(Theta0)^2/cos(theta0)^2*(Ex+Ey)-2*a^2*cos(2*Theta0)/cos(theta0)^2 ...
    +2*a*a_st*cos(Theta0)/cos(theta0)+8*k_theta/k_l ...
    -sqrt(2)*m*B*Vol/2/k_l*(cos(Theta0)+sin(Theta0)))/F;
C2 = (-3*a^2*sin(2*Theta0)/4/cos(theta0)^2*(Ex+Ey)+2*a^2*sin(2*Theta0)/cos(theta0)^2 ...
    -a*a_st*sin(Theta0)/cos(theta0) ...
    -sqrt(2)*m*B*Vol/4/k_l*(cos(Theta0)-sin(Theta0)))/F;
C3 = (a^2*(1-7*cos(2*Theta0))/12/cos(theta0)^2*(Ex+Ey)+4*a^2*cos(2*Theta0)/3/cos(theta0)^2 ...
    -a*a_st*cos(Theta0)/3/cos(theta0) ...
    +sqrt(2)*m*B*Vol/12/k_l*(cos(Theta0)+sin(Theta0)))/F;

C1 = double(C1);
C2 = double(C2);
C3 = double(C3);

D1 = -C2/3/C1;
D2 = sqrt(C2^2/9/C1^2-C3/2/C1);
W = 1/sqrt(C1);
%%%% A+
A = 1/(D1+D2);
% %%%% A-
% A = 1/(D1-D2);

increment = 0.5*a;
X = 0:increment:500*a;
%%%% A+
theta_out = 1./(D1+D2*cosh((X-phase_shift)./W));
% %%%% A-
% theta_out = 1./(D1-D2*cosh((X-phase_shift)./W));
        
disp(['width=',num2str(W/a)]);
disp(['amplitude=',num2str(A)]);
disp(['c=',num2str(c)]); 


H = 6*cos(Theta0)*D1*(-D1^2+D2^2)+sin(Theta0)*(2*D1^2*(1+12*D2^2) ...
    +D2^2*(1-12*D2^2)-12*D1^4);

%%%% A+
Ux_out = a*W*Ex/12/a_st/cos(fai)/cos(theta0)/(-D1^2+D2^2)^(5/2)*(theta_out.*sin(Theta0).*tanh((X-phase_shift)./W)/(D1./cosh((X-phase_shift)./W)+D2)*D2*(-D1^2+D2^2)^(3/2) ...
        +2*H*(atan((D1-D2)/sqrt(-D1^2+D2^2))-atan((D1-D2)/sqrt(-D1^2+D2^2)*tanh((X-phase_shift)/2/W))) ...
        -3*(1-tanh((X-phase_shift)./W)*D2./(D1./cosh((X-phase_shift)./W)+D2))*(-sin(Theta0)*D1+2*cos(Theta0)*(D1^2-D2^2))*sqrt(-D1^2+D2^2));
    
Uy_out = a*W*Eyy/12/a_st/cos(theta0)/(-D1^2+D2^2)^(5/2)*(theta_out.*sin(Theta0).*tanh((X-phase_shift)./W)/(D1./cosh((X-phase_shift)./W)+D2)*D2*(-D1^2+D2^2)^(3/2) ...
        +2*H*(atan((D1-D2)/sqrt(-D1^2+D2^2))-atan((D1-D2)/sqrt(-D1^2+D2^2)*tanh((X-phase_shift)/2/W))) ...
        -3*(1-tanh((X-phase_shift)./W)*D2./(D1./cosh((X-phase_shift)./W)+D2))*(-sin(Theta0)*D1+2*cos(Theta0)*(D1^2-D2^2))*sqrt(-D1^2+D2^2));
    
%     %%%% A-
%     Ux_out = a*W*Ex/12/a_st/cos(fai)/cos(theta0)/(-D1^2+D2^2)^(5/2)*(-theta_out.*sin(Theta0).*tanh((X-phase_shift)./W)/(D1./cosh((X-phase_shift)./W)-D2)*D2*(-D1^2+D2^2)^(3/2) ...
%         +2*H*(atan((D1+D2)/sqrt(-D1^2+D2^2))-atan((D1+D2)/sqrt(-D1^2+D2^2)*tanh((X-phase_shift)/2/W))) ...
%         -3*(1+tanh((X-phase_shift)./W)*D2./(D1./cosh((X-phase_shift)./W)-D2))*(-sin(Theta0)*D1+2*cos(Theta0)*(D1^2-D2^2))*sqrt(-D1^2+D2^2));
%     
%     Uy_out = a*W*Eyy/12/a_st/cos(theta0)/(-D1^2+D2^2)^(5/2)*(-theta_out.*sin(Theta0).*tanh((X-phase_shift)./W)/(D1./cosh((X-phase_shift)./W)-D2)*D2*(-D1^2+D2^2)^(3/2) ...
%         +2*H*(atan((D1+D2)/sqrt(-D1^2+D2^2))-atan((D1+D2)/sqrt(-D1^2+D2^2)*tanh((X-phase_shift)/2/W))) ...
%         -3*(1+tanh((X-phase_shift)./W)*D2./(D1./cosh((X-phase_shift)./W)-D2))*(-sin(Theta0)*D1+2*cos(Theta0)*(D1^2-D2^2))*sqrt(-D1^2+D2^2));

    Ux_out = double(-Ux_out/a);
    Uy_out = double(-Uy_out/a);
    Ux_out = Ux_out-Ux_out(1);
    Uy_out = Uy_out-Uy_out(1);

% for ii = 1:length(X)
%     Ux_out(ii) = a*W*Ex/12/a_st/cos(fai)/cos(theta0)/(-D1^2+D2^2)^(5/2)*( ...
%         +theta_out(ii)^2*sin(Theta0)*sinh(X(ii)/W)*D2*(-D1^2+D2^2)^(3/2) ...
%         +2*H*(atan((D1-D2)/sqrt(-D1^2+D2^2))-atan((D1-D2)/sqrt(-D1^2+D2^2)*tanh(X(ii)/2/W))) ...
%         -3*(1-theta_out(ii)*D2*sinh(X(ii)/W))*(-sin(Theta0)*D1+2*cos(Theta0)*(D1^2-D2^2))*sqrt(-D1^2+D2^2));
% end
T_out = X/c*sqrt(k_l/M);
% taking the real part of solutions
% theta_out = real(theta_out);
% U_out = real(U_out);
end
