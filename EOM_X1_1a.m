function dX = EOM_X1_1a(t,X,p)

% (t+2)/4*100
 
dX = zeros(4,1);

 
T_A = p.T_A0*(1-usgn(p,t))+p.T_Af*usgn(p,t);
T_B = p.T_B0;

% Exponential decay of input torque
% T_A = p.T_A0*(1-usgn(t)+exp(-t/p.beta)*usgn(t));


 
Phi_AB =    p.k_AB*(X(1)-X(2)+p.TH_AB)*(1-u(p,X(1)-X(2)+p.TH_AB))+...
            p.k_AB*(X(1)-X(2)-p.TH_AB)*u(p,X(1)-X(2)-p.TH_AB);       


%Phi_AB = p.k_AB*(X(1)-X(2)) + sign(X(1)-X(2))*p.k_AB*(abs(X(1)-X(2))-p.TH_AB)*sign(abs(X(1)-X(2))-p.TH_AB);
 
dX(1) = X(3);
dX(2) = X(4);
%% Benchmark Case

dX(3) = (T_A-Phi_AB-p.h_A*tanh(50*X(3))-p.h_AB*tanh(50*(X(3)-X(4))))/p.J_A;
dX(4) = (T_B+Phi_AB-p.k_BC*X(2)-p.h_AB*tanh(50*(X(4)-X(3)))-p.h_BC*tanh(50*X(4)))/p.J_B; 

%% Using Smoothening Functions
% dX(3) = (T_A-Phi_AB-p.h_A*ufric(p,X(3))-p.h_AB*ufric(p,X(3)-X(4)))/p.J_A;
% dX(4) = (T_B+Phi_AB-p.k_BC*X(2)-p.h_AB*ufric(p,X(4)-X(3))-p.h_BC*ufric(p,X(4)))/p.J_B; 

%% Impact Damping
% 
% dX(3) = (T_A-Phi_AB*(1+p.idbeta*(X(3)-X(4)))-p.h_A*ufric(p,X(3)))/p.J_A;
% dX(4) = (T_B+Phi_AB*(1+p.idbeta*(X(3)-X(4)))-p.k_BC*X(2)-p.h_BC*ufric(p,X(4)))/p.J_B; 

%% Impact Damping with Hysteresis
% dX(3) = (T_A-Phi_AB*(1+p.beta*(X(3)-X(4)))-p.h_AB*ufric(p,X(3)-X(4))-p.h_A*ufric(p,X(3)))/p.J_A;
% dX(4) = (T_B+Phi_AB*(1+p.beta*(X(3)-X(4)))+p.h_AB*ufric(p,X(3)-X(4))-p.k_BC*X(2)-p.h_BC*ufric(p,X(4)))/p.J_B;

%% Viscous Damping
% dX(3) = (T_A-Phi_AB - p.C_v*(X(3)-X(4))-p.h_AB*tanh(50*(X(3)-X(4)))-p.h_A*ufric(p,X(3)))/p.J_A;
% dX(4) = (T_B+Phi_AB+ p.C_v*(X(3)-X(4))+p.h_AB*tanh(50*(X(3)-X(4)))-p.k_BC*X(2)-p.h_BC*ufric(p,X(4)))/p.J_B; 

%% Viscous Damping without Hysteresis
% dX(3) = (T_A-Phi_AB - p.C_v*(X(3)-X(4))-p.h_A*ufric(p,X(3)))/p.J_A;
% dX(4) = (T_B+Phi_AB+ p.C_v*(X(3)-X(4))-p.k_BC*X(2)-p.h_BC*ufric(p,X(4)))/p.J_B; 

%% Viscous Damping at AB and BC with hysteresis

% dX(3) = (T_A-Phi_AB - p.C_v*(X(3)-X(4))-p.h_AB*tanh(50*(X(3)-X(4)))-p.h_A*ufric(p,X(3)))/p.J_A;
% dX(4) = (T_B+Phi_AB + p.C_v*(X(3)-X(4))+p.h_AB*tanh(50*(X(3)-X(4)))-p.k_BC*X(2)-p.C_v1*(X(4)))/p.J_B; 

%% Impact, Viscous and Hysteresis damping at AB and Viscous+Hysteresis at BC
% dX(3) = (T_A-Phi_AB*(1+p.idbeta*(X(3)-X(4))) - p.C_v*(X(3)-X(4))-p.h_AB*tanh(50*(X(3)-X(4)))-p.h_A*ufric(p,X(3)))/p.J_A;
% dX(4) = (T_B+Phi_AB*(1+p.idbeta*(X(3)-X(4))) + p.C_v*(X(3)-X(4))+p.h_AB*tanh(50*(X(3)-X(4)))-p.k_BC*X(2)-p.C_v1*(X(4)))/p.J_B; 

%% Impact and VIscous without Hysteresis 
% dX(3) = (T_A-Phi_AB*(1+p.beta*(X(3)-X(4))) - p.C_v*(X(3)-X(4))-p.h_A*ufric(p,X(3)))/p.J_A;
% dX(4) = (T_B+Phi_AB*(1+p.beta*(X(3)-X(4))) + p.C_v*(X(3)-X(4))-p.k_BC*X(2)-p.C_v1*(X(4)))/p.J_B; 
end



