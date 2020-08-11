function dX = EOM_X1_1actual(t,X,p)
 
% (t+2)/4*100
 
dX = zeros(4,1);
 
 
T_A = p.T_A0*(1-u(t))+p.T_Af*u(t);
T_B = p.T_B0;
 
Phi_AB =    p.k_AB*(X(1)-X(2)+p.TH_AB)*(1-u(X(1)-X(2)+p.TH_AB))+...
            p.k_AB*(X(1)-X(2)-p.TH_AB)*u(X(1)-X(2)-p.TH_AB);
 
dX(1) = X(3);
dX(2) = X(4);
dX(3) = (T_A-Phi_AB-p.h_A*tanh(50*X(3))-p.h_AB*tanh(50*(X(3)-X(4))))/p.J_A;
dX(4) = (T_B+Phi_AB-p.k_BC*X(2)-p.h_AB*tanh(50*(X(4)-X(3)))-p.h_BC*tanh(50*X(4)))/p.J_B;    
 
end
