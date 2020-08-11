format compact; clear all; close all; clc;
 
p.J_A = 0.005080424;    %torsional inertia [kgm^2] (measured from step-response)
p.J_B = 0.000397876;    %torsional inertia [kgm^2] (measured from impulse-response)
p.J_C = p.J_B;
 
p.k_AB = 8600;          %torsional stiffness [Nm/rad] (static measurement)
p.k_BC = 47.99;          %torsional stiffness [Nm/rad] (static measurement)
p.k_C = p.k_AB;
 
p.TH_AB = 1.5*0.1*pi/180;  %clearance transition [rad] (measured from geometry)
p.TH_C = 0.1*pi/180;
 
p.h_A = 0.033078;       %Coulomb friction [Nm] (measured from step-response)
p.h_AB = (0.2803-p.h_A)/1.5;  %Coulomb friction [Nm] (0.375-0.625 kg at 0.0762 m)
p.h_BC = 0.028075;       %Coulomb friction [Nm] (measrued from impulse-response)
p.h_C = p.h_AB;
 
p.T_A0 = 2.877*9.81*5.125*2.54/100;  %initial external torque [Nm]
p.T_Af = 0;                          %final external torque [Nm]
 
p.T_B0 = 0;
 
tau_s = 1/12.8e3;                   %measurement period
 
J = diag([p.J_A p.J_B p.J_C]);
K = [   p.k_AB      -p.k_AB         0;
        -p.k_AB     p.k_AB+p.k_BC   -p.k_BC;
        0           -p.k_BC         p.k_BC+p.k_C;   ];
lambda = eig(K,J);
w1 = lambda(1).^0.5;
t_S = 2*pi/w1;
th_S = 4.5*pi/180;
dth_S = th_S/t_S;
d2th_S = dth_S/t_S;
 
p.t_S = t_S;
p.th_S = th_S;
p.dth_S = dth_S;
p.d2th_S = d2th_S;
 
% w_n = sqrt(p.k_BC/(p.J_A+p.J_B));
% tau_n = 2*pi/w_n;
% t_S = tau_n;
% th_S = p.TH_AB;
% dth_S = th_S/t_S;
% d2th_S = dth_S/t_S;
 
th_B0 = p.T_A0/p.k_BC;
th_A0 = th_B0+p.TH_AB+p.T_A0/p.k_AB;
 
options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Maxstep',tau_s);
t_span = -2:tau_s:2;
X_0 = [th_A0 th_B0 0 0];    %initial conditions
% X_0 = [0 0 0 0];    %initial conditions
a = cputime;
[t,X] = ode45(@EOM_X1_1actual,t_span,X_0,options,p);
cputime-a
 
%%
LW = 2;
 
%Results
th_A = X(:,1);          %angular displacement [rad]
th_B = X(:,2);          %angular displacement [rad]+
th_AB = th_A-th_B;      %relative angular displacement [rad]
dth_A = X(:,3);         %angular velcoity [rad/s]
dth_B = X(:,4);         %angular velocity [rad/s]
dth_AB = dth_A-dth_B;   %relative angular velocity [rad/s]
 
T_A = p.T_A0*(1-u(t))+p.T_Af*u(t);    %external torque [Nm]
T_B = p.T_B0;
% for i=1:length(t)
% if t(i)<0;
%     T_A(i) =(p.T_A0/2);
%     T_B(i) =0;
% elseif 0<=t(i) && p.t_S>=t(i)
% T_A(i) = (p.T_A0/2)+(1.75*(cos((2*pi*t(i)/p.t_S)+(pi/2))));
% T_B(i) = 0;
% elseif t(i)>p.t_S 
%     T_A(i) = (p.T_A0/2);
%     T_B(i) = 0;
% end
% end
 
Phi_AB =    p.k_AB.*(th_AB+p.TH_AB).*(1-u(th_AB+p.TH_AB))+...
            p.k_AB.*(th_AB-p.TH_AB).*u(th_AB-p.TH_AB);          %elastic torque through clearance [Nm]
 
d2th_A = (T_A-Phi_AB-p.h_A*tanh(50*dth_A)-p.h_AB*tanh(50*dth_AB))/p.J_A;             %angular acceleration [rad/s^2]
d2th_B = (T_B+Phi_AB-p.k_BC*th_B+p.h_AB*tanh(50*dth_AB)-p.h_BC*tanh(50*dth_B))/p.J_B;  %angular acceleration [rad/s^s]
%% Peaks
 
[A_pks_hi,A_loc_hi] = findpeaks(d2th_A./d2th_S,t./t_S,'MinPeakHeight',10,'MinPeakDistance',0.42);
[A_pks_low,A_loc_low] = findpeaks(-d2th_A./d2th_S,t./t_S,'MinPeakHeight',10,'MinPeakDistance',0.42);
[B_pks_hi,B_loc_hi] = findpeaks(d2th_B./d2th_S,t./t_S,'MinPeakHeight',50,'MinPeakDistance',0.30);
[B_pks_low,B_loc_low] = findpeaks(-d2th_B./d2th_S,t./t_S,'MinPeakHeight',70,'MinPeakDistance',0.30);
 
figure(1)
plot(t./p.t_S,T_A,'g-','LineWidth',2);
title('Input: Step Down Torque');
xlabel('$\overline{t}[s]$','interpreter','latex','FontSize',15);
ylabel('T [Nm]');
 
figure(2)
plot(t./p.t_S,d2th_A./p.d2th_S,'g-','LineWidth',2);
hold on
plot(A_loc_hi,A_pks_hi,'gO');
plot(A_loc_low,-A_pks_low,'rO');
title('Acceleration at inertia A');
xlabel('$\overline{t}[s]$','interpreter','latex','FontSize',15);
ylabel('$\overline{\ddot{\theta}}_{A}$','interpreter','latex','FontSize',15);
% xlim([-0.5 10]);
 
figure(3)
plot(t./p.t_S,d2th_B./p.d2th_S,'g-','LineWidth',2);
hold on
plot(B_loc_hi,B_pks_hi,'gO');
plot(B_loc_low,-B_pks_low,'rO');
title('Acceleration at inertia B');
xlabel('$\overline{t}[s]$','interpreter','latex','FontSize',15);
ylabel('$\overline{\ddot{\theta}}_{B}$','interpreter','latex','FontSize',15);
% xlim([-0.5 10]);
 
%% Time Periods
 
% ind = find(A_loc_hi > 0)
% A_loc_hi = A_loc_hi(ind(1):end);
% A_loc_low = A_loc_low(ind(1):end);
% 
% for i = 1:(length(A_loc_hi(:,1))-1)
%                 
%                 pos_tau_A(i) = A_loc_hi(i+1)-A_loc_hi(i);
%                 pos_Tosc_A(i) = (A_loc_hi(i+1)+A_loc_hi(i))/2;
%                 
% end
% 
% for i = 1:(length(A_loc_low(:,1))-1)
%                 
%                 neg_tau_A(i) = A_loc_low(i+1)-A_loc_low(i);
%                 neg_Tosc_A(i) = (A_loc_low(i+1)+A_loc_low(i))/2;
%                 
% end
% 
% A_tau = [neg_tau_A pos_tau_A];
% A_tosc = [neg_Tosc_A pos_Tosc_A];
% 
% ind = find(B_loc_hi > 0)
% B_loc_hi = B_loc_hi(ind(1):end);
% B_loc_low = B_loc_low(ind(1):end);
% 
% for i = 0:(length(B_loc_hi(:,1))-4)/2
%                 
%                 pos_tau_B(i+1) = B_loc_hi(2*i+4)-B_loc_hi(2*i+2);
%                 pos_Tosc_B(i+1) = (B_loc_hi(2*i+4)+B_loc_hi(2*i+2))/2;
%                 
% end
% 
% for i = 0:(length(B_loc_low(:,1))-3)/2
%                 
%                 neg_tau_B(i+1) = B_loc_low(2*i+3)-B_loc_low(2*i+1);
%                 neg_Tosc_B(i+1) = (B_loc_low(2*i+3)+B_loc_low(2*i+1))/2;
%                 
% end
% 
% B_tau = [neg_tau_B pos_tau_B];
% B_tosc = [neg_Tosc_B pos_Tosc_B];
% 
% figure(4)
% plot(A_tosc,A_tau,'O','MarkerFaceColor','r');
% hold on
% plot(B_tosc,B_tau,'O','MarkerFaceColor','g');
% % xlim([0 7]);
% ylim([0.95 1.2]);
% xlabel('$\overline{t}$','interpreter','latex','FontSize',15);
% ylabel('$\overline{\tau}$','interpreter','latex','FontSize',15);
% 
x_lim = [-0.5 20];
% LW = 2;
 
figure(5)
hold on
plot(t./t_S,th_AB/th_S,'-g','LineWidth',LW)
ylabel('$\overline{{\theta}}_{AB}$', 'interpreter','latex','FontSize',15);
xlabel('$\overline{t}$','interpreter','latex','FontSize',15);
plot(x_lim,p.TH_AB/th_S*[1 1],'--k')
plot(x_lim,-p.TH_AB/th_S*[1 1],'--k')
xlim([-0.5 20])
box on
 
%%
 
E_h_A = trapz(t,p.h_A*tanh(50*dth_A).*dth_A)
E_h_AB = trapz(t,p.h_AB*tanh(50*dth_AB).*dth_AB)
E_h_BC = trapz(t,p.h_BC*tanh(50*dth_B).*dth_B)
 
 
 
% [expt.B_pks_hi,expt.B_loc_hi] = findpeaks(expt.d2th_Bm./d2th_S,expt.tm./t_S,'MinPeakHeight',100,'MinPeakDistance',0.35);
% [expt.B_pks_low,expt.B_loc_low] = findpeaks(-expt.d2th_Bm./d2th_S,expt.tm./t_S,'MinPeakHeight',50,'MinPeakDistance',0.35);
 
for i = 1:length(dth_A)
   
    KE(i) = abs((0.5*p.J_A*dth_A(i)*dth_A(i)*sign(dth_A(i)))-(0.5*p.J_B*dth_B(i)*dth_B(i)*sign(dth_B(i))));
    
end
 
x_lim = [-0.5 20];
LW = 2;
figure(1)
subplot(2,1,1)
hold on
plot(t/t_S,d2th_A/d2th_S,'-g','LineWidth',LW)
hold on
plot(A_loc_hi,A_pks_hi,'go');
plot(A_loc_low,-A_pks_low,'ro');
xlim(x_lim)
box on
subplot(2,1,2)
hold on
plot(t/t_S,d2th_B/d2th_S,'-g','LineWidth',LW)
hold on
plot(B_loc_hi,B_pks_hi,'go');
plot(B_loc_low,-B_pks_low,'ro');
xlim(x_lim)
box on
 
figure(2)
hold on
plot(t/t_S,th_AB/th_S,'-g','LineWidth',LW)
ylabel('$\overline{{\theta}}_{AB}$', 'interpreter','latex','FontSize',15);
plot(x_lim,p.TH_AB/th_S*[1 1],'--k')
plot(x_lim,-p.TH_AB/th_S*[1 1],'--k')
xlim(x_lim)
xlabel('$\overline{t}$','interpreter','latex','FontSize',15);
 
figure(3)
subplot(211)
plot(t/t_S,dth_A/dth_S,'-g','LineWidth',LW)
ylabel('$\overline{\dot{\theta}}_{A}$', 'interpreter','latex','FontSize',15);
xlim([-1 11])
subplot(212)
plot(t/t_S,dth_B/dth_S,'-g','LineWidth',LW)
ylabel('$\overline{\dot{\theta}}_{B}$', 'interpreter','latex','FontSize',15);
xlim([-1 11])
 
ind = find(A_loc_hi > 0);
A_loc_hi = A_loc_hi(ind(1):end);
A_loc_low = A_loc_low(ind(1):end);
 
for i = 1:(length(A_loc_hi(:,1))-1)
                
                pos_tau_A(i) = A_loc_hi(i+1)-A_loc_hi(i);
                pos_Tosc_A(i) = (A_loc_hi(i+1)+A_loc_hi(i))/2;
                
end
 
for i = 1:(length(A_loc_low(:,1))-1)
                
                neg_tau_A(i) = A_loc_low(i+1)-A_loc_low(i);
                neg_Tosc_A(i) = (A_loc_low(i+1)+A_loc_low(i))/2;
                
end
 
ind = find(B_loc_hi > 0)
B_loc_hi = B_loc_hi(ind(1):end);
B_loc_low = B_loc_low(ind(1):end);
 
for i = 0:(length(B_loc_hi(:,1))-4)/2
                
                pos_tau_B(i+1) = B_loc_hi(2*i+4)-B_loc_hi(2*i+2);
                pos_Tosc_B(i+1) = (B_loc_hi(2*i+4)+B_loc_hi(2*i+2))/2;
                
end
 
for i = 0:(length(B_loc_low(:,1))-3)/2
                
                neg_tau_B(i+1) = B_loc_low(2*i+3)-B_loc_low(2*i+1);
                neg_Tosc_B(i+1) = (B_loc_low(2*i+3)+B_loc_low(2*i+1))/2;
                
end
 
% B_tau = [1.0222,1.026,1.033,1.039,1.057,1.097,1.0208,1.024,1.029,1.036,1.047,1.071,1.161];
% B_Tosc = [1.85,2.876,3.909,4.948,6.005,7.102,1.338,2.362,3.391,4.427,5.474,6.545,7.706];
 
 
 
 
% figure(3)
% plot(neg_Tosc_B,neg_tau_B,'kX');
% hold on
% plot(pos_Tosc_B,pos_tau_B,'kO');
% legend('Predicted : Odd','Predicted : Even');
% title('Time Periods at B');
% 
% B_tau = [neg_tau_B pos_tau_B];
% B_tosc = [neg_Tosc_B pos_Tosc_B];
 
% figure(4)
% plot(neg_Tosc_A,neg_tau_A,'kX');
% hold on
% plot(pos_Tosc_A,pos_tau_A,'kO');
% legend('Predicted : Odd','Predicted : Even');
% title('Time Periods at A');
 
figure(101)
subplot(211)
hold on
plot(t/t_S,-d2th_A/d2th_S,'-g','LineWidth',LW)
ylabel('$\overline{\ddot{\theta}}_{A}$', 'interpreter','latex','FontSize',15);
legend('Experiment','Prediction');
xlim([-1 10])
ylim([-500 500]);
set(gca,'YTick',[-500 0 500]);
subplot(212)
hold on
plot(t/t_S,d2th_B/d2th_S,'-g','LineWidth',LW)
ylabel('$\overline{\ddot{\theta}}_{B}$', 'interpreter','latex','FontSize',15);
xlabel('$\overline{t}$','interpreter','latex','FontSize',15);
xlim([-1 10])
ylim([-2200 2200]);
set(gca,'YTick',[-2000 0 2000]);
 
 
figure(102)
subplot(311)
plot(t./t_S,th_B./th_S,'g','LineWidth',LW);
hold on
ylabel('$\overline{{\theta}}_{B}$', 'interpreter','latex','FontSize',15);
xlim([-1 10]);
subplot(312)
plot(t./t_S,dth_B./dth_S,'g','LineWidth',LW);
xlim([-1 10]);
ylabel('$\overline{\dot{\theta}}_{B}$', 'interpreter','latex','FontSize',15);
subplot(313)
plot(t./t_S,d2th_B./d2th_S,'g','LineWidth',LW);
xlim([-1 10]);
ylabel('$\overline{\ddot{\theta}}_{B}$', 'interpreter','latex','FontSize',15);
 
figure(103)
subplot(211)
plot(t./t_S,d2th_B./d2th_S,'g','LineWidth',LW);
xlim([-1 10]);
ylabel('$\overline{\ddot{\theta}}_{B}$', 'interpreter','latex','FontSize',15);
subplot(212)
plot(t./t_S,KE,'g','LineWidth',LW);
xlim([-1 10]);
ylabel('|\Delta KE|_{AB}');
 
%%
 
norm_time = t./t_S;
imp1.ind_start = find( norm_time < 0.264 );
imp1.ind_stop = find(norm_time > 0.5326);
 
imp1.t = norm_time(imp1.ind_start(end):imp1.ind_stop(1));
imp1.d2th_B = d2th_B(imp1.ind_start(end):imp1.ind_stop(1));
imp1.d2th_B = imp1.d2th_B./d2th_S;
 
figure(4)
plot(imp1.t,imp1.d2th_B,'g','LineWidth',LW);
ylabel('$\overline{\ddot{\theta}}_{B}$', 'interpreter','latex','FontSize',15);
xlabel('$\overline{t}$','interpreter','latex','FontSize',15);
legend('Prediction');
 
tp = 0.3126;
tw1 = 0.15;
tw2 = 0.04;
 
imp1.ind_on = find(imp1.t == tp);
imp1.ind_off = find(imp1.t == tp+tw1+tw2);
 
eps = log(0.1/tw2);
 
for count = 1:length(imp1.t)
if imp1.t(count)<tp
    w(count) = 0;
elseif tp <= imp1.t(count) && imp1.t(count)< tp+tw1
    w(count) = 1;
elseif imp1.t(count)>= tp+tw1
    w(count) = exp(eps*(imp1.t(count)-tp-tw1));
end
end
 
cap_d2th_B = imp1.d2th_B.*w';
cap_d2th_B_2 = cap_d2th_B.^2;
 
figure(401)
plot(imp1.t,cap_d2th_B_2,'g','LineWidth',LW);
ylabel('$\hat{\overline{\ddot{\theta}}}_{B}$', 'interpreter','latex','FontSize',15);
xlabel('$\overline{t}$','interpreter','latex','FontSize',15);
 
 
Q7 = 0;
 
for i = 1:length(imp1.t)-1
   
    Q7 = Q7+(1/(tw1+tw2))*(cap_d2th_B_2(i)^2)*(imp1.t(i+1)-imp1.t(i));
end
 
figure(402)
subplot(211)
plot(imp1.t,imp1.d2th_B,'g','LineWidth',LW);
ylabel('$\overline{\ddot{\theta}}_{B}$', 'interpreter','latex','FontSize',15);
xlabel('$\overline{t}$','interpreter','latex','FontSize',15);
hold on
plot(tp*[1 1],[-1500 500],'--k')
plot((tw1+tp)*[1 1],[-1500 500],'--k')
plot((tw1+tp+tw2)*[1 1],[-1500 500],'--k')
subplot(212)
plot(imp1.t,cap_d2th_B_2,'g','LineWidth',LW);
ylabel('${\hat{\overline{\ddot{\theta^2}}}_{B}}$', 'interpreter','latex','FontSize',15);
xlabel('$\overline{t}$','interpreter','latex','FontSize',15);
 
%%
 
figure(600)
plot(th_AB(25601:26465),dth_AB(25601:26465),'b','LineWidth',LW);
xlabel('\psi_{AB}','FontSize',15)
ylabel('$\dot\psi_{AB}$','interpreter','latex','FontSize',15);
title('$\overline{t} : 0 - 1$','interpreter','latex');
figure(601)
% hold on
plot(th_AB(26465:27329),dth_AB(26465:27329),'b','LineWidth',LW);
xlabel('\psi_{AB}','FontSize',15)
ylabel('$\dot\psi_{AB}$','interpreter','latex','FontSize',15);
title('$\overline{t} : 1 - 2$','interpreter','latex');
 
%%
 
% DC = load('DC.mat');
% 
% figure(600)
% hold on
% plot(DC.th_AB(25601:26465),DC.dth_AB(25601:26465),'r','LineWidth',LW);
% legend('Single Clearance','Dual Clearance');
% figure(601)
% hold on
% plot(DC.th_AB(26465:27329),DC.dth_AB(26465:27329),'r','LineWidth',LW);
% legend('Single Clearance','Dual Clearance');
% 
% figure(102)
% subplot(311)
% plot(t./t_S,th_B./th_S,'b','LineWidth',LW);
% hold on
% plot(DC.t./t_S,DC.th_B./th_S,'r','LineWidth',LW);
% ylabel('$\overline{{\theta}}_{B}$', 'interpreter','latex','FontSize',15);
% xlim([-1 10]);
% subplot(312)
% plot(t./t_S,dth_B./dth_S,'b','LineWidth',LW);
% hold on
% plot(DC.t./t_S,DC.dth_B./dth_S,'r','LineWidth',LW);
% xlim([-1 10]);
% ylabel('$\overline{\dot{\theta}}_{B}$', 'interpreter','latex','FontSize',15);
% subplot(313)
% plot(t./t_S,d2th_B./d2th_S,'b','LineWidth',LW);
% hold on
% plot(DC.t./t_S,DC.d2th_B./d2th_S,'r','LineWidth',LW);
% xlim([-1 10]);
% ylabel('$\overline{\ddot{\theta}}_{B}$', 'interpreter','latex','FontSize',15);
% legend('Single Clearance','Dual Clearance');
% 
% figure(2)
% hold on
% plot(t/t_S,th_AB/th_S,'-b','LineWidth',LW)
% ylabel('$\overline{{\theta}}_{AB}$', 'interpreter','latex','FontSize',15);
% hold on
% plot(DC.t/t_S,DC.th_AB/th_S,'-r','LineWidth',LW)
% plot(x_lim,p.TH_AB/th_S*[1 1],'--k')
% plot(x_lim,-p.TH_AB/th_S*[1 1],'--k')
% xlim([-0.5 10]);
% legend('Single Clearance','Dual Clearance');
% box on
% xlabel('$\overline{t}$','interpreter','latex','FontSize',15);
% 
% figure(103)
% subplot(211)
% plot(t./t_S,d2th_B./d2th_S,'b','LineWidth',LW);
% xlim([-0.5 8]);
% hold on
% plot(DC.t./t_S,DC.d2th_B./d2th_S,'r','LineWidth',LW);
% ylabel('$\overline{\ddot{\theta}}_{B}$', 'interpreter','latex','FontSize',15);
% box on
% subplot(212)
% plot(t./t_S,KE,'b','LineWidth',LW);
% hold on
% plot(DC.t./t_S,DC.KE_AB,'r','LineWidth',LW);
% xlim([-0.5 8]);
% ylabel('|\Delta KE|_{AB}');
% xlabel('$\overline{t}$','interpreter','latex');
% box on
% legend('Single Clearance','Dual Clearance');
% 
% figure(110)
% subplot(211)
% hold on
% plot(t/t_S,-d2th_A/d2th_S,'-b','LineWidth',LW)
% plot(DC.t/t_S,-DC.d2th_A/d2th_S,'-r','LineWidth',LW)
% ylabel('$\overline{\ddot{\theta}}_{A}$', 'interpreter','latex','FontSize',15);
% legend('Single Clearance ','Dual Clearance');
% xlim([-1 10])
% ylim([-100 100]);
% set(gca,'YTick',[-100 0 100]);
% box on
% subplot(212)
% hold on
% plot(t/t_S,d2th_B/d2th_S,'-b','LineWidth',LW)
% plot(DC.t/t_S,DC.d2th_B/d2th_S,'-r','LineWidth',LW)
% ylabel('$\overline{\ddot{\theta}}_{B}$', 'interpreter','latex','FontSize',15);
% xlabel('$\overline{t}$','interpreter','latex','FontSize',15);
% xlim([-1 10])
% ylim([-1500 1500]);
% set(gca,'YTick',[-1500 0 1500]);
% box on
% 
% figure(402)
% subplot(211)
% plot(imp1.t,imp1.d2th_B,'b','LineWidth',LW);
% hold on
% plot(DC.imp1.t,DC.imp1.d2th_B,'r','LineWidth',LW);
% ylabel('$\overline{\ddot{\theta}}_{B}$', 'interpreter','latex','FontSize',15);
% xlabel('$\overline{t}$','interpreter','latex','FontSize',15);
% xlim([0.264 0.5326]);
% subplot(212)
% plot(imp1.t,cap_d2th_B_2,'b','LineWidth',LW);
% hold on
% plot(DC.imp1.t,DC.cap_d2th_B_2,'r','LineWidth',LW);
% ylabel('${\hat{\overline{\ddot{\theta^2}}}_{B}}$', 'interpreter','latex','FontSize',15);
% xlabel('$\overline{t}$','interpreter','latex','FontSize',15);
% xlim([0.264 0.5326]);
% legend('Single Clearance ','Dual Clearance');
Td_ABbm = Phi_AB + p.h_AB*tanh(50*dth_AB);
th_Abm = th_A;
th_Bbm = th_B;
dth_Abm = dth_A;
dth_Bbm = dth_B;
th_ABbm = th_AB;
d2th_Abm = d2th_A;
d2th_Bbm = d2th_B;
save('values.mat','th_ABbm','d2th_Abm','d2th_Bbm','Td_ABbm','th_Abm','th_Bbm','dth_Abm','dth_Bbm');


