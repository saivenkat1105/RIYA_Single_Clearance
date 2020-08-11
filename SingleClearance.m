format compact;clear all; clc;
close all;

load values.mat  % Load the Benchmark Values

% global A_pks_hi A_pks_low B_pks_hi B_pks_low
% global Ahipks Alowpks Bhipks Blowpks
% 
% Ahipks = zeros(60,5);
% Alowpks = zeros(60,5);
% Bhipks = zeros(60,5);
% Blowpks = zeros(60,5);

%% Multiple values for each parameter

% Uncomment whichever multilier to be needed

% multiplier = [0.1,0.5,2,10,1]; %% Multiplier for all design parameters
stiffsigma = [100000,10000,1000,100]; %% Sigma for Stiffness Smoothening Function
% fricsigma = [1000,100,50,10]; %% Sigma for friction smoothening function
% idbeta = [0.01 0.1]; %% Beta for impact damping
% C_v = [0.001 0.001 0.001]; %% Viscous Damping at AB
% C_v1 = [0.01 0.5 0.1]; %% Viscous Damping at BC

freq = [];
q3B = [];

markers = {'-m','-b','-g','-r','-c'}; %% Colors for the plots

% Loop through all parameters
for i = 1 : size(stiffsigma,2)
   %% All the constants
   %Assign all the constants to the struct p. 
    %m = multiplier(i);  % Muliply the design paramter with this m to vary
    %it and plot
    
    j = i;
    
   % p.beta = betam(i);
    %p.sigma = sigmam(i);
    p.stiffsigma = stiffsigma(i);
    %p.fricsigma = fricsigma(i);
    %p.idbeta = idbeta(i);
    %p.C_v = C_v(i);
    %p.C_v1 = C_v1(i);

    p.J_A = 0.0005080424;    %torsional inertia [kgm^2] (measured from step-response)
    p.J_B = 0.000397876;    %torsional inertia [kgm^2] (measured from impulse-response)
    p.J_C = 0.000397876;
    p.J_e = p.J_A + p.J_B;

    p.k_AB = 8600;          %torsional stiffness [Nm/rad] (static measurement)
    p.k_BC = 47.99;          %torsional stiffness [Nm/rad] (static measurement)
    p.k_C = 8600;

    p.TH_AB = 0.15*pi/180;  %clearance transition [rad] (measured from geometry)
    p.TH_C = 0.1*pi/180;
    p.TH_ABbm = 0.15*pi/180;
    
    p.h_A = 0.033078;       %Coulomb friction [Nm] (measured from step-response)
    p.h_AB = (0.2803-p.h_A)/1.5;  %Coulomb friction [Nm] (0.375-0.625 kg at 0.0762 m*)
    p.h_BC = 0.028075;       %Coulomb friction [Nm] (measured from impulse-response)
    p.h_C = 0.1648;


    p.T_A0 = 2.877*9.81*5.125*2.54/100;  %initial external torque [Nm]
    p.T_Af = 0;                          %final external torque [Nm]

    p.T_B0 = 0;

    tau_s = 1/12.8e3;                   %measurement period
    tau_sbm = 1/12.8e3;                   %measurement period 
    J = diag([p.J_A p.J_B p.J_C]);
    K = [   p.k_AB      -p.k_AB         0;
            -p.k_AB     p.k_AB+p.k_BC   -p.k_BC;
            0           -p.k_BC         p.k_BC+p.k_C;   ];
    lambda = eig(K,J);
    w1 = lambda(1).^0.5; 
    freq = [freq w1];
    t_S = 0.0675;
    
    th_S = p.TH_ABbm;
    dth_S = th_S/t_S;
    d2th_S = dth_S/t_S;

    t_Sbm = 0.0675;
    dth_Sbm = th_S/t_Sbm;
    d2th_Sbm = dth_Sbm/t_Sbm;
    
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
    td = -2:tau_s:2;
    th = -0.5*pi/180:pi/180000:0.5*pi/180;
    t_cpu = 7.40625;

    LW = 2;

    %% Solve the equations 
    options = odeset('RelTol',1e-3,'AbsTol',1e-3,'Maxstep',tau_s);
    t_span1 = -2:tau_s:2;
    tbm = -2:tau_sbm:2;
    t_span2 = -2:tau_s/10:2;

    X_0 = [th_A0 th_B0 0 0];    %initial conditions
    %  X_0 = [0 0 0 0];         %initial conditions
    
    %% ode45
    a = cputime;
    [t45,X45] = ode45(@EOM_X1_1a,t_span1,X_0,options,p);    
    ode45time = (cputime - a)/t_cpu;
    
    %% ode23s
    % a = cputime;
    % [t23s,X23s] = ode23s(@EOM_X1_1a,t_span1,X_0,options,p);
    % ode23stime = (cputime - a)/t_cpu;
    
    %% ode23t
    
    % a = cputime;
    % [t23t,X23t] = ode23t(@EOM_X1_1a,t_span1,X_0,options,p);
    % ode23ttime = cputime - a;
    
    %% ode23tb 
    % a = cputime;
    % [t23tb,X23tb] = ode23tb(@EOM_X1_1a,t_span1,X_0,options,p);
    % ode23tbtime = cputime - a;
    
    %% ode15s
    % a = cputime;
    % [t15s,X15s] = ode15s(@EOM_X1_1a,t_span2,X_0,options,p);
    % ode23tbtime = cputime - a;

    % X = [X45 , X23s, X23t, X23tb, X15s];
    X = X45;

    %times = [ode45time, ode23stime, ode23ttime, ode23tbtime, ode15stime];
    %% Results
    T_A = p.T_A0*(1-usgn(p,t45))+p.T_Af*usgn(p,t45); %external torque [Nm]
%     T_A = p.T_A0*(1-usgn(t45)+exp(-t45/p.beta).*usgn(t45));
    
    T_B = p.T_B0;
    th_A = [];
    th_B = [];
    dth_A = [];
    dth_B = [];

    %Results
    % Loop through this if multiple solvers are used
    for k = 0 
        th_A = [th_A,X(:,4*k+1)];          %angular displacement [rad]
        th_B = [th_B,X(:,4*k+2)];          %angular displacement [rad]+
        th_AB = th_A-th_B;      %relative angular displacement [rad]
        dth_A = [dth_A,X(:,4*k+3)];         %angular velcoity [rad/s]
        dth_B = [dth_B,X(:,4*k+4)];         %angular velocity [rad/s]
        dth_AB = dth_A-dth_B;   %relative angular velocity [rad/s]
    end


    Phi_AB =    p.k_AB.*(th_AB+p.TH_AB).*(1-u(p,th_AB+p.TH_AB))+...
                p.k_AB.*(th_AB-p.TH_AB).*u(p,th_AB-p.TH_AB);          %elastic torque through clearance [Nm]


    Td_AB = Phi_AB + p.h_AB*ufric(p,dth_AB);
    T0 = p.T_A0;
    d2th_A = (T_A-Phi_AB-p.h_A*ufric(p,dth_A)-p.h_AB*ufric(p,dth_AB))/p.J_A;             %angular acceleration [rad/s^2]
    d2th_B = (T_B+Phi_AB-p.k_BC*th_B+p.h_AB*ufric(p,dth_AB)-p.h_BC*ufric(p,dth_B))/p.J_B;  %angular acceleration [rad/s^s]
    t = t45;

    % Impact damping torque with hysteresis at AB
    %Tg_h = Phi_AB.*(1+p.beta*(dth_A-dth_B))+p.h_AB*ufric(p,dth_A-dth_B);
    
    % Impact damping torque at AB
    %Tg = Phi_AB.*(1+p.beta*(dth_A-dth_B));
    
    % Viscous dapming with hysteresis at AB
    %Tv = Phi_AB+ p.C_v*(dth_A-dth_B)+p.h_AB*ufric(p,(dth_A-dth_B));
    
    %% Plot the input torque
    
%     figure(1)
%     plot(t./p.t_S,T_A,'g-','LineWidth',2);
%     title('Input: Step Down Torque');
%     xlabel('$\overline{t}[s]$','interpreter','latex','FontSize',15);
%     ylabel('T [Nm]');


 %%   
     figure(100)
     hold on
     plot(th_AB/th_S,Td_AB,markers{i},'LineWidth',LW/2)
     xlabel('$\overline{{\theta}}_{AB}$', 'interpreter','latex','FontSize',30);
     ylabel('$\overline{Td}$','interpreter','latex','FontSize',30);
     title('Damping Torque vs $\overline{{\theta}}_{AB}$','interpreter','latex','FontSize',30)
     %legend('ode45','ode23s','interpreter','latex','FontSize',30)
     lgd = legend('1000','100','50','10','interpreter','latex','FontSize',25);
     title(lgd,'$\sigma$','interpreter','latex','FontSize',25)
     ax = gca;
     ax.YAxis.FontSize = 20;
     ax.XAxis.FontSize = 20;
     
     x_lim = [-0.5 20];
    
    %% Plot relative angular dispacemet with time
    figure(60)
    subplot(211)
    hold on
    plot(tbm./t_Sbm,th_ABbm/th_S,'-b','LineWidth',LW)
    % plot(t./t_S,th_AB(:,2)/th_S,'LineWidth',LW)
    xlim([-1 20])
    yline(p.TH_ABbm/th_S)
    yline(-p.TH_ABbm/th_S)
    legend('Benchmark')
    
    subplot(212)
    hold on
    plot(tbm./t_Sbm,th_ABbm/th_S,'-g','LineWidth',LW)
    
    p1 = plot(t./t_S,th_AB(:,1)/th_S,markers{i},'LineWidth',LW);
    %     plot(t./t_S,th_AB(:,3)/th_S,'-y','LineWidth',LW)
    %     plot(t./t_S,th_AB(:,4)/th_S,'-k','LineWidth',LW)
    %     plot(t./t_S,th_AB(:,5)/th_S,'-m','LineWidth',LW)
    ylabel('$\overline{{\theta}}_{AB}$', 'interpreter','latex','FontSize',30);
    xlabel('$\overline{t}$','interpreter','latex','FontSize',30);
    plot(x_lim,p.TH_AB/th_S*[1 1],'--k')
    plot(x_lim,-p.TH_AB/th_S*[1 1],'--k')
    xlim([-1 20])
    %xlim([-0.00005 0.00005])
    box on
    lgd = legend('100000','10000','1000','100','interpreter','latex','FontSize',25);
    title(lgd,'\sigma','interpreter','tex','FontSize',25)
    title('$\overline{\theta}_{AB}$ vs $\overline{t} $ - Arc Tangent Approximation','interpreter','latex','FontSize',25)
    ax = gca;
    ax.YAxis.FontSize = 20;
    ax.XAxis.FontSize = 20;

    %% Find teh Peaks in the accelearion plot
    %  
    [A_pks_hi,A_loc_hi] = findpeaks(d2th_A./d2th_S,t./t_S,'MinPeakHeight',0,'MinPeakDistance',0.42);
    [A_pks_low,A_loc_low] = findpeaks(-d2th_A./d2th_S,t./t_S,'MinPeakHeight',0,'MinPeakDistance',0.42);
    [B_pks_hi,B_loc_hi] = findpeaks(d2th_B./d2th_S,t./t_S,'MinPeakHeight',50,'MinPeakDistance',0.3);
    [B_pks_low,B_loc_low] = findpeaks(-d2th_B./d2th_S,t./t_S,'MinPeakHeight',70,'MinPeakDistance',0.3);
 

    a = linspace(1,15,15);
    b = linspace(2,16,15);

    A_ind_hi = find(A_loc_hi>0);
    A_ind_low = find(A_loc_low>0);
    B_ind_hi = find(B_loc_hi>0);
    B_ind_low = find(B_loc_low>0);

    A_loc_hi = A_loc_hi(A_loc_hi>0);
    A_loc_low = A_loc_low(A_loc_low>0);
    B_loc_hi = B_loc_hi(B_loc_hi>0);
    B_loc_low = B_loc_low(B_loc_low>0);

    A_pks_hi = A_pks_hi(A_ind_hi);
    A_pks_low = A_pks_low(A_ind_low);
    B_pks_hi = B_pks_hi(B_ind_hi);
    B_pks_low = B_pks_low(B_ind_low);
    
    % Calculate q3B values (Use the indexfor A_pks and B_pks appropriately)
    val1 = bdb(B_pks_hi(1),B_pks_low(2));
    val2 = bdb(B_pks_hi(1),B_pks_low(1));
    val = max(val1,val2);
    q3B = [q3B, val];

%% Only Acceleration plots at A and B
%     
    figure(20)

    subplot(211)
    hold on
 
    plot(tbm./t_Sbm,d2th_Bbm(:,1)./d2th_Sbm,'LineWidth',LW);
    plot(t./t_S,d2th_B(:,1)./d2th_S,markers{i},'LineWidth',LW);
    
    % Plotting the peaks and lows
    % plot(B_loc_hi,B_pks_hi,'bO');
    % plot(B_loc_low,-B_pks_low,'rO');
    
    xlim([-1 10]);
    ylabel('$\overline{\ddot{\theta}}_{B}$', 'interpreter','latex','FontSize',30);
    %legend('Benchmark','Current Case','interpreter','latex','FontSize',20)
    %legend('Current Case','Benchmark','interpreter','latex','FontSize',20,'Location','best')
    title('$\overline{\ddot{\theta}}_{B}$ vs $\overline{t}$')
    ax = gca;
    ax.YAxis.FontSize = 20;
    ax.XAxis.FontSize = 20;
    
    subplot(212)
    hold on
    plot(tbm./t_Sbm,d2th_Abm(:,1)./d2th_Sbm,'LineWidth',LW);
    plot(t./t_S,d2th_A(:,1)./d2th_S,markers{i},'LineWidth',LW);
    plot(A_loc_hi,A_pks_hi,'bO');
    plot(A_loc_low,-A_pks_low,'rO');
    xlim([-1 10]);
    ylabel('$\overline{\ddot{\theta}}_{A}$', 'interpreter','latex','FontSize',30);
    x_lim = [-0.5 20];
    lgd = legend('1000','100','50','10','interpreter','latex','FontSize',25);
    title(lgd,'$\sigma$','interpreter','latex','FontSize',25)
    ax = gca;
    ax.YAxis.FontSize = 20;
    ax.XAxis.FontSize = 20;
    title('$\overline{\ddot{\theta}}_{A}$ vs $\overline{t}$')
    
    
    %% All theta plots
%     figure(20+j)
%     subplot(231)
%     hold on
%     % plot(t./t_S,th_B./th_S,'g','LineWidth',LW);
%     plot(t./t_S,th_Bbm./th_S,'b','LineWidth',LW);
%     plot(t./t_S,th_B./th_S,'g','LineWidth',LW);
%     
%     ylabel('$\overline{{\theta}}_{B}$', 'interpreter','latex','FontSize',30);
%     xlim([-1 10]);
%     subplot(232)
%     hold on
%     % plot(t./t_S,dth_B./dth_S,'g','LineWidth',LW);
%     plot(t./t_S,dth_Bbm./dth_S,'b','LineWidth',LW);
%     plot(t./t_S,dth_B./dth_S,'g','LineWidth',LW);
%     xlim([-1 10]);
%     ylabel('$\overline{\dot{\theta}}_{B}$', 'interpreter','latex','FontSize',30);
%     subplot(233)
% %     subplot(211)
%     hold on
% %     plot(t./t_S,d2th_B(:,1)./d2th_S,'g','LineWidth',LW);
%     plot(tbm./t_Sbm,d2th_Bbm(:,1)./d2th_Sbm,'b','LineWidth',LW);
%     plot(t./t_S,d2th_B(:,1)./d2th_S,'g','LineWidth',LW);
%     plot(B_loc_hi,B_pks_hi,'bO');
%      plot(B_loc_low,-B_pks_low,'rO');
%     xlim([-1 10]);
%     ylabel('$\overline{\ddot{\theta}}_{B}$', 'interpreter','latex','FontSize',30);
%     legend('Benchmark','Current Case','interpreter','latex','FontSize',20)
% %     legend('Current Case','Benchmark','interpreter','latex','FontSize',20,'Location','best')
% 
%     subplot(234)
%     hold on
%     % plot(t./t_S,th_A./th_S,'g','LineWidth',LW);
%     plot(t./t_S,th_Abm./th_S,'b','LineWidth',LW);
%     plot(t./t_S,th_A./th_S,'g','LineWidth',LW);
%     
%     ylabel('$\overline{{\theta}}_{A}$', 'interpreter','latex','FontSize',30);
%     xlim([-1 10]);
%     subplot(235)
%     hold on
%     % plot(t./t_S,dth_A./dth_S,'g','LineWidth',LW);
%     plot(t./t_S,dth_Abm./dth_S,'b','LineWidth',LW);
%     plot(t./t_S,dth_A./dth_S,'g','LineWidth',LW);
%     xlim([-1 10]);
%     ylabel('$\overline{\dot{\theta}}_{A}$', 'interpreter','latex','FontSize',30);
%     subplot(236)
% %     subplot(212)
%     hold on
% %     plot(t./t_S,d2th_A(:,1)./d2th_S,'g','LineWidth',LW);
%     plot(tbm./t_Sbm,d2th_Abm(:,1)./d2th_Sbm,'b','LineWidth',LW);
%     plot(t./t_S,d2th_A(:,1)./d2th_S,'g','LineWidth',LW);
%     plot(A_loc_hi,A_pks_hi,'bO');
%     plot(A_loc_low,-A_pks_low,'rO');
%     xlim([-1 10]);
%     ylabel('$\overline{\ddot{\theta}}_{A}$', 'interpreter','latex','FontSize',30);
%     x_lim = [-0.5 20];


 %% Plot Phi_AB
    % 
    % figure(40)
    % hold on
    % plot(t./t_S,Phi_AB,'-g','LineWidth',LW)
    % ylabel('$\overline{{\phi}}_{AB}$', 'interpreter','latex','FontSize',30);
    % xlabel('$\overline{t}$','interpreter','latex','FontSize',30);
    % xlim([-0.5 0.8])
    % 
%% Plot  k_AB
   
    % k1_AB =    p.k_AB.*(th+p.TH_AB).*(1-u(th+p.TH_AB))+...
    %             p.k_AB.*(th-p.TH_AB).*u(th-p.TH_AB);
    %         
    % %k_AB = p.k_AB.*uhc(th,p.TH_AB);     
    % %Phi_AB =  sign(th).*p.k_AB.*(abs(th)-p.TH_AB).*uhc(abs(th)-p.TH_AB);
    % 
    % figure(50)
    % hold on
    % %plot(th,Phi_AB,'-b','LineWidth',2)
    % plot(th,k1_AB,'-g','LineWidth',2)
    % 
    % ylabel('$\overline{{k}}_{AB}$', 'interpreter','latex','FontSize',30);
    % xlabel('$\overline{theta}$','interpreter','latex','FontSize',30);
    % xline(p.TH_AB)
    % xline(-p.TH_AB)
    
%% Comparison with benchmark
    % t = t45;
    % x_lim = [-0.5 20];
    % figure(60)
    % % subplot(411)
    % hold on
    % plot(tbm./t_S,th_ABbm/th_S,'-g','LineWidth',LW)
    % plot(t./t_S,th_AB(:,1)/th_S,'-b','LineWidth',LW)
    % plot(t./t_S,th_AB(:,2)/th_S,'-r','LineWidth',LW)
    % plot(t./t_S,th_AB(:,3)/th_S,'-y','LineWidth',LW)
    % plot(t./t_S,th_AB(:,4)/th_S,'-k','LineWidth',LW)
    % plot(t./t_S,th_AB(:,5)/th_S,'-m','LineWidth',LW)
    % ylabel('$\overline{{\theta}}_{AB}$', 'interpreter','latex','FontSize',30);
    % xlabel('$\overline{t}$','interpreter','latex','FontSize',30);
    % plot(x_lim,p.TH_AB/th_S*[1 1],'--k')
    % plot(x_lim,-p.TH_AB/th_S*[1 1],'--k')
    % xlim([-0.5 20])
    % box on
    % title('Friction Non Linearity - Hyp Tan with \sigma = 1000')
    % legend('Benchmark','ode45','ode15s','ode23s','ode23s','ode23t','ode23tb')
    % 
    % %%
    % subplot(412)
    % hold on
    % plot(th_ABbm/th_S,Td_ABbm/T0,'-g','LineWidth',LW)
    % plot(th_AB/th_S,Td_AB/T0,'-b','LineWidth',LW)
    % 
    % xlabel('$\overline{{\theta}}_{AB}$', 'interpreter','latex','FontSize',30);
    % ylabel('$\overline{T}_{dAB}$','interpreter','latex','FontSize',30);
    % 
    % plot(p.TH_AB/th_S*[1 1],[-3 3],'--k')
    % plot(-p.TH_AB/th_S*[1 1],[-3 3],'--k')
    % 
    % subplot(413)
    % hold on
    % plot(tbm./t_S,d2th_Bbm./d2th_S,'g','LineWidth',LW);
    % plot(t./t_S,d2th_B./d2th_S,'b','LineWidth',LW);
    % 
    % xlim([-1 10]);
    % % % ylim([-10e4 10e4]);
    % ylabel('$\overline{\ddot{\theta}}_{B}$', 'interpreter','latex','FontSize',30);
    % xlabel('$\overline{t}$','interpreter','latex','FontSize',30);
    % 
    % 
    % subplot(414)
    % hold on
    % plot(tbm./t_S,d2th_Abm./d2th_S,'g','LineWidth',LW);
    % plot(t./t_S,d2th_A./d2th_S,'b','LineWidth',LW);
    % 
    % xlim([-1 10]);
    % ylabel('$\overline{\ddot{\theta}}_{A}$', 'interpreter','latex','FontSize',30);
    % xlabel('$\overline{t}$','interpreter','latex','FontSize',30);
    % 
    % 
 %  close all
end
%% Plot variation of q3B values with various parameters.
% q3Bhyp =  q3B;
% betahyp = betam;
% figure()
% hold on
% plot(stiffsigma,q3B,'.-r','MarkerSize',15,'LineWidth',LW)
% xlabel('$\sigma$','interpreter','latex','FontSize',30)
% ylabel('$Q_{3}B$','interpreter','latex','FontSize',30)
% save('Q3Bhyp.mat','q3Bhyp','betahyp');
