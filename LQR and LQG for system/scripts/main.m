clear all; clc;
% Defining variables for the system
syms M m_1 l_1 m_2 l_2 g 

%Values for system variables
M = 1000;
m_1= 100;
l_1= 100;
m_2 = 20;
l_2= 10;
g = 9.8;
tspan = 0:0.1:200; % time span for simulation
y0 = [1; 0; 0; 0; 0; 0];  % Initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FIRST COMPONENT %%%%

%% PART A, B and C ==> Calculations in Report 

%% PART D 
% Linearized model state variable matrix
A=[0 1 0 0 0 0 ; 0 0 -(m_1*g)/M 0 -(m_2*g)/M 0 ; 0 0 0 1 0 0 ; 0 0 -(M + m_1)*g/(M*l_1) 0 -(m_2*g)/(M*l_1) 0 ; 0 0 0 0 0 1;
    0 0 -(m_1*g)/(M*l_2) 0 -(M + m_2)*g/(M*l_2) 0];
B=[0 ; 1/M; 0; 1/(M*l_1) ; 0 ; 1/(M*l_2)];
C=[1 0 0 0 0 0 ; 0 0 1 0 0 0 ; 0 0 0 0 1 0];
D=[0 ; 0 ; 0];
sys_1 = ss(A,B,C,D);

% Controllability Matrix
ctrb_matrix = ctrb(A,B);
rank_ctrb_matrix = rank(ctrb_matrix);  % Rank = number of columns(6) ==> controllable

% Observability Matrix
obsv_matrix = obsv(A,C);
rank_obsv_matrix = rank(obsv_matrix);  % Rank = number of columns(6) ==> observable

%% LQR Controller implementation 

% R=0.0001;
% Q = diag([100,10,5000000,1000,1000000,1000]); %Most stable 
F =0.5*ones(size(tspan));
Q = [5 0 0 0 0 0; 0 0 0 0 0 0; 0 0 5000 0 0 0; 0 0 0 0 0 0; 0 0 0 0 5000 0; 0 0 0 0 0 0];
R = 0.001;
% Q(1,1)=1000;
% Q(3,3)=50;
% Q(5,5)=30;

%% LQR controller for Linear System 
[K_closed,S,e_closed] = lqr(A,B,Q,R); % Returns closed loop gains, S (in ricatti eqn) and closed loop eigen values
A_closed = A - B*K_closed;

[t,y1] = ode45(@(t,y)linear_system(m_1, m_2, M, l_1, l_2, y, -K_closed*y) ,tspan,y0);
figure;
hold on
plot(t,y1(:,1),'g')
plot(t,y1(:,3),'b')
plot(t,y1(:,5),'r')
ylabel('state variables')
xlabel('time(sec)')
title('LQR response for linear system')
legend('x_{cart}','theta1','theta2')


%% LQR controller for Non-Linear System 
[t,y2] = ode45(@(t,y)nonlinear_system(t, y, -K_closed*y, m_1, m_2, l_1, l_2, M),tspan,y0);
figure(2);
hold on
plot(t,y2(:,1))
plot(t,y2(:,3))
plot(t,y2(:,5))
ylabel('state variables')
xlabel('time (sec)')
title('LQR response for non - linear system')
legend('x_{cart}','theta1','theta2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SECOND COMPONENT %%%%
%% PART E 
C1 = [1 0 0 0 0 0];    %For output (x(t))
C2 = [0 0 1 0 0 0; 0 0 0 0 1 0];   %For output (theta_1(t),theta_2(t))       
C3 = [1 0 0 0 0 0;  0 0 0 0 1 0];  %For output (x(t),theta_2(t))        
C4 = [1 0 0 0 0 0;  0 0 1 0 0 0;    0 0 0 0 1 0]; %For output (x(t),theta_1(t),theta_2(t))

% Check observability for each case of C
obsv_matrix_C1 = obsv(A,C1);  
rank_obsv_C1 = rank(obsv_matrix_C1); % output = 6 ==> observable
obsv_matrix_C2 = obsv(A,C2);
rank_obsv_C2 = rank(obsv_matrix_C2); % output = 4 ==> Not observable
obsv_matrix_C3 = obsv(A,C3);
rank_obsv_C3 = rank(obsv_matrix_C3); % output = 6 ==> observable
obsv_matrix_C4 = obsv(A,C4);
rank_obsv_C4 = rank(obsv_matrix_C4); % output = 6 ==> observable

%% PART F: Plot the Luenberger Response for the output vectors which were observable
%% Observer based Control
%% Kalman Estimator Design
%% Obtaining Luenberger Observers for output vectors of each observable(Using Kalman Bucy Filter (Lqe))
%Introducing Noise and Disturbances in the system
Bd = 0.01.*eye(6);             %input disturbance covariance
Bn = 0.1;
Bn1 = 0.1;                      %output measurement noise for each case
Bn3 = 0.1*[0,1;0,1];
Bn4 = 0.1*[0,0,1;0,0,1;0,0,1];

%%
% system
[Lb1] = lqe(A,Bd,C1,Bd,Bn);
[Lb3] = lqe(A,Bd,C3,Bd,Bn*eye(2));
[Lb4] = lqe(A,Bd,C4,Bd,Bn*eye(3));

% Creating Augmented Matrices for Simulation
uD = randn(6,size(tspan,2));      %input for disturbance
uN = randn(size(tspan));          %input for noise
u = 0*tspan;
u(200:length(tspan)) = 1;      % Step input at t = 10
u1 = [u; Bd*Bd*uD; uN];

Be = [B,Bd,zeros(size(B))];  %Augmented B matrix

%% Luenberger Observer with State Estimator for C1 = [1,0,0,0,0,0] Kalman Filter

sysLO1 = ss(A-Lb1*C1,[B Lb1],C1,zeros(1,2));     %State Estimator system

%Obtaining Y values for a system simulated with noise and disturbance. 
De1 = [0, 0 ,0, 0, 0, 0, 0, Bn1];                     %Augmented D matrix

sys1 = ss(A,Be,C1,De1);
[y1,t] = lsim(sys1,u1,tspan);

%Simulating the States of the output variables 
[x1,t] = lsim(sysLO1,[u; y1'],tspan);

figure();
hold on
plot(t,y1(:,1),'r','Linewidth',1)
plot(t,x1(:,1),'y--','Linewidth',1)
ylabel('x-position of cart')
xlabel('time in s')
legend('Output obtained from noisy system','Estimated output of the system')
title('Estimated Response for C1: output vector x(t) - Linear System')
hold off

% opt = simset('solver','ode45','SrcWorkspace','Current');
% [tout2]=sim('nonlinerLO',tspan,opt);
% figure();
% hold on
% plot(tout2,out1(:,1),'r','Linewidth',2)
% plot(tout2,states1(:,1),'k--','Linewidth',2)
% ylabel('x-position of cart')
% xlabel('time in s')
% legend('Output obtained from noisy system','Estimated output of the system')
% title('Estimated Response for C1: output vector x(t) - Nonlinear System')
% hold off

%% Luenberger Observer with State Estimator for C3: For output (x(t),t2(t))- Kalman Filter

sysLO3 = ss(A-Lb3*C3,[B Lb3],C3,zeros(2,3));     %State Estimator system

%Obtaining Y values for a system simulated with noise and disturbance. 
De3 = [zeros(size(C3)),Bn3];                     %Augmented D matrix

sys3 = ss(A,Be,C3,De3);
[y3,t] = lsim(sys3,u1,tspan);

%Simulating the States of the output variables 
[x3,t] = lsim(sysLO3,[u; y3'],tspan);

figure();
hold on
plot(t,y3(:,1),'g','Linewidth',2)
plot(t,y3(:,2),'b','Linewidth',2)
plot(t,x3(:,1),'k--','Linewidth',1)
plot(t,x3(:,2),'r--','Linewidth',1)
ylabel('State Variables ')
xlabel('time in s')
legend('Noisy output x(t)','Noisy output theta2(t)','Estimated x(t)','Estimated theta2(t)')
title('Estimated Response for C3: output vector (x(t),t2(t))')
hold off

opt = simset('solver','ode45','SrcWorkspace','Current');
[tout3]=sim('nonlinearLO3',tspan,opt);
figure();
hold on
plot(tout3,out3(:,1),'r','Linewidth',2)
plot(tout3,out3(:,2),'g','Linewidth',2)
plot(t,states3(:,1),'k--','Linewidth',1)
plot(t,states3(:,2),'r--','Linewidth',1)
ylabel('x-position of cart')
xlabel('time in s')
legend('Noisy output x(t)','Noisy output theta2(t)','Estimated x(t)','Estimated theta2(t)')
title('Estimated Response for C3: output vector (x(t),t2(t)) - Nonlinear System')
hold off

%% Luenberger Observer with State Estimator for C4: For output (x(t),t1(t),t2(t)) - Kalman Filter

sysLO4 = ss(A-Lb4*C4,[B Lb4],C4,zeros(3,4));     %State Estimator system

%Obtaining Y values for a system simulated with noise and disturbance. 
De4 = [zeros(3,5),Bn4];                     %Augmented D matrix

sys4 = ss(A,Be,C4,De4);
[y4,t] = lsim(sys4,u1,tspan);

%Simulating the States of the output variables 
[x4,t] = lsim(sysLO4,[u;y4'],tspan);

figure();
hold on
plot(t,y4(:,1),'g','Linewidth',2)
plot(t,y4(:,2),'b','Linewidth',3)
plot(t,y4(:,3),'c','Linewidth',2)
plot(t,x4(:,1),'m--','Linewidth',1)
plot(t,x4(:,2),'r--','Linewidth',2)
plot(t,x4(:,3),'k--','Linewidth',1)
ylabel('State Variables ')
xlabel('time in s')
legend('Noisy output x(t)','Noisy output theta1(t)','Noisy output theta2(t)','Estimated x(t)','Estimated theta1(t)','Estimated theta2(t)')
title('Estimated Response for C4: output vector (x(t),t1(t),t2(t))')
hold off

opt = simset('solver','ode45','SrcWorkspace','Current');
[tout4]=sim('nonlinearLO4',tspan,opt);

figure();
hold on
plot(tout4,out4(:,1),'g','Linewidth',2)
plot(tout4,out4(:,2),'b','Linewidth',3)
plot(tout4,out4(:,3),'c','Linewidth',2)
plot(tout4,states4(:,1),'m--','Linewidth',1)
plot(tout4,states4(:,2),'r--','Linewidth',2)
plot(tout4,states4(:,3),'k--','Linewidth',1)
ylabel('State Variables ')
xlabel('time in s')
legend('Noisy output x(t)','Noisy output theta1(t)','Noisy output theta2(t)','Estimated x(t)','Estimated theta1(t)','Estimated theta2(t)')
title('Estimated Response for C4: output vector (x(t),t1(t),t2(t))')
hold off

%% LQG Controller for Output Vector C1 = [1,0,0,0,0,0]

Ac = A-Lb1*C1;
Bc = [B Lb1];
Cc = eye(6);
Dc = 0*[B Lb1];

opt = simset('solver','ode45','SrcWorkspace','Current');
[tout]=sim('nonlinearlqg',tspan,opt);
