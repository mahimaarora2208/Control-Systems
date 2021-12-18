
%{
The code for C part to check the conditions of controllability
%}
syms M m1 m2 l1 l2 g

A=[0 1 0 0 0 0;0 0 -m1*g/M 0 -m2*g/M 0;0 0 0 1 0 0;0 0 -(M+m1)*g/(M*l1) 0 -m2*g/(M*l1) 0;0 0 0 0 0 1;
    0 0 (-m1*g)/(M*l2) 0 -(M+m2)*g/(M*l2) 0];
B=[0;1/M;0;1/(M*l1);0;1/(M*l2)];

%{
The code for D part, to check controllability, to plot the response for LQR
controller and to check stability using Lyapunov's Indirectt method.
%}
M=1000;
m1=100;
m2=100;
l1=20;
l2=10;
g=9.8;

A=[0 1 0 0 0 0;0 0 -m1*g/M 0 -m2*g/M 0;0 0 0 1 0 0;0 0 -(M+m1)*g/(M*l1) 0 -m2*g/(M*l1) 0;0 0 0 0 0 1;
    0 0 (-m1*g)/(M*l2) 0 -(M+m2)*g/(M*l2) 0];
B=[0;1/M;0;1/(M*l1);0;1/(M*l2)];
Co=[B A*B A*A*B A*A*A*B A*A*A*A*B A*A*A*A*A*B];
det(Co);


Co=[B A*B A*A*B A*A*A*B A*A*A*A*B A*A*A*A*A*B];
rank(Co);

C=[1 0 0 0 0 0;0 0 1 0 0 0;0 0 0 0 1 0];
D=[0;0;0];

Q=C'*C;
R=1;

Q(1,1)=9000;
Q(3,3)=30;
Q(5,5)=30;

Kc=lqr(A,B,Q,R);

Ac = [(A-B*Kc)];
Bc = [B];
Cc = [C];
Dc = [D];
states = {'x' 'x_dot' 'phi1' 'phi1_dot' 'phi2' 'phi2_dot'};
inputs = {'F'};
outputs = {'x'; 'phi1'; 'phi2'};
sys_cl = ss(Ac,Bc,Cc,Dc,'statename',states,'inputname',inputs,'outputname',outputs);
t = 0:0.01:2000;
F =0.2*ones(size(t));
[y,t,x]=lsim(sys_cl,F,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','cart position (m)')
set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)')
title('Step Response with LQR Control')

Ac = [(A-B*Kc)];
Bc = [B];
Cc = [C];
Dc = [D];
states = {'x' 'x_dot' 'phi1' 'phi1_dot' 'phi2' 'phi2_dot'};
inputs = {'F'};
outputs = {'x'; 'phi1'; 'phi2'};
sys_cl = ss(Ac,Bc,Cc,Dc,'statename',states,'inputname',inputs,'outputname',outputs);
t = 0:0.01:1200;
F =0.2*ones(size(t));
[y,t,x]=lsim(sys_cl,F,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,3),'plot');
set(get(AX(1),'Ylabel'),'String','cart position (m)')
set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)')
title('Step Response with LQR Control')

poles=eig(Ac);



%{ 
The code for E part to check the observability for each of the output vectors
%}
%}
C1=C(1:1,:);
O1=[C1' A'*C1' (A'^2)*C1' (A'^3)*C1' (A'^4)*C1' (A'^5)*C1'];
rank(O1);
C2=C(2:3,:);
O2=[C2' A'*C2' (A'^2)*C2' (A'^3)*C2' (A'^4)*C2' (A'^5)*C2'];
rank(O2);
C3=[C(1:1,:);C(3:3,:)];
O3=[C3' A'*C3' (A'^2)*C3' (A'^3)*C3' (A'^4)*C3' (A'^5)*C3'];
rank(O3);
C4=C;
O4=[C4' A'*C4' (A'^2)*C4' (A'^3)*C4' (A'^4)*C4' (A'^5)*C4'];
rank(O4);

%{
The code for F part, to plot the Luenberger Response for the output vectors
which were observable.
%}

P=[-0.04 -0.05 -0.06 -0.07 -0.08 -0.09];

%{
Case 1: When the output vector is X.
%}

L1= place(A',C1',P)'

Ace = [(A-B*Kc) (B*Kc);
       zeros(size(A)) (A-L1*C1)];
Bce = [B;
       zeros(size(B))];
Cce = [C1 zeros(size(C1))];
Dce = [0];

states = {'x' 'x_dot' 'phi1' 'phi1_dot' 'phi2' 'phi2_dot' 'e1' 'e2' 'e3' 'e4' 'e5' 'e6'};
inputs = {'F'};
outputs = {'x'};

sys_est_cl = ss(Ace,Bce,Cce,Dce,'statename',states,'inputname',inputs,'outputname',outputs);

t = 0:0.01:1500;
F = 0.2*ones(size(t));
[y,t,x]=lsim(sys_est_cl,F,t);
plot(t,y(:,1))
title('Step Response with Observer-Based State-Feedback Control')

%{
Case 2: When the output vector is X and theta 2
%}

L3= place(A',C3',P)'
Ace = [(A-B*Kc) (B*Kc);
       zeros(size(A)) (A-L3*C3)];
Bce = [B;
       zeros(size(B))];
Cce = [C3 zeros(size(C3))];
Dce = [0;0];

states = {'x' 'x_dot' 'phi1' 'phi1_dot' 'phi2' 'phi2_dot' 'e1' 'e2' 'e3' 'e4' 'e5' 'e6'};
inputs = {'F'};
outputs = {'x'; 'phi2'};

sys_est_cl = ss(Ace,Bce,Cce,Dce,'statename',states,'inputname',inputs,'outputname',outputs);

t = 0:0.01:1500;
F = 0.2*ones(size(t));
[y,t,x]=lsim(sys_est_cl,F,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','cart position (m)')
set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)')
title('Step Response with Observer-Based State-Feedback Control')

%{
Case 3: When the output vector is X, theta1 and theta 2
%}
L4= place(A',C4',P)'
Ace = [(A-B*Kc) (B*Kc);
       zeros(size(A)) (A-L4*C4)];
Bce = [B;
       zeros(size(B))];
Cce = [C4 zeros(size(C4))];
Dce = [0;0;0];

states = {'x' 'x_dot' 'phi1' 'phi1_dot' 'phi2' 'phi2_dot' 'e1' 'e2' 'e3' 'e4' 'e5' 'e6'};
inputs = {'F'};
outputs = {'x'; 'phi1'; 'phi2'};

sys_est_cl = ss(Ace,Bce,Cce,Dce,'statename',states,'inputname',inputs,'outputname',outputs);

t = 0:0.01:1500;
F = 0.2*ones(size(t));
[y,t,x]=lsim(sys_est_cl,F,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','cart position (m)')
set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)')
title('Step Response with Observer-Based State-Feedback Control')

L4= place(A',C4',P)'
Ace = [(A-B*Kc) (B*Kc);
       zeros(size(A)) (A-L4*C4)];
Bce = [B;
       zeros(size(B))];
Cce = [C4 zeros(size(C4))];
Dce = [0;0;0];

states = {'x' 'x_dot' 'phi1' 'phi1_dot' 'phi2' 'phi2_dot' 'e1' 'e2' 'e3' 'e4' 'e5' 'e6'};
inputs = {'F'};
outputs = {'x'; 'phi1'; 'phi2'};

sys_est_cl = ss(Ace,Bce,Cce,Dce,'statename',states,'inputname',inputs,'outputname',outputs);

t = 0:0.01:1500;
F = 0.2*ones(size(t));
[y,t,x]=lsim(sys_est_cl,F,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,3),'plot');
set(get(AX(1),'Ylabel'),'String','cart position (m)')
set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)')
title('Step Response with Observer-Based State-Feedback Control')

%{
The code for F part, to plot the Luenberger Observer Response(non-linear system) for the output vectors
which were observable.
%}

R = 0.00001;
Q= C'*C;
Q(1,1) =450;
Q(2,2) =500;
Q(3,3) =480;
Q(4,4) =550;
Q(5,5) =460;
Q(6,6) =500;
K = lqr(A,B,Q,R);

r = 0.00001;
w = 250*eye(6);
v = 0.00005*eye(6);
u = w*v*w';
sol = care(A',C1',u,r);
L1 = sol*C1'*(1/r);
L3 = sol*C3'*(1/r);
L4 = sol*C4'*(1/r);

 
OA1 = A - L1*C1; 
OB1 = [B L1]; 
OD1  = [0 0 ;0 0 ;0 0 ;0 0 ;0 0 ;0 0 ];


OA3 = A - L3*C3; 
OB3 = [B L3];
OD3 = [0 0 0; 0 0 0; 0 0 0;0 0 0; 0 0 0; 0 0 0];  

OA4 = A - L4*C4; 
OB4 = [B L4]; 
OD4 = [0 0 0 0; 0 0 0 0; 0 0 0 0;0 0 0 0; 0 0 0 0; 0 0 0 0];


OA = OA1; 
OB = OB1;   
OC = eye(6);
OD = OD1; 

sstem = ss(OA,OB,OC,OD);
step(sstem);

%{
The code for G part, to plot the output feedback controller (non-linear system) for the smallest output vector x.
%}


K = lqr(A,B,Q,R);

Nn=0.1;
Qn=20;
Rn=0.01;
D1=[0];

system=ss(A,B,C1,D1);

%{
For reference tracking we use kalman estimator
%}
[kalman_est,L,~]=kalman(system,Qn,Rn,Nn);

Ac = [A-B*K B*K;zeros(size(A)) A-L*C1];
Bc = [B;zeros(size(B))];
Cc = [C1 zeros(size(C1))];
finn = ss(Ac,Bc,Cc,D1);



