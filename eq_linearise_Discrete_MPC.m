clear all
clc

syms z1 z2 z3 u1 u2 u3 u4 d1 d2
m = 1500; g = 9.81; rw = 0.28; x1 = 1.3; x2 = 1.3;
x3 = -1.7; x4 = -1.7; y1 = 0.8; y2 = -0.8; y3 = 0.8;
y4 = -0.8; Izz = 2000; d3 = 0; d4 = 0;

% calculating reaction forces
Fz12 = 0.5* m * g *( abs ( x3 ) /( x1 + abs ( x3 ) ) ) ;
Fz34 = 0.5* m * g *( x1 /( x1 + abs ( x3 ) ) ) ;
% calculating longitudinal forces
F1t = u1 / rw ; F2t = u2 / rw ; F3t = u3 / rw ; F4t = u4 / rw ;

alpha_p = 6* pi /180; % peak alpha value
mu_p = 0.9; % peak mu value


% calculating lateral forces
alpha1s = d1 - atan (( z2 + z3 * x1 ) /( z1 - z3 * y1 ) ) ;
alpha2s = d2 - atan (( z2 + z3 * x2 ) /( z1 - z3 * y2 ) ) ;
alpha3s = 0 - atan (( z2 + z3 * x3 ) /( z1 - z3 * y3 ) ) ;
alpha4s = 0 - atan (( z2 + z3 * x4 ) /( z1 - z3 * y4 ) ) ;

F1s = Fz12 * ((2* alpha_p * mu_p ) /( alpha_p ^2 + alpha1s ^2) ) * alpha1s ; % tyre ( d1 , x1 , y1 , z1 , z2 , z3 ) ;
F2s = Fz12 * ((2* alpha_p * mu_p ) /( alpha_p ^2 + alpha2s ^2) ) * alpha2s ;
F3s = Fz34 * ((2* alpha_p * mu_p ) /( alpha_p ^2 + alpha3s ^2) ) * alpha3s ;
F4s = Fz34 * ((2* alpha_p * mu_p ) /( alpha_p ^2 + alpha4s ^2) ) * alpha4s ;

fz1 = (1/m)*( F1t*cos(d1) - F1s*sin(d1) + F2t*cos(d2)- F2s *sin (d2) + F3t + F4t + m*z2*z3);
fz2 = (1/m)*( F1t*sin(d1) + F1s*cos(d1) + F2t*sin(d2)+ F2s*cos(d2) + F3s + F4s - m*z1*z3);
fz3 = (1/ Izz)*( x1 *( F1t*sin(d1) + F1s*cos(d1)) - y1 *( F1t *cos (d1) - F1s *sin (d1)) + x2 *( F2t*sin(d2)...
+ F2s*cos(d2)) - y2 *( F2t*cos(d2) - F2s* sin(d2))+ x3*F3s - y3*F3t + x4*F4s - y4*F4t);

Ac = jacobian ([ fz1; fz2; fz3], [z1 z2 z3 ]);
Bc = jacobian ([ fz1; fz2; fz3], [u1 u2 u3 u4 d1 d2]);
Cc = diag ( ones (3 ,1));


syms d
w = (y1 + abs(y2)); % axle track
l = (x2 + abs(x3)); % wheelbase
d1 = acot (cot (d) - w /(2* l));
d2 = acot (cot (d) + w /(2* l));

% substitute one steering angle input
% B(: ,5) = B(: ,5) + B(: ,6); B(: ,6) = []; % merge d1 ,d2
% A = eval (A); B = eval (B);

% feed equilibrium point in
% find an equilibrium point
assume (z1 > 0)
u1 = 0; u2 = 0; u3 = 0; u4 = 0; d=0; d1=0; d2=0;
f1 = eval (fz1 ); f2 = eval (fz2 ); f3 = eval (fz3);
[r1 , r2 , r3] = solve ([0 == f1 , 0 == f2 , 0 == f3],[z1 , z2 , z3 ]);

z1 = r1; z2 = r2; z3 = r3; u1 = 0; u2 = 0; u3 = 0; u4 = 0; d = 0;
Ac = eval (Ac); 
Bc = eval (Bc);

Bc(: ,5) = Bc(: ,5) + Bc(: ,6); Bc(: ,6) = []; % merge d1 ,d2
Ts = 1.e-4;
[A, B]=conti2disc(Ac,Bc,Ts)
A=double(A)
B=double(B)
C=eye(3);
D=0*eye(3,5);
%%%%%%%%%%%%Also try with linmod

% clearvars -except A B m rw
B(: ,5) = [];
% computing predictive control matrices
n = 4; % horizon parameter
% P array
[rA , cA] = size (A);
P = zeros (rA*n, cA); % pre - allocating
subP = A;
for i = 1:n
for j = 1: rA
P(j+(i -1)*rA ,:) = subP (j ,:);
end
subP = subP *A; % A^n
end
% H matrix
[rB , cB] = size (B);
H = zeros (rB*n, cB*n); % pre - allocating ;
subH = B;
for l = 1:n
for i = (1:n+1-l)
for j = 1: rB
for k = 1: cB
H(j+(i -1)*rB +(l -1)*rB ,k+(i -1)*cB) = subH (j,k);
end
end
end
subH = A* subH ;
end


% creating system control structure
% clearvars -except P H Ts m rw n
rw=0.28
Ts=1.e-4
l = 3; % wheelbase
% sys = struct ('P',P,'H',H,'Ts',Ts ,'m',m,'rw',rw,'l',l);
% sysBus = Simulink .Bus . createObject (sys);

% defining Q and R weighting matrices
r = ones (4*n ,1); R = diag (r);
q = 5e6* ones (3*n ,1); Q = diag (q);

% creating control structure
% wgt = struct ('R',R,'Q',Q);
% busWgt = Simulink .Bus . createObject (wgt );