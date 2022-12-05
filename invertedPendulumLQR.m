clear; 
clc;
close all;

disp('Inverted Pendulum on a Cart')
% Part A: Equations of motion

syms m l M g t F

syms x(t) theta(t) 

dx = diff(x);
dth = diff(theta);

%---Solve for Kinetic Energy
T = 0.5 * M * dx^2 + 0.5*m* (dx - l*dth*cos(theta)).^2 + 0.5* m * (l*dth^2 * sin(theta)).^2;

%---Solve for Potential Energy
V = m * g * l * cos(theta);

%--Solve for Langrangian of the System
L = T - V;

%%%%%%%%%%%%%%%%%%%%%%%

%---Euler-Lagrange Method for Equations of Motion

dL_x = diff(L, x);
dL_dx = diff(L, dx);
term_x = diff(dL_dx) - dL_x - F

dL_theta = diff(L, theta);
dL_dtheta = diff(L, dth);
term_theta = diff(dL_dtheta) - dL_theta

ddx = diff(dx);
ddth = diff(dth);

%-------Solve for DDx, DDtheta----------
%---Make Equations Symbolic for solver---
syms Dx Dth DDx DDtheta x_a th
term_x = subs (term_x, {x, theta, dx, dth ddx, ddth},...
    {x_a, th, Dx, Dth, DDx, DDtheta});

term_theta = subs (term_theta, {x, theta, dx, dth ddx, ddth},...
    {x_a, th, Dx, Dth, DDx, DDtheta});

eqns = [term_x; term_theta];

vars = [DDx; DDtheta];

S = solve(eqns, vars);
disp('Equations of Motion')
DDx = S.DDx    %x DD = doubledot
DDtheta = S.DDtheta

%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Solve for Matrix A and B using Jacobian---
X = [x_a; Dx; th; Dth];
Xdot = [Dx; DDx; Dth; DDtheta];

U = F;  %Input

J_U = jacobian(Xdot, U)
J_X = jacobian(Xdot, X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------Part B: Plug in Equilibrium Points
% x_a, th Dx, Dth
disp('Linearized System about the equilibrium Point')
A = subs (J_X, {sym('x_a'),sym('th'), sym('Dx'),sym('Dth')}, ...
    {0,0,0,0})

B = subs (J_U, {sym('x_a'),sym('th'), sym('Dx'),sym('Dth')}, ...
    {0,0,0,0})

% Create Controllability Matrix/conditions for M, m l

% C = [C = B AB A2B A3B ]

C = B;
temp = B;

% size(A) %(4, 4)
% rank(A) %3

n = size(A)-1;

for i = 1:n
%     temp = (A.^i)*temp;
    temp = A*temp;
    C = horzcat(C, temp);
end

disp('Controllability Matrix (C) is ')
C

% size(C) %(4,4)
C_det = simplify (det(C))

%non-singular if M and l are non-zero

% Plug in M, m, l, g values into C matrix
C = subs(C, {M, m, l, g},{1000, 100, 10, 10})

%Check for controllability
C_det = sprintf('%10e',det(C)) %'-1.000000e-14'

if rank(C) == size(A,1)
    disp('System is controllable!')
else
    disp('System is NOT controllable!!!')
end

%Obtain LQR controller
syms F X_1 X_2 X_3 X_1d
tspan = 0:0.1:200;
% 10 degrees = pi/18
s0 = [0;  0; pi/18; 0]; %1 meter push to serve as step input
X = [X_1; X_2; X_3; X_1d;];

A = subs(A, {M, m, l, g},{1000, 100, 10, 10});
B = subs(B, {M, m, l, g},{1000, 100, 10, 10});

%Convert sym matrixes to num
A = double (A);
B = double (B);

R = .00001;
Q = diag([10 1000 1000 10000]);

[K,S,e] = lqr(A,B,Q,R);
disp('LQR Controller control Matrix (K) is ')
K
%K = [-0.0100   -0.1212    2.4575    2.4658]

[t,state_history] = ode45(@(t,state)linearized_model(state,t,A,B,K),tspan,s0);

figure('Name', 'Linearized Model with Force');
subplot(2,2,1);
plot(t,state_history(:,1),'k')
grid on;
ylabel('x position of the cart')
xlabel('time in s')

subplot(2,2,2);
plot(t,state_history(:,2),'r')
grid on;

subplot(2,2,3);
plot(t,state_history(:,3),'b')
grid on;

%----------------------------
% Original Non-linear system
[t,state_history] = ode45(@(t,state)og_model(state,t,K),tspan,s0);

figure('Name', 'Non-Linearized Model with Force');
subplot(2,2,1);
plot(t,state_history(:,1),'k')
grid on;
ylabel('x position of the cart')
xlabel('time in s')

subplot(2,2,2);
plot(t,state_history(:,2),'r')
grid on;

subplot(2,2,3);
plot(t,state_history(:,3),'b')
grid on;

function sdot = linearized_model(s,t,A,B,K)
%% Linearized Model
sdot = (A-B*K)*s;
end


function sdot = og_model(s,t,K)
%% Non-Linearized Model
F = -K*s;
sdot = [s(4); 
    (F + 500.0*sin(2.0*s(2)) + 500.0*sin(2.0*s(3)))/(100*sin(s(2))^2 + 100*sin(s(3))^2 + 1000);
    -(10*(100*sin(s(3)))^2 + 1100)*sin(s(2)) + (F)*cos(s(2)) - 250.0*sin(s(2) - 2.0*s(3)) + 250.0*sin(s(2) + 2.0*s(3))/(20*(100*sin(s(2))^2 + 100*sin(s(3))^2 + 1000));
    -(10*(100*sin(s(2))^2 + 1100)*sin(s(3)) + (F)*cos(s(3)) + 250.0*sin(2.0*s(2) - s(3)) + 250.0*sin(2.0*s(2) + s(3)))/(10*(100*sin(s(2))^2 + 100*sin(s(3))^2 + 1000));     
];
end 


function sdot = no_input_lin_model(s,t,A)
%% Model without forces acting on it
sdot = (A)*s;
end
