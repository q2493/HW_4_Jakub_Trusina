% Program for plotting displacement, velocity and acceleration of 
% connecting rod and piston of a Short Crank Mechanism
clc, clear

a = 0.1; b = 0.2;           % parameters
w = 1;                      % angular velocity

T = 10;                     % Time end
points = 101;               % number of division on time interval
t = linspace(0,T,points);
phi = pi/6 + w*t;            % angle 

J = @(x) [-b*sin(x(1))  -1 ;
          -b*cos(x(1))    0];   % Jacobian
tol = 1e-5;                     % tolerance
guess = [0 ; 0];              % Initial Guess

% Inicialization of vectors
x=zeros(length(t),2);   % displacement
v=zeros(length(t),2);   % velocity
aa=zeros(length(t),2);  % acceleration

% Utilization of Newton method for solving system of equations
for ii = 1:length(t)
    F = @(u) [a*cos(phi(ii)) + b*cos(u(1)) - u(2);
             a*sin(phi(ii)) - b*sin(u(1))        ];
    [x(ii,:), n] = NR_method(F, J, guess, tol);

    Fv=@(u)[-a*w*sin(phi(ii))-b*u(1)*sin(x(ii,1))-u(2);
             a*w*cos(phi(ii))-b*u(1)*cos(x(ii,1))];
    [v(ii,:), nv] = NR_method(Fv, J, guess, tol);    

    Fa=@(u)[-a*w^2*cos(phi(ii)) - b*u(1)*sin(x(ii,1)) - b*v(ii,1)^2*cos(x(ii,1)) - u(2);
            -a*w^2*sin(phi(ii)) - b*u(1)*cos(x(ii,1)) + b*v(ii,1)^2*sin(x(ii,1))];
    [aa(ii,:), na] = NR_method(Fa, J, guess, tol);  
end

theta = x(:,1);     % angle displacement
d = x(:,2);         % slider displacement
theta_dot = v(:,1);     % angular velocity
d_dot = v(:,2);         % slider velocity
theta_ddot = aa(:,1);   % angular acceleration
d_ddot = aa(:,2);       % slider acceleration
    
% Plotting
subplot(2,1,1)
plot(t,theta, 'k-', t, theta_dot, 'b-',t, theta_ddot, 'r-','Linewidth', 1.2);
lgnd1 = legend('$\theta$', '$\dot{\theta}$', '$\ddot{\theta}$', 'Location','northeast');
set(lgnd1, 'Interpreter', 'latex');
xlabel('t [s]');
ylabel('Quantity [Base unit]');
title('Conecting Rod Position vs Time');
grid on

subplot(2,1,2)
plot(t,d, 'k-', t, d_dot, 'b-',t, d_ddot, 'r-','Linewidth', 1.2);
lgnd2 = legend('x', 'v', 'a', 'Location','northeast');
set(lgnd2, 'Interpreter', 'latex');
xlabel('t [s]');
ylabel('Quantity [Base unit]');
title('Piston Position vs Time');
grid on

