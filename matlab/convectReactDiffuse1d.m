function U = convectReactDiffuse1d(th)

D = th(2); % diffisivity
a = th(1); % convection velocity
L = 10; m = 100; % domain = [0,L]
T = 0.85; n = 100; % time = [0,T]
dx = L*2/m;
dt = T/n;

U = zeros(m+1,n+1); % store solutions

% form iteration matrix
al = a*dt/dx^2;
be = D*dt/dx^2;
c = [1-al-2*be, al+be, zeros(1,m-3)];
r = [1-al-2*be, be, zeros(1,m-3)];
A = toeplitz(c,r);
f = zeros(m-1, 1);

x = linspace(0,L,m+1)';
U(:,1) = exp(-x); % initial condition

for i = 1:n
    f(1) = (al+be)*U(1,i); f(end) = be*U(end,i);
    R = react(U(2:end-1,i),th);
    U(2:end-1,i+1) = A*U(2:end-1,i) + f + dt * R;
    U(1,i+1) = exp((i+1)*dt);
    U(end,i+1) = U(end-1,i+1)/(1+dx);
end

% surf(U,'edgealpha',0)



function R = react(u,th)
% reaction term
R = th(3) * u + th(4) * u.^2;