function U = reactDiffuse1d(th)
% Reaction-Diffusion equation 1-D
% u_t = D*u_xx + R(u)

D = 1; % diffisivity
L = 5; m = 100; % domain = [-L,L]
T = 5; n = 30; % time = [0,T]
dx = L*2/m;
dt = T/n;

U = zeros(m+1,n+1); % store solutions

% form iteration matrix
a = D*dt/dx^2;
r = [1+2*a, -a, zeros(1,m-1)];
A = toeplitz(r);
A(1,2) = -2*a; A(m+1,m) = -2*a;

x = linspace(-L,L,m+1)';
U(:,1) = 0.05*exp(-5*x.^2); % initial condition

for i = 1:n
    R = react(U(:,i),th);
    U(:,i+1) = A\(U(:,i) + dt*R);
end

function R = react(u,th)
% reaction term
R = th(1) * u - th(2) * u.^2 + th(3) * u.^3 + 0 * u.^4;