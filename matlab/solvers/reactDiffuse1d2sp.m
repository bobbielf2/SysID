function [U, V, UV0, t, x] = reactDiffuse1d2sp(UV0,th_val,th_ind)
% Reaction-Diffusion equation 1-D, 2 species.
% u_t = D1*u_xx + R1(u,v)
% v_t = D2*v_xx + R2(u,v)

th = [0.1, 4, 0.1, -1, 1, 0.9, -1]; % default param
if nargin >= 2 && ~isempty(th_val)
    if nargin < 3 || isempty(th_ind)
        th_ind = 1:numel(th_val);
    elseif numel(th_val)~=numel(th_ind)
        error('th_val and th_ind need to have same length!')
    end
    th(th_ind) = th_val;
end

D1 = th(1); % diffisivity
D2 = th(2); % diffisivity
L = 40; m = ceil(L/0.1); % domain = [-L,L]
dx = L*2/m;
T = 15; n = ceil(T/0.01)*4; % time = [0,T]
dt = T/n;

% store solutions
U = zeros(m+1,n+1); % 1st component
V = zeros(m+1,n+1); % 2nd component

% form iteration matrix
al1 = D1*dt/dx^2;
c = [1-2*al1, al1, zeros(1,m-3)];
A1 = toeplitz(c);

al2 = D2*dt/dx^2;
c = [1-2*al2, al2, zeros(1,m-3)];
A2 = toeplitz(c);

% vectors for boundary condition
f = zeros(m-1, 1);

x = linspace(-L,L,m+1)';
t = linspace(0,T,n+1)';
if nargin > 1 && ~isempty(UV0)
    U(:,1) = UV0(:,1);
    V(:,1) = UV0(:,2);
else
    %rng(1)
    U(:,1) = 0.5 + 0.2 * (rand(size(x))-0.5); % initial condition
    V(:,1) = 0.5 + 0.2 * (rand(size(x))-0.5); % initial condition
    UV0 = [U(:,1), V(:,1)];
end

for i = 1:n
    [R1, R2] = react(U(2:end-1,i),V(2:end-1,i),th);
    f(1) = al1*U(1,i); f(end) = al1*U(end,i);
    U(2:end-1,i+1) = A1*U(2:end-1,i) + f + dt * R1;
    U(1,i+1) = U(2,i+1); U(end,i+1) = U(end-1,i+1);
    f(1) = al2*V(1,i); f(end) = al2*V(end,i);
    V(2:end-1,i+1) = A2*V(2:end-1,i) + f + dt * R2;
    V(1,i+1) = V(2,i+1); V(end,i+1) = V(end-1,i+1);
end


% concaternate for system identification
if nargout == 1
    U = [U(1:5:end,1:100:end); V(1:5:end,1:100:end)];
    %U = [U; V];
end

% for i = 1:round(n/300):n+1
%     plot(x, U(:,i)); hold on
%     plot(x, V(:,i)); hold off
%     legend({'u','v'})
%     ylim([0,4])
%     title(['time = ',num2str(dt*(i-1))])
%     drawnow
% end



function [R1, R2] = react(u,v,th)
% reaction term
R1 = th(3) + th(4) * u + th(5) * u.^2 .* v;
R2 = th(6) + th(7) * u.^2 .* v;
