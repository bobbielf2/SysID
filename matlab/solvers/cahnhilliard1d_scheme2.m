k = 0.5;
dt = 0.001*1000; T = 1000; m = round(T / dt); % time, m pts in [0,T]
dx = 0.1; L = 5; n = 2*L / dx; % space, n pts in [-L,L]
t = 0:dt:T;
x = -L:dx:L;

[D1,D2,D4] = DMAT(k,dt,dx,n);
A = speye(n+1) + D4; condest(A)

%%
start_fresh = 0;
if start_fresh
    C = zeros(n+1,m+1);
    C(:,1) = (rand(n+1,1)-0.5)*0.1; % initial condition
    C(:,1) = C(:,1) - mean(C(:,1)); % make sure zero mean
else
    C(:,1) = C(:,end); % initial condition
    C(:,1) = C(:,1) - mean(C(:,1)); % make sure zero mean
    C(:,m+2:end) = [];
end

for j = 1:m
    %j
    c0 = C(:,j);
    if 1
        c1 = A\RHSimp(c0,D1,D2); % implicit scheme
    else
        c1 = c0 + RHSexp(c0,D1,D2,D4); % explicit scheme
    end
    C(:,j+1) = c1;
end

for j = 1:m/100:m+1
    plot(x,C(:,j))
    ylim([-1,1]*2)
    xlim([-1,1]*L)
    title(sprintf('T = %.1f',(j-1)*dt))
    drawnow
end

disp(norm(C(:,end) - tanh(x.')))

%%
function [D1,D2,D4] = DMAT(k,dt,dx,n)


% 4th derivative matrix, modified by BCs
a4 = k * dt / dx^4;
D4 = toeplitz([6, -4, 1, zeros(1, n-2)] * a4);
D4(1,2:3) = [-8, 2] * a4; D4(end,end-2:end-1) = [2, -8] * a4;
D4(2,2) = 7 * a4; D4(end-1,end-1) = 7 * a4;
D4 = sparse(D4);

% 2nd derivative matrix, modified by BCs
a2 = dt / dx^2;
D2 = toeplitz([-2, 1, zeros(1, n-1)] * a2);
D2(1,2) = 2 * a2; D2(end,end-1) = 2 * a2;
D2 = sparse(D2);

% 1st derivative matrix, modified by BCs
a1 = dt / (2*dx);
D1 = toeplitz([0, -1, zeros(1, n-1)] * a1 ,[0, 1, zeros(1, n-1)] * a1);
D1(1,2) = 0; D1(end,end-1) = 0;
D1 = sparse(D1);
end

function r = RHSexp(c,D1,D2,D4)
r = 6*c .* (D1*c).^2 + (3*c.^2 - 1) .* (D2*c) - (D4*c);
end

function r = RHSimp(c,D1,D2)
r = c + 6*c .* (D1*c).^2 + (3*c.^2 - 1) .* (D2*c);
end