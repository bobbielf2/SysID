k = 0.5;
dt = 0.001*1; T = 1000; m = T / dt; % time, m pts in [0,T]
dx = 0.1; L = 5; n = 2*L / dx; % space, n pts in [-L,L]
t = 0:dt:T;
x = -L:dx:L;

[D1,D2,D4] = DMAT(k,dt,dx,n);
A = speye(n-1) + D4; condest(A)

%%
start_fresh = 0;
if start_fresh
    C = zeros(n+1,m+1);
    C(:,1) = (rand(n+1,1)-0.5)*0.1; % initial condition
else
    C(:,1) = C(:,end); % initial condition
    C(:,m+2:end) = [];
end
tic
for j = 1:m
    %j
    c0 = C(2:end-1,j);
    if 1
        c1 = A\RHSimp(c0,D1,D2); % implicit scheme
    else
        c1 = c0 + RHSexp(c0,D1,D2,D4); % explicit scheme
    end
    C(:,j+1) = [(15*c1(1)-6*c1(2)+c1(3))/10;
                c1;
                (15*c1(end)-6*c1(end-1)+c1(end-2))/10];
end
tt = toc
tt/m

for j = 1:m/100:m+1
    plot(x,C(:,j))
    ylim([-1,1]*2)
    xlim([-1,1]*L)
    title(sprintf('T = %.1f',(j-1)*dt))
    drawnow
end

%%
function [D1,D2,D4] = DMAT(k,dt,dx,n)


% 4th derivative matrix, modified by BCs
a4 = k * dt / dx^4;
D4 = toeplitz([6, -4, 1, zeros(1, n-4)] * a4);
D4(1:2,1:3) = [1, -8/5, 3/5; -5/2, 27/5, -39/10] * a4;
D4(end-1:end,end-2:end) = rot90(D4(1:2,1:3),2);
D4 = sparse(D4);

% 2nd derivative matrix, modified by BCs
a2 = dt / dx^2;
D2 = toeplitz([-2, 1, zeros(1, n-3)] * a2);
D2(1,1:3) = [-1/2, 2/5, 1/10] * a2;
D2(end,end-2:end) = fliplr(D2(1,1:3));
D2 = sparse(D2);

% 1st derivative matrix, modified by BCs
a1 = dt / (2*dx);
D1 = toeplitz([0, -1, zeros(1, n-3)] * a1 ,[0, 1, zeros(1, n-3)] * a1);
D1(1,1:3) = [3/2, -8/5, 1/10] * a1;
D1(end,end-2:end) = fliplr(D1(1,1:3));
D1 = sparse(D1);
end

function r = RHSexp(c,D1,D2,D4)
r = 6*c .* (D1*c).^2 + (3*c.^2 - 1) .* (D2*c) - D4*c;
end

function r = RHSimp(c,D1,D2)
r = c + 6*c .* (D1*c).^2 + (3*c.^2 - 1) .* (D2*c);
end