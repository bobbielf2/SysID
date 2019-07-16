k = 0.5;
dt = 0.01; T = 500; m = round(T / dt); % time, m pts in [0,T]
dx = 0.1; L = 5; n = 2*L / dx; % space, n pts in [-L,L]
t = 0:dt:T;
x = -L:dx:L;

[D1,D2,D4] = DMAT(k,dt,dx,n);

%%
start_fresh = 0;
if start_fresh
    C = zeros(n+1,m+1);
    C(:,1) = (rand(n+1,1)-0.5)*0.1; % initial condition
    C(:,1) = C(:,1) - mean(C(:,1)); % make sure zero mean
else
    C(:,1) = C(:,end); % initial condition
%      C(:,1) = C(:,1) - mean(C(:,1)); % make sure zero mean
    C(:,m+2:end) = [];
end

options = optimset('Display','off');
for j = 1:m
    if mod(j,m/50)==0, disp([num2str(j/m*100),'% done.']), end
    c0 = C(:,j);
    % newton
    y0 = c0;
    count = 1;
    while 1
        dy = JF(y0,c0,D2,D4)\F(y0,c0,D2,D4);
        err = norm(dy);
        y1 = y0 - dy;
        if count > 10, fprintf('iter=%2d, res=%.2e\n',count,err), end
        if err < 1e-12
            break
        end
        y0 = y1;
        count = count + 1;
    end
    C(:,j+1) = y1;
%     % matlab fsolve
%     fun = @(y) F(y,c0,D2,D4);
%     [c1,fval,exitflag] = fsolve(fun,c0,options);
%     if exitflag ~= 1, disp('flag!'), end
%     C(:,j+1) = c1;
end


for j = 1:m/200:m+1
    plot(x,C(:,j))
    ylim([-1,1]*2)
    xlim([-1,1]*L)
    title(sprintf('T = %.1f',(j-1)*dt))
    drawnow
end

disp(sum(C(:,end)*sign(C(end,end)) - tanh(x.'))*dx)

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

function out = F(c1,c0,D2,D4)
out = c1-c0 + D2*((c1+c0)/2 - (c1.^3+c1.^2.*c0+c1.*c0.^2+c0.^3)/4) + D4 * ((c0+c1)/2);
end

function J = JF(c1,c0,D2,D4)
nn = size(D2,1);
J = speye(nn) + (D2 + D4)/2 - D2 * spdiags(3*c1.^2 + 2*c1.*c0+c0.^2,0,nn,nn)/4;
end


function r = RHSexp(c,D1,D2,D4)
r = 6*c .* (D1*c).^2 + (3*c.^2 - 1) .* (D2*c) - (D4*c);
end

function r = RHSimp(c,D1,D2)
r = c + 6*c .* (D1*c).^2 + (3*c.^2 - 1) .* (D2*c);
end