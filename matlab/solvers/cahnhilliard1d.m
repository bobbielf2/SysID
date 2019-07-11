function [C,C0,t,x] = cahnhilliard1d(C0,th_val,th_ind)
v = 0; % verbose
th = [0.5,-1,0,1];
if nargin >= 2 && ~isempty(th_val)
    if nargin < 3 || isempty(th_ind)
        th_ind = 1:numel(th_val);
    elseif numel(th_val)~=numel(th_ind)
        error('th_val and th_ind need to have same length!')
    end
    th(th_ind) = th_val;
end
k = th(1); % 0.05
dt = 0.01; T = 50; n = round(T / dt); % time, m pts in [0,T]
dx = 0.1; L = 60; m = 2*L / dx; % space, n pts in [-L,L]
t = 0:dt:T;
x = -L:dx:L;

[D2,D4] = DMAT(k,dt,dx,m);

%%

C = zeros(m+1,n+1);
if nargin > 0 && ~isempty(C0)
    C(:,1) = C0;
else
    C(:,1) = (rand(m+1,1)-0.5)*0.1; % initial condition
    C(:,1) = C(:,1) - mean(C(:,1)); % make sure zero mean
    %C(:,1) = C(:,1) + 0.6;
    C0 = C(:,1);
end


tic
for j = 1:n
    if v && mod(j,n/5)==0, disp([num2str(j/n*100),'% done.']), end
    c0 = C(:,j);
    % newton solve
    c1 = newton(@(y) F(y,c0,D2,D4,th), @(y) JF(y,c0,D2,D4,th), c0);
    C(:,j+1) = c1;
end
if v
    toc
end

if v % visualize
    for j = 1:round(n/200):n+1
        plot(x,C(:,j))
        ylim([-1,1]*2)
        xlim([-1,1]*L)
        title(sprintf('T = %.1f',(j-1)*dt))
        drawnow
    end
end

end


function [D2,D4] = DMAT(k,dt,dx,n)


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

% % 1st derivative matrix, modified by BCs
% a1 = dt / (2*dx);
% D1 = toeplitz([0, -1, zeros(1, n-1)] * a1 ,[0, 1, zeros(1, n-1)] * a1);
% D1(1,2) = 0; D1(end,end-1) = 0;
% D1 = sparse(D1);
end

function out = F(c1,c0,D2,D4,th)
out = c1-c0 + D4 * ((c0+c1)/2) - D2*(...
                                    th(2)/2*(c1+c0)...
                                  + th(3)/3*(c1.^2+c1.*c0+c0.^2)...
                                  + th(4)/4*(c1.^3+c1.^2.*c0+c1.*c0.^2+c0.^3));
end

function J = JF(c1,c0,D2,D4,th)
nn = size(D2,1);
J = speye(nn) + D4/2 - D2 * spdiags(th(2)/2 + th(3)/3*(2*c1+c0) + th(4)/4*(3*c1.^2+2*c1.*c0+c0.^2),0,nn,nn);
end


function y1 = newton(F, JF, y0)
dy = JF(y0)\F(y0);
y1 = y0 - dy;
res0 = norm(dy); % initial residual
res = res0;
relres = 1; % relative residual
count = 1;
while (relres > 1e-9 && res > 1e-11) && count <= 10
    count = count + 1;
    y0 = y1;
    dy = JF(y0)\F(y0);
    y1 = y0 - dy;
    res = norm(dy); % abs residual
    relres = res/res0; % rel residual
    %if count > 10, fprintf('iter=%2d, relres=%.2e, res=%.2e\n',count,relres,res), end
end
end