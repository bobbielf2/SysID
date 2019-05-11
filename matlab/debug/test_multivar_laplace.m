s = linspace(-2,2,101);
[x, y] = meshgrid(s);
p = zeros(size(x));
for i = 1:numel(x)
    p(i) = mvlaplaceSym([x(i); y(i)]);
end

figure
subplot(1,2,1)
surf(p,'edgealpha',0.1)
subplot(1,2,2)
z = exp(-(abs(x).^1.0+abs(y).^1.0));
surfc(x,y,z,'edgealpha',0.1)
%contour(x,y,z)


function p = mvlaplaceSym(x)

n = length(x);
nu = (2-n)/2;
SIGMA = eye(n); % diag of covariance mtx

xSIGx = dot(x, x ./ diag(SIGMA));
lambda = det(SIGMA);
p = 2/sqrt(2*pi)^n/sqrt(lambda) * besselk(nu,sqrt(2*xSIGx)) * sqrt(xSIGx/2).^nu ;

end