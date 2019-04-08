s = linspace(-2,2,100);
[x, y] = meshgrid(s);
p = zeros(size(x));
for i = 1:numel(x)
    p(i) = mvlaplaceSym([x(i); y(i)]);
end

surf(p)


function p = mvlaplaceSym(x)

n = length(x);
nu = (2-n)/2;
SIGMA = eye(n); % diag of covariance mtx

xSIGx = dot(x, x ./ diag(SIGMA));
lambda = det(SIGMA);
p = 2/sqrt(2*pi)^n/sqrt(lambda) * besselk(nu,sqrt(2*xSIGx)) * sqrt(xSIGx/2).^nu ;

end