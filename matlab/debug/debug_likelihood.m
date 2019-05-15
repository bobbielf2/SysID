% System identification using Bayesian inference

% MCMC initialize

test_case       = 7;    % test case num
niter           = 1000; % num of iterations
noise_level     = 0.01; % noise (used as std of pointwise gaussian noise)
sparse_prior    = 0;    % use Laplace prior?   
likelihood_type = 0;    % type of likelihood? 0=multivar Gaussian, 1=square of 2-norm of err (Gamma), 2=mean of err (Normal)

% generate test case
model = buildTestCase(test_case,niter,noise_level,sparse_prior,likelihood_type);

% (relative) posterior(=likelihood*prior) on a regular theta-grid
n = 10;
thth = linspace(-1,1,n+1)*100;
p = zeros(n,n);
for i = 1:n+1
    disp(i)
    for j = 1:n+1
        for k = n/2+1
            th = model.mu_th + [thth(i)/40,thth(j)/10,thth(k)/50];
            p(i,j) = posterior(th,model);
        end
    end
end
%%
subplot(1,2,1)
[U, V, ~, x] = reactDiffuse1d2sp([model.mu_th(1)/10,model.mu_th(2),0.1, model.mu_th(3), 1, 0.9, -1],model.UV0);
plot(x,[U(:,end),V(:,end)])
legend({'species 1','species 2'})
xlabel('space')
ylabel('concentration')
title('2-species reaction-diffusion (snapshot)')
subplot(1,2,2)
[x,y] = meshgrid(model.mu_th(1)+thth/40,model.mu_th(2)+thth/10);
% surf(x,y,log10(p)); colorbar
surf(x,y,(p)); colorbar
title('$\log_{10}p(y|\theta)p(\theta)$, ($\theta_1$-$\theta_2$ slice at $\theta_3^*$)','interpreter','latex')