% System identification using Bayesian inference

% MCMC initialize

test_case       = 8;    % test case num
niter           = 1000; % num of iterations
noise_level     = 0.01; % noise (used as std of pointwise gaussian noise)
sparse_prior    = 0;    % use Laplace prior?   
likelihood_type = 0;    % type of likelihood? 0=multivar Gaussian, 1=square of 2-norm of err (Gamma), 2=mean of err (Normal)
proposal_anneal = 0;    % use annealing for the proposal?
fixInit         = 1;    % remember initial condition?

% generate test case
model = buildTestCase(test_case,niter,noise_level,sparse_prior,likelihood_type,proposal_anneal,fixInit);

% (relative) posterior(=likelihood*prior) on a regular theta-grid
th_true = [log(0.1), log(4), -1];
n = 50;
thth = linspace(-1,1,n+1).';
thth = [thth*5, thth*5, thth];
p = zeros(n,n);
p_eps = zeros(n,n);
for i = 1:n+1
    disp(i)
    for j = 1:n+1
        for k = n/2+1
            th = th_true + [thth(i,1),thth(j,2),thth(k,3)];
            [p(i,j),~,~,p_eps(i,j)] = model.posterior(th);
        end
    end
end
%%
figure(1)
subplot(1,2,1)
[U, V, ~, x] = reactDiffuse1d2sp([exp(th_true(1)),exp(th_true(2)),0.1, th_true(3), 1, 0.9, -1],model.UV0);
plot(x,[U(:,end),V(:,end)])
legend({'species 1','species 2'})
xlabel('space')
ylabel('concentration')
title('2-species reaction-diffusion (snapshot)')
subplot(1,2,2)
[x,y] = meshgrid(th_true(1)+thth(:,1),th_true(2)+thth(:,2));
% surf(x,y,-log10(-p)); colorbar
surfc(x,y,(p)); colorbar
title('$\log_{10}p(y|\theta)p(\theta)$, ($\theta_1$-$\theta_2$ slice at $\theta_3^*$)','interpreter','latex')
figure(2)
subplot(1,2,1)
surf(x,y,p_eps); colorbar
set(gca,'SortMethod','ChildOrder')
title('$\log_{10}p(y|\theta)$, ($\theta_1$-$\theta_2$ slice at $\theta_3^*$)','interpreter','latex')
subplot(1,2,2)
surf(x,y,(p)); colorbar
set(gca,'SortMethod','ChildOrder')
title('$\log_{10}p(y|\theta)p(\theta)$, ($\theta_1$-$\theta_2$ slice at $\theta_3^*$)','interpreter','latex')
