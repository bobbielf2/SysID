% System identification using Bayesian inference
setup

% MCMC initialize

test_case       = 8;    % test case num
niter           = 500; % num of iterations
noise_level     = 0.01; % noise (used as std of pointwise gaussian noise)
sparse_prior    = 0;    % use Laplace prior?   
likelihood_type = 0;    % type of likelihood? 0=multivar Gaussian, 1=square of 2-norm of err (Gamma), 2=mean of err (Normal)
proposal_anneal = 0;    % use annealing for the proposal?
fixInit         = 0;    % remember initial condition?

% generate test case
model = buildTestCase(test_case,niter,noise_level,sparse_prior,likelihood_type,proposal_anneal,fixInit);

% Metropolis-Hastings Iteration
model = metropolis_hastings(model);
%model = metropolis_hastings(model,1000); % continue for 1000 more iters

%% plot resulted posterior distr
th      = cell2mat(model.th); % all the theta tried
p       = model.p; p(isnan(p)) = 0; % posterior prob corresp to th
i_burn  = floor(numel(p)/2):numel(p); % burn-in
th_T    = cell2mat(model.th_T); % Markov chain of theta's

switch test_case
    case 1
        subplot(1,2,1),cla
        scatter(th,p,5,'r')
        %scatter(th(i_burn),p(i_burn),5,'r')
        hold on, fplot(@sin,[-10,10])
        fplot(@(x) 0*x+sin(2),[-10,10])
        % find modes
        [th_mode,i_sort] = sort(th); p_mode = p(i_sort);                    % sort p in ascending order of th
        [p_mode ,loc]    = findpeaks(p_mode); th_mode = th_mode(loc);       % find modes of th based on local max of p
        [p_mode,i_sort] = sort(p_mode,'desc'); th_mode = th_mode(i_sort);   % sort modes in descending order
        th_mode = th_mode(1:min(end,6)); % pick 6 modes max
        title(sprintf(['predicted modes (in desc order): ', repmat('%.2f, ', 1,length(th_mode))],th_mode))
        subplot(1,2,2),cla
        plot(th_T) % plot Markov chain
    case 2
        figure(1)
        subplot(2,2,3)
        scatter3(th(i_burn,1),th(i_burn,2),p(i_burn),5,'r')
        subplot(2,2,4)
        plot(th_T) % plot Markov chain
        th_m = mean(th_T(i_burn,:));
        [~, mu_eps] = model.posterior(th_m);
        title(['predicted mean: ',num2str(mean(th_T(i_burn,:)))])
        xlabel(['likelihood mean $\mu_\epsilon=$',num2str(mu_eps)],'interpreter','latex')
    case {3, 4, 5}
        figure(1)
        subplot(2,3,2)
        scatter3(th(i_burn,1),th(i_burn,2),p(i_burn),5,'r')
        title('$\theta_1$-$\theta_2$ distribution','interpreter','latex')
        subplot(2,3,5)
        scatter3(th(i_burn,2),th(i_burn,3),p(i_burn),5,'r')
        title('$\theta_2$-$\theta_3$ distribution','interpreter','latex')
        subplot(2,3,3)
        scatter3(th_T(i_burn,1),th_T(i_burn,2),th_T(i_burn,3),5,'r')
        axis equal
        title('uncertainty window')
        xlabel('\theta_1');ylabel('\theta_2');zlabel('\theta_3')
        subplot(2,3,6)
        plot(th_T) % plot Markov chain
        legend({'\theta_1','\theta_2','\theta_3'})
        th_m = mean(th_T(i_burn,:));
        title(sprintf('predicted mean: %.3f, %.3f, %.3f',th_m),'interpreter','latex')
        [~, mu_eps] = model.posterior(th_m);
        xlabel(['likelihood mean $\mu_\epsilon=$',num2str(mu_eps)],'interpreter','latex')
        %ylim([-1,1]*0.002-1); ylim([-1,1]*0.05+4); % case 7 debug
        
        % subplot(2,3,5)
        % U = reactDiffuse1d(th_m);
        % surf(U), title('predicted solution')
    case 6
        figure(1)
        subplot(2,3,2)
        scatter3(th(i_burn,1),th(i_burn,2),p(i_burn),5,'r')
        title('$\theta_1$-$\theta_2$ distribution','interpreter','latex')
        subplot(2,3,5)
        scatter3(th(i_burn,2),th(i_burn,3),p(i_burn),5,'r')
        title('$\theta_2$-$\theta_3$ distribution','interpreter','latex')
        subplot(2,3,3)
        scatter3(th_T(i_burn,1),th_T(i_burn,2),th_T(i_burn,3),5,'r')
        axis equal
        title('uncertainty window')
        xlabel('\theta_1');ylabel('\theta_2');zlabel('\theta_3')
        subplot(2,3,6)
        plot(th_T) % plot Markov chain
        legend({'\theta_1','\theta_2','\theta_3','\theta_4','\theta_5','\theta_6','\theta_7'})
        th_m = mean(th_T(i_burn,:));
        title(sprintf('predicted mean: %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f',th_m),'interpreter','latex')
        [~, mu_eps] = model.posterior(th_m);
        xlabel(['likelihood mean $\mu_\epsilon=$',num2str(mu_eps)],'interpreter','latex')
    case 7
        figure(1)
        subplot(2,3,2)
        scatter3(th(i_burn,1),th(i_burn,2),p(i_burn),5,'r')
        title('$\theta_1$-$\theta_2$ distribution (log)','interpreter','latex')
        subplot(2,3,5)
        scatter3(th(i_burn,2),th(i_burn,3),p(i_burn),5,'r')
        title('$\theta_2$-$\theta_3$ distribution (log)','interpreter','latex')
        subplot(2,3,3)
        scatter3(th_T(i_burn,1),th_T(i_burn,2),th_T(i_burn,3),5,'r')
        axis equal
        title('uncertainty window')
        xlabel('\theta_1');ylabel('\theta_2');zlabel('\theta_3')
        subplot(2,3,6)
        plot(th_T) % plot Markov chain
        legend({'\theta_1','\theta_2','\theta_3'})
        th_m = mean(th_T(i_burn,:));
        title(sprintf('predicted mean: %.3f, %.3f, %.3f',th_m),'interpreter','latex')
        [~, mu_eps] = model.posterior(th_m);
        xlabel(['likelihood mean $\mu_\epsilon=$',num2str(mu_eps)],'interpreter','latex')
        %ylim([-1,1]*0.002-1); ylim([-1,1]*0.05+4); % case 7 debug
        
        % subplot(2,3,5)
        % U = reactDiffuse1d(th_m);
        % surf(U), title('predicted solution')
    case 8
        figure(1)
        subplot(2,3,2)
        scatter3(th(i_burn,1),th(i_burn,2),p(i_burn),5,'r')
        title('$\theta_1$-$\theta_2$ distribution (log)','interpreter','latex')
        subplot(2,3,5)
        scatter3(th(i_burn,2),th(i_burn,3),p(i_burn),5,'r')
        title('$\theta_2$-$\theta_3$ distribution (log)','interpreter','latex')
        subplot(2,3,3)
        scatter3(th_T(i_burn,1),th_T(i_burn,2),th_T(i_burn,3),5,'r')
        hold on
        th_true = [log(0.1), log(4), -1];
        scatter3(th_true(1),th_true(2),th_true(3),100,'k')
        hold off
        axis equal
        title(['uncertainty window',sprintf('solution, (true $\\theta$ = %.3f, %.3f, %.3f)',[log(0.1), log(4), -1])],'interpreter','latex')
        xlabel('\theta_1');ylabel('\theta_2');zlabel('\theta_3')
        legend({'MC \theta''s','true \theta'})
        subplot(2,3,6)
        plot(th_T) % plot Markov chain
        legend({'\theta_1','\theta_2','\theta_3'})
        th_m = mean(th_T(i_burn,:));
        title(sprintf('predicted mean: %.3f, %.3f, %.3f',th_m),'interpreter','latex')
end