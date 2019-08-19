% System identification using Bayesian inference
setup

% MCMC initialize

test_case       = 13;    % test case num
niter           = 200; % num of iterations
noise_level     = 0.01; % noise (used as std of pointwise gaussian noise)
sparse_prior    = 0;    % use Laplace prior?   
likelihood_type = 0;    % type of likelihood? 0=multivar Gaussian, 1=square of 2-norm of err (Gamma), 2=mean of err (Normal)
fixInit         = 1;    % remember initial condition?
%%%%%% tune proposal sigma
tuneMixing      = 0.5;    % tune proposal sigma to improve mixing


% generate test case
model = buildTestCase(test_case,niter,noise_level,sparse_prior,likelihood_type,fixInit,tuneMixing);

% Metropolis-Hastings Iteration
model = metropolis_hastings(model);

%% plot resulted posterior distr

% N = 500; model = metropolis_hastings(model,N); % continue for N more iters

th      = cell2mat(model.th); % all the theta tried
p       = model.p; p(isnan(p)) = 0; % posterior prob corresp to th
i_burn  = min(floor(numel(p)/4),500):numel(p); % burn-in
th_T    = cell2mat(model.th_T); % Markov chain of theta's
th_m    = mean(th_T(i_burn,:)); % predicted mean
th_true = model.th_true; % true param

figure(1)
switch test_case
    case 1
        subplot(1,2,1), scatter(th,p,5,'r')
        %scatter(th(i_burn),p(i_burn),5,'r')
        hold on, fplot(@sin,[-10,10])
        fplot(@(x) 0*x+sin(2),[-10,10])
        % find modes
        [th_mode,i_sort] = sort(th); p_mode = p(i_sort);                    % sort p in ascending order of th
        [p_mode ,loc]    = findpeaks(p_mode); th_mode = th_mode(loc);       % find modes of th based on local max of p
        [p_mode,i_sort] = sort(p_mode,'desc'); th_mode = th_mode(i_sort);   % sort modes in descending order
        th_mode = th_mode(1:min(end,6)); % pick 6 modes max
        title(sprintf(['predicted modes (in desc order): ', repmat('%.2f, ', 1,length(th_mode))],th_mode))
        subplot(1,2,2), plot(th_T,'.-') % plot Markov chain
    case 2
        subplot(2,2,3)
        scatter3(th(i_burn,1),th(i_burn,2),p(i_burn),5,'r')
        subplot(2,2,4), plot(th_T,'.-') % plot Markov chain
        [~, mu_eps] = model.posterior(th_m);
        title(['predicted mean: ',num2str(mean(th_T(i_burn,:)))])
        xlabel(['likelihood mean $\mu_\epsilon=$',num2str(mu_eps)],'interpreter','latex')
    case {3, 4, 5}
        subplot(2,3,2)
        scatter3(th(i_burn,1),th(i_burn,2),p(i_burn),5,'r')
        title('$\theta_1$-$\theta_2$ distribution','interpreter','latex')
        subplot(2,3,5)
        scatter3(th(i_burn,2),th(i_burn,3),p(i_burn),5,'r')
        title('$\theta_2$-$\theta_3$ distribution','interpreter','latex')
        subplot(2,3,3)
        scatter3(th_T(i_burn,1),th_T(i_burn,2),th_T(i_burn,3),5,'r'); axis equal;
        title('uncertainty window')
        xlabel('\theta_1');ylabel('\theta_2');zlabel('\theta_3')
        subplot(2,3,6), plot(th_T,'.-') % plot Markov chain
        legend({'\theta_1','\theta_2','\theta_3'})
        title(sprintf('predicted mean: %.3f, %.3f, %.3f',th_m),'interpreter','latex')
        [~, mu_eps] = model.posterior(th_m);
        xlabel(['likelihood mean $\mu_\epsilon=$',num2str(mu_eps)],'interpreter','latex')
        %ylim([-1,1]*0.002-1); ylim([-1,1]*0.05+4); % case 7 debug
        
        % subplot(2,3,5)
        % U = reactDiffuse1d(th_m);
        % surf(U), title('predicted solution')
    case 6
        subplot(2,3,2)
        scatter3(th(i_burn,1),th(i_burn,2),p(i_burn),5,'r')
        title('$\theta_1$-$\theta_2$ distribution','interpreter','latex')
        subplot(2,3,5)
        scatter3(th(i_burn,2),th(i_burn,3),p(i_burn),5,'r')
        title('$\theta_2$-$\theta_3$ distribution','interpreter','latex')
        subplot(2,3,3)
        scatter3(th_T(i_burn,1),th_T(i_burn,2),th_T(i_burn,3),5,'r'); axis equal;
        title('uncertainty window')
        xlabel('\theta_1');ylabel('\theta_2');zlabel('\theta_3')
        subplot(2,3,6), plot(th_T,'.-') % plot Markov chain
        legend({'\theta_1','\theta_2','\theta_3','\theta_4','\theta_5','\theta_6','\theta_7'})
        title(sprintf('predicted mean: %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f',th_m),'interpreter','latex')
        [~, mu_eps] = model.posterior(th_m);
        xlabel(['likelihood mean $\mu_\epsilon=$',num2str(mu_eps)],'interpreter','latex')
    case 7
        subplot(2,3,2)
        scatter3(th(i_burn,1),th(i_burn,2),p(i_burn),5,'r')
        title('$\theta_1$-$\theta_2$ distribution (log)','interpreter','latex')
        subplot(2,3,5)
        scatter3(th(i_burn,2),th(i_burn,3),p(i_burn),5,'r')
        title('$\theta_2$-$\theta_3$ distribution (log)','interpreter','latex')
        subplot(2,3,3)
        scatter3(th_T(i_burn,1),th_T(i_burn,2),th_T(i_burn,3),5,'r'); axis equal;
        title('uncertainty window')
        xlabel('\theta_1');ylabel('\theta_2');zlabel('\theta_3')
        subplot(2,3,6), plot(th_T,'.-') % plot Markov chain
        legend({'\theta_1','\theta_2','\theta_3'})
        title(sprintf('predicted mean: %.3f, %.3f, %.3f',th_m),'interpreter','latex')
        [~, mu_eps] = model.posterior(th_m);
        xlabel(['likelihood mean $\mu_\epsilon=$',num2str(mu_eps)],'interpreter','latex')
        %ylim([-1,1]*0.002-1); ylim([-1,1]*0.05+4); % case 7 debug
        
        % subplot(2,3,5)
        % U = reactDiffuse1d(th_m);
        % surf(U), title('predicted solution')
    case 8
        subplot(2,3,2)
        scatter3(th(i_burn,1),th(i_burn,2),p(i_burn),5,'r')
        title('$\theta_1$-$\theta_2$ distribution (log)','interpreter','latex')
        subplot(2,3,5)
        scatter3(th(i_burn,2),th(i_burn,3),p(i_burn),5,'r')
        title('$\theta_2$-$\theta_3$ distribution (log)','interpreter','latex')
        subplot(2,3,3)
        scatter3(th_T(i_burn,1),th_T(i_burn,2),th_T(i_burn,3),5,'r'); hold on;
        scatter3(th_true(1),th_true(2),th_true(3),100,'k'); hold off; axis equal;
        title(['uncertainty window',sprintf('solution, (true $\\theta$ = %.3f, %.3f, %.3f)',[log(0.1), log(4), -1])],'interpreter','latex')
        xlabel('\theta_1');ylabel('\theta_2');zlabel('\theta_3')
        legend({'MC \theta''s','true \theta'})
        subplot(2,3,6), plot(th_T,'.-') % plot Markov chain
        legend({'\theta_1','\theta_2','\theta_3'})
        title(sprintf('predicted mean: %.3f, %.3f, %.3f',th_m),'interpreter','latex')
    case {9, 10}
        if 0
            [Ct, ~, ~, x] = Solvers(model.th_true,model);
            subplot(1,3,2), plot(x,Ct(:,model.qoi.i_snap))
            title(sprintf('solution, true $\\theta$ = %.3f, %.3f, %.3f',model.th_true),'interpreter','latex')
            U = QoIs(Ct,model.qoi);
            edges = model.qoi.edges(2:end) - model.qoi.dx; % size distr QoI
            subplot(1,3,1), for j = 1:size(U,2), bar(edges,U(:,j)); hold on, end
            legend(strcat("species",string(1:2))); xlim([0,10]); hold off;
        end
        subplot(1,3,3), plot(th_T,'.-','markersize',5) % plot Markov chain
        legend(strcat("\theta_",string(1:numel(th_m))))
        %title(sprintf('predicted mean: %.3f, %.3f, %.1f, %.1f',th_m),'interpreter','latex')
        title(sprintf('predicted mean: %.3f, %.3f, %.3f',th_m),'interpreter','latex')
    case 11
        if 0
            [Ct, ~, ~, x] = Solvers(model.th_true,model);
            subplot(1,3,1), plot(x,Ct(:,model.qoi.i_snap))
            title(sprintf('solution, true $\\theta$ = %.3f, %.2f, %.2f',model.th_true),'interpreter','latex')
            U = QoIs(Ct,model.qoi);
            edges = model.qoi.edges(2:end) - model.qoi.dx; % size distr QoI
            subplot(1,3,2), for j = 1:size(U,2), bar(edges,U(:,j)); hold on, end
            legend(strcat("snap",string(1:5))); xlim([0,10]); hold off;
        end
        subplot(1,3,3), plot(th_T,'.-') % plot Markov chain
        legend(strcat("\theta_",string(1:numel(th_m))))
        title(sprintf('predicted mean: %.3f, %.2f, %.2f, %.2f',th_m),'interpreter','latex')
    case 12
        if 0
            [Ct, ~, ~, x] = Solvers(model.th_true,model);
            U = QoIs(Ct,model.qoi);
            subplot(1,2,1), plot(x,Ct(:,end)); hold on;
            title(sprintf('solution, true $\\theta$ = %.3f, %.2f, %.2f',model.th_true),'interpreter','latex')
            if model.qoi.imax == 1, plot(x,U(end-1:end)+0*x,'--','linewidth',1);
            else, plot(x,U(end)+0*x,'--','linewidth',1);  end; hold off
            legend({'snapshot','loc min','loc max'})
        end
        subplot(1,2,2), plot(th_T,'.-') % plot Markov chain
        legend(strcat("\theta_",string(1:numel(th_m))))
        title(sprintf('predicted mean: %.3f, %.2f, %.2f, %.2f',th_m),'interpreter','latex')
    otherwise
        if 0
            [Ct, ~, ~, x] = Solvers(model.th_true,model);
            U = QoIs(Ct,model.qoi);
            subplot(1,3,1), 
            for i = model.qoi.i_snap
                plot(x,Ct(:,i)); hold on;
            end
            title(sprintf('solution, true $\\theta$ = %.3f, %.2f, %.2f',model.th_true),'interpreter','latex')
            %if model.qoi.imax == 1, plot(x,U(end-1:end)+0*x,'--','linewidth',1);
            %else, plot(x,U(end)+0*x,'--','linewidth',1);  end;
            %legend({'snapshot','loc min','loc max'})
            hold off
        end
        subplot(1,3,2), plot(th_T,'.-') % plot Markov chain
        legend(strcat("\theta_",string(1:numel(th_m))))
        title(sprintf('predicted mean: %.3f, %.2f, %.2f, %.2f',th_m),'interpreter','latex')
        t0 = 1; accept = model.accept(t0:end);
        subplot(1,3,3), plot((accept(2:end)-accept(1))./(1:numel(accept)-1)') % plot acceptance rate
        %subplot(1,3,3), plot(accept,'-') % plot acceptance number
        title('acceptance rate','interpreter','latex')
        
end
%%
figure(2)
theta = model.th_true;

[Ct, ~, ~, x] = Solvers(theta,model);
subplot(2,3,2), plot(x,Ct(:,model.qoi.i_snap))
title(sprintf('true $\\theta^*$ = %.3f, %.3f, %.3f',model.th_true),'interpreter','latex')
U = QoIs(Ct,model.qoi);
edges = model.qoi.edges(2:end) - model.qoi.dx; % size distr QoI
subplot(2,3,1), for j = 1:size(U,2), bar(edges,U(:,j)); hold on, end
legend(strcat("species",string(1:2))); xlim([0,20]); hold off;
title('likelihood pick $\sigma = 0.1$ ','interpreter','latex')

subplot(2,3,3), plot(th_T,'.-','markersize',5) % plot Markov chain
legend(strcat("\theta_",string(1:numel(th_m))))
%title(sprintf('predicted mean: %.3f, %.3f, %.1f, %.1f',th_m),'interpreter','latex')
title(sprintf('stuck at: $\\theta$ =%.3f, %.3f, %.3f',th_T(end,:)),'interpreter','latex')

theta = th_T(end,:);

[Ct, ~, ~, x] = Solvers(theta,model);
subplot(2,3,5), plot(x,Ct(:,model.qoi.i_snap))
title(sprintf('trap $\\theta$ = %.3f, %.3f, %.3f',theta),'interpreter','latex')
U = QoIs(Ct,model.qoi);
edges = model.qoi.edges(2:end) - model.qoi.dx; % size distr QoI
subplot(2,3,4), for j = 1:size(U,2), bar(edges,U(:,j)); hold on, end
legend(strcat("species",string(1:2))); xlim([0,20]); hold off;
%title('','interpreter','latex')
