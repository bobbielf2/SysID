%% MCMC initialize

test_case   = 5;                                % test case num
noise_level = 0.02;                              % noise (percentage)
model       = testCase(test_case,noise_level);  % generate test case

% Metropolis-Hastings Iteration

model = metropolis_hastings(model);


% plot resulted posterior distr
th      = cell2mat(model.th); % all the theta tried
p       = model.p; p(isnan(p)) = 0; % posterior prob corresp to th
i_burn  = floor(numel(p)/4):numel(p); % burn-in
th_T    = cell2mat(model.th_T); % Markov chain of theta's
%close all
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
        title(['predicted modes (in desc order): ',num2str(th_mode(1:6)')])
        subplot(1,2,2),cla
        plot(th_T) % plot Markov chain
    case 2
        subplot(1,2,1)
        scatter3(th(i_burn,1),th(i_burn,2),p(i_burn),5,'r')
        title(['predicted mean: ',num2str(mean(th_T(i_burn,:)))])
        subplot(1,2,2)
        plot(th_T) % plot Markov chain
    case 3
        subplot(2,2,1)
        scatter3(th(i_burn,1),th(i_burn,2),p(i_burn),5,'r')
        title(['predicted mean: ',num2str(mean(th_T(i_burn,:)))])
        subplot(2,2,2)
        scatter3(th(i_burn,2),th(i_burn,3),p(i_burn),5,'r')
        subplot(2,2,3)
        scatter3(th_T(i_burn,1),th_T(i_burn,2),th_T(i_burn,3),5,'r')
        axis equal
        xlabel('\theta_1');ylabel('\theta_2');zlabel('\theta_3')
        subplot(2,2,4)
        plot(th_T) % plot Markov chain
    case 5
        subplot(1,2,1)
        scatter3(th(i_burn,1),th(i_burn,2),p(i_burn),5,'r')
        title(['predicted mean: ',num2str(mean(th_T(i_burn,:)))])
        subplot(1,2,2)
        plot(th_T) % plot Markov chain
end