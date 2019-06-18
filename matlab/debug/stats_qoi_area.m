%% Experiment with Statistical QoI
% generate data

[U, V, ~, x, ~] = reactDiffuse1d2sp;
species{1} = U(:,end);
species{2} = V(:,end);

figure(1); set(gcf,'position',[200 200 800 300]);
for s = 1:2
    sample = species{s};
    subplot(1,2,1)
    plot(x,sample); hold on
    
    % normalize concentration
    
    cmin = min(sample); cmax = max(sample);
    sample_sc = (sample - cmin)/(cmax - cmin); % normalized sample to [0,1]
    
    % area QoI
    
    Nt = 101; threshold = linspace(0,1,Nt); % sample thresholds
    
    area = zeros(Nt, 1);
    for j = 1:Nt
        thd = threshold(j);
        area(j) = sum(sample_sc >= thd);
    end
    area = area/numel(sample);
    
    subplot(1,2,2)
    plot(threshold,area); hold on
end
subplot(1,2,1)
hold off; title('concentration')
legend({'species 1', 'species 2'})
subplot(1,2,2)
hold off; title('normalized area QoI')
legend({'species 1', 'species 2'})

%% Stability of QoIs

figure(2); clf
Nt = 101; threshold = linspace(0,1,Nt); % sample thresholds
[U, V, ~, x, dt] = reactDiffuse1d2sp;
i_snap = round(linspace(1,size(U,2),5));

N = 10; % num of samples
for i = 1:N
    [U, V, ~, x, dt] = reactDiffuse1d2sp;
    
    for it = 1:numel(i_snap)
        species1 = U(:,i_snap(it));
        species2 = V(:,i_snap(it));
        
        sample = species1;
        cmin = min(sample); cmax = max(sample);
        sample_sc = (sample - cmin)/(cmax - cmin); % normalized sample to [0,1]
        area = zeros(Nt, 1);
        for j = 1:Nt
            thd = threshold(j);
            area(j) = sum(sample_sc >= thd);
        end
        area = area/numel(sample);
        subplot(2,5,it)
        plot(threshold,area); hold on
        title(['snapshot',num2str(it)])
        
        sample = species2;
        cmin = min(sample); cmax = max(sample);
        sample_sc = (sample - cmin)/(cmax - cmin); % normalized sample to [0,1]
        area = zeros(Nt, 1);
        for j = 1:Nt
            thd = threshold(j);
            area(j) = sum(sample_sc >= thd);
        end
        area = area/numel(sample);
        subplot(2,5,it+5)
        plot(threshold,area); hold on
        title(['snapshot',num2str(it)])
    end
end

subplot(2,5,1)
ylabel('species 1')
subplot(2,5,6)
ylabel('species 2')
% subplot(1,2,1)
% hold off; title('species 1')
% legend(strcat("sample",string(1:N)))
% subplot(1,2,2)
% hold off; title('species 2')
% legend(strcat("sample",string(1:N)))

%% plot Markov chain of sysid via the stats QoIs

load('sysid_statsQoI_testcase.mat','model')
th_T = cell2mat(model.th_T); % Markov chain of theta's
i_burn  = floor(model.niter/2):model.niter; % burn-in
%i_burn = 100:model.niter;
th_true = model.th_true; % actual theta

figure; set(gcf,'position',[200 200 800 300]);
subplot(1,2,1)
plot(th_T) % plot Markov chain
legend({'\theta_1','\theta_2','\theta_3'})
th_m = mean(th_T(i_burn,:));
title(sprintf('predicted mean: %.3f, %.3f, %.3f',th_m),'interpreter','latex')
subplot(1,2,2)
scatter3(th_T(i_burn,1),th_T(i_burn,2),th_T(i_burn,3),5,'r')
hold on
scatter3(th_true(1),th_true(2),th_true(3),100,'k')
axis equal
title(['uncertainty window',sprintf('solution, (true $\\theta$ = %.3f, %.3f, %.3f)',th_true)],'interpreter','latex')
xlabel('\theta_1');ylabel('\theta_2');zlabel('\theta_3')
legend({'MC \theta''s','true \theta'})

