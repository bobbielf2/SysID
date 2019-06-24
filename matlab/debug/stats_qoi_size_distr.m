%% Experiment with Statistical QoI
% generate data

[U, V, ~, x, t] = reactDiffuse1d2sp;
dx = mean(diff(x));
species{1} = U(:,end);
species{2} = V(:,end);

%%
figure(1);clf; set(gcf,'position',[200 200 800 300]);

s = 1;
sample = species{s};
sample = sample + 0.01*randn(size(sample));
subplot(1,2,1)
plot(x,sample); hold on

% normalize concentration

cmin = min(sample); cmax = max(sample);
sample_sc = (sample - cmin)/(cmax - cmin); % normalized sample to [0,1]
plot(x,cmin+0*x,'--k')
plot(x,cmax+0*x,'--k')

% size distribution QoI

threshold = 0.25;
plot(x,(1-threshold)*cmin+threshold*cmax+0*x,'k')

A = (sample_sc >= threshold); % find continuous components
A = diff(A); % starting & end points of components (starting labeled 1, ending labeled -1)
startpt = find(A == 1); endpt = find(A == -1);
% ignore incomplete components
if endpt(1) < startpt(1), endpt(1) = []; end % remove incomplete component at left bdry
if startpt(end) > endpt(end), startpt(end) = []; end % remove incomplete component  at right bdry
%%%%%% alternatively, include the incomplete components %%%%%%%%%%%
% if endpt(1) < startpt(1), startpt = [x(1); startpt]; end % include incomplete component at left bdry
% if startpt(end) > endpt(end), endpt = [endpt;x(end)]; end % include incomplete component at right bdry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = (endpt - startpt)*dx; % sizes

edges = [0:2*dx:max(x)-min(x)]; % bins of histogram
subplot(1,2,2)
histogram(S,edges)
[vals, edges] = histcounts(S,edges);

%%
load('sysid_statsQoI_sizDis_testcase1.mat','model')
th_T = cell2mat(model.th_T); % Markov chain of theta's
i_burn  = floor(model.niter/2):model.niter; % burn-in
%i_burn = 100:model.niter;
th_true = model.th_true; % actual theta
x = model.x; dx = mean(diff(x)); % spatial domain
edges = [0:2*dx:max(x)-min(x)].'; edges = edges(2:end)-dx; % bins for the size distribution
thd = model.threshold;        

figure; set(gcf,'position',[200 200 800 800]);

subplot(2,2,1)
bar(edges,model.y(:,1)), xlim([0, 10])
title('species 1')

subplot(2,2,2)
bar(edges,model.y(:,2),'r'), xlim([0, 10])
title('species 2')


[U, V] = reactDiffuse1d2sp([1/10, 40/10, 0.1, -1, 1, 0.9, -1],model.UV0);

subplot(2,2,3)
sample = U(:,end);
cmax = max(sample); cmin = min(sample); threshold = (1-thd)*cmin+thd*cmax;
plot(x,sample), hold on
plot(x,0*x + threshold, '--k')
legend({'species 1', 'threshold'})
title(['species 1',sprintf('solution, (true $\\theta$ = %.3f, %.3f, %.3f)',th_true)],'interpreter','latex')

subplot(2,2,4)
sample = V(:,end);
cmax = max(sample); cmin = min(sample); threshold = (1-thd)*cmin+thd*cmax;
plot(x,sample,'r'), hold on
plot(x,0*x + threshold, '--k')
legend({'species 2', 'threshold'})
title(['species 2',sprintf('solution, (true $\\theta$ = %.3f, %.3f, %.3f)',th_true)],'interpreter','latex')

% subplot(2,2,4)
% plot(th_T) % plot Markov chain
% legend({'\theta_1','\theta_2','\theta_3'})
% th_m = mean(th_T(i_burn,:));
% title(sprintf('predicted mean: %.3f, %.3f, %.3f',th_m),'interpreter','latex')
