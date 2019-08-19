%% Experiment with Statistical QoI
% generate data
model.th_true   = [log(0.1), log(4), -1];
model.th_ind    = [1,2,4];
model.eqtype    = 'reactdiffuse1d2sp';
model.n_spec    = 2;
[samps, ~, t, x] = Solvers(model.th_true,model);

% setup qoi
model.n_snap    = 1;
model.t = t; model.x = x;
model.qoiType   = 'SizeDistr';
model.threshold = 0.5;
model.k_dx = 5;
model = qoiInit(samps,model);


%%%
figure(1);clf; %set(gcf,'position',[200 200 800 300]);
% normalize concentration

for i = 1:2
sample = samps(:,model.qoi.i_snap(i));
cmin = min(sample); cmax = max(sample);
subplot(1,3,1)
plot(x,sample); hold on
plot(x,cmin+0*x,'--k')
plot(x,cmax+0*x,'--k')
threshold = model.threshold;
plot(x,(1-threshold)*cmin+threshold*cmax+0*x,'m')
end
hold off

% plot size distribution
subplot(1,3,2)
U = QoIs(samps,model.qoi);
for i = 1:2
    bar(model.qoi.edges(1:end-1),U(:,i)), xlim([0, 10]); hold on
end
%% hellinger distance
model.k_dx = 2; % bin size = k_dx * dx
model = qoiInit(samps,model);
% base case
s = 1; % s=1 or 2, which species to pick
C0 = Solvers(model.th_true,model); U0 = QoIs(C0,model.qoi); sample0 = U0(:,s);
dists = []; % save the distance for analysis
% plot it
figure(1);clf; %set(gcf,'position',[200 200 800 300]);
edges = (model.qoi.edges(1:end-1) + model.qoi.edges(2:end))/2;
subplot(1,3,1); bar(edges,sample0), xlim([0, 10]);
title(sprintf('true distribution (binsize $%d\\Delta x$)',model.k_dx),'interpreter','latex')
% generate n samples to compare with base case
nsam = 2; % num of samples
C = cell(nsam,1);
U = cell(nsam,1); 
sample = cell(nsam,1);

for i = 1:nsam
    C{i} = Solvers(model.th_true,model);
    U{i} = QoIs(C{i},model.qoi);
    sample{i} = U{i}(:,s);
end


for i = 1:nsam
    dist = norm(sqrt(sample{i}) - sqrt(sample0))/sqrt(2); % Hellinger distance
    dists = [dists;dist]; % save for analysis
    if i <= 2
        subplot(1,3,i+1)
        bar(edges,sample{i}), xlim([0, 10]);
        title(sprintf('Hellinger dist $ = %.4g$',dist),'interpreter','latex')
    end
end
%%
histogram(dists,50)
title(sprintf('distribution of Hellinger distance'),'interpreter','latex')