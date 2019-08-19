%% Experiment with Statistical QoI
sm = [];
sv = [];

%%
% generate data
model.th_true   = [log(0.1), log(4), -1];
model.th_ind    = [1,2,4];
model.eqtype    = 'reactdiffuse1d2sp';
model.n_spec    = 2;
[samps, ~, t, x] = Solvers(model.th_true,model);

% setup qoi
model.threshold = 0.5;
model.n_snap    = 1;
model.t = t; model.x = x;


%%%
figure(1);clf; set(gcf,'position',[200 200 800 300]);
% normalize concentration

for j = 1:2
sample = samps(:,model.qoi.i_snap(j));
cmin = min(sample); cmax = max(sample);
subplot(1,3,1)
plot(x,sample); hold on
plot(x,cmin+0*x,'--k')
plot(x,cmax+0*x,'--k')
threshold = model.threshold;
plot(x,(1-threshold)*cmin+threshold*cmax+0*x,'m')
end
hold off

% size distribution
subplot(1,3,2)
model.qoiType   = 'SizeDistr';
model = qoiInit(samps,model);
U = QoIs(samps,model.qoi);
for j = 1:2
    bar(model.qoi.edges(1:end-1),U(:,j)), xlim([0, 10]); hold on
end

% size mean and var
model.qoiType   = 'SizeMeanVar';
model = qoiInit(samps,model);
U = QoIs(samps,model.qoi);
for j = 1:2
    m = U(1,j);
    plot([m,m],[0,1])
end
hold off

% repeat to get more samples of size means and vars
sm = [sm;U(1,:)]; % stack sample size means
sv = [sv;U(2,:)]; % stack sample size vars
disp(['sig of mean: ',num2str(std(sm))]) % estimate sigma for the size means,~ [0.096, 0.36]
disp(['sig of var: ',num2str(std(sv))])  % estimate sigma for the size vars, ~ [0.050, 0.63]
