%% Experiment with Statistical QoI (matrix composition)
% generate data
[C, C0, t, x] = cahnhilliard1d([],[0.5,-1,0,1]); % get initial condition
sample = C(:,end);
%%
% find local max & min
variation = sign(diff(sample));
extrema = sign(diff(variation));
imax = find(extrema == -1) + 1;
imin = find(extrema == 1) + 1;

matcomp = mean(sample(imin)); % mean of loc mins
matcomp2 = mean(sample(imax)); % mean of loc maxs

% plot
figure(1)
plot(x,sample); hold on
scatter(x(imax),sample(imax))
scatter(x(imin),sample(imin))
plot(x,matcomp+0*x,'k--')
hold off

% use a likelihood to ensure accuracy up to 2 digits
err = 0.01;
sig = 0.005;
p = 1/sqrt(2*pi*sig^2)*exp(-err^2/2/sig^2)