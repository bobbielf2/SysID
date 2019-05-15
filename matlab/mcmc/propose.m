function th = propose(th_t,model)
% proposal algorithm, given the current theta_t, propose the next theta
if nargin == 0, testPropose; return; end

sig = model.sig_th/1.5*0+0.15; % std deviation for the proposal distr
%sig = ones(size(th_t))*model.noise*10; % proposal distribution determined by noise_level
th = mvnrnd(th_t,sig.^2,1);

function testPropose
model.sig_th = 1;
x = zeros(1,1000);
for i = 1:1000
    x(i) = propose(3,model); % normal dist rand num with mean 3
end
figure
hist(x,50); % histogram to verify if normally distributed
title('Normal distr, 1000 rand samples (mean=3, sigma=1)')