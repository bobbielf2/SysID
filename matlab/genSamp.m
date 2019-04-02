function th = genSamp(th_t,model)
% proposal algorithm, given the current theta_t, propose the next theta
if nargin == 0, testGenSamp; return; end

sig = model.sig_th/1.5; % std deviation for the proposal distr
th = mvnrnd(th_t,sig.^2,1);

function testGenSamp
model.sig_th = 1;
x = zeros(1,1000);
for i = 1:1000
    x(i) = genSamp(3,model); % normal dist rand num with mean 3
end
figure
hist(x,50); % histogram to verify if normally distributed
title('Normal distr, 1000 rand samples (mean=3, sigma=1)')