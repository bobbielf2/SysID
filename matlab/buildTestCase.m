function model = buildTestCase(testCase,niter,noise_level,ifSparPrior,likelihoodType)

if nargin < 5, likelihoodType = 1; end
if nargin < 4, ifSparPrior = 0; end
if nargin < 3, noise_level = 0.01; end
if nargin < 2, niter = 3000; end

% Initialization
model.testCase          = testCase;
model.noiseLevel        = noise_level; % noise level
model.isSparse          = ifSparPrior; % use sparse inducing prior?
model.likelihoodType    = likelihoodType; % what type of likelihood? see postProb.m
model.modelFun          = @(th) myModel(th,testCase); % model
model.y                 = genData(model); % data
model.niter             = niter; % num of MCMC interations

switch testCase
    case 1 % solve equation: y = sin(theta)
        model.mu_th     = 0; % mean & std for prior
        model.sig_th    = pi;
        model.mu_eps    = 0; % mean & std for epsilon
        model.sig_eps   = max(noise_level,1e-4);
    case 2  % system id: react-diffuse eq. th=prefactors in reaction terms to be identified
        model.mu_th     = [0, 0];   % mean & std for prior
        model.sig_th    = [1, 1]*0.1;
        model.mu_eps    = 0;   % mean & std for epsilon
        model.sig_eps   = max(noise_level,1e-4);
    case 3  % system id: react-diffuse eq. th=prefactors in reaction terms to be identified
        model.mu_th     = [0, 0, 0];   % mean & std for prior
        model.sig_th    = [1, 1, 1]*0.2;
        model.mu_eps    = 0;   % mean & std for epsilon
        model.sig_eps   = max(noise_level,1e-4);
    case 5 % sys id: convect-react-diffuse eq. th(1)=diffusivity, th(2)=convect vel, th(3:4)=react terms
        model.mu_th     = [0, 0, 0];   % mean & std for prior
        model.sig_th    = [1, 1, 1]*0.3; % make convection velocity have smaller variance
        model.mu_eps    = 0;   % mean & std for likelihood
        model.sig_eps   = max(noise_level,1e-4);
end

model.propose   = @(th_t) propose(th_t,model);  % proposal algorithm
model.posterior = @(th) posterior(th,model);    % posterior

% Initialize
model.th        = cell(niter,1);    % theta samples
model.p         = zeros(niter,1);   % posterior probs

th0             = model.mu_th;      % pick a starting point
model.th{1}     = th0;              % initial theta
model.p(1)      = model.posterior(th0); % initial posterior


function U = myModel(th,test_case)
% model functions

switch test_case
    case 1
        U = sin(th);
    case 2
        U = reactDiffuse1d([th(1), th(1), th(2)]);
    case {3, 4}
        U = reactDiffuse1d(th);
    case 5
        U = convectReactDiffuse1d([th,0]);
end