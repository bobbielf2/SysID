function model = buildTestCase(testCase,niter,noise_level,ifSparPrior,likelihoodType,annealingProposal,fixInit)

if nargin < 7, fixInit = 1; end
if nargin < 6, annealingProposal = 0; end
if nargin < 5, likelihoodType = 1; end
if nargin < 4, ifSparPrior = 0; end
if nargin < 3, noise_level = 0.01; end
if nargin < 2, niter = 3000; end

% Initialization
model.testCase          = testCase;
model.noiseLevel        = noise_level; % noise level
model.fixInit           = fixInit; % remember initial condition?
model.isSparse          = ifSparPrior; % use sparse inducing prior?
model.likelihoodType    = likelihoodType; % what type of likelihood? see mcmc/posterior.m
model.annealingProposal = annealingProposal; % use annealing for the proposal algorithm? see mcmc/propose.m
model.niter             = niter; % num of MCMC interations

switch testCase
    case 1 % solve equation: y = sin(theta)
        model.th_true   = 2;
        model.mu_th     = 0; % mean & std for prior
        model.sig_th    = pi;
        model.mu_eps    = 0; % mean & std for epsilon
        model.sig_eps   = max(noise_level,1e-4);
    case 2  % system id: react-diffuse eq. th=prefactors in reaction terms to be identified
        model.th_true   = [1,1,0];
        model.mu_th     = [0, 0];   % mean & std for prior
        model.sig_th    = [1, 1]*0.1;
        model.mu_eps    = 0;   % mean & std for epsilon
        model.sig_eps   = max(noise_level,1e-4);
    case 3  % system id: react-diffuse eq. th=prefactors in reaction terms to be identified
        model.th_true   = [1,1,0];
        model.mu_th     = [0, 0, 0];   % mean & std for prior
        model.sig_th    = [1, 1, 1]*0.2;
        model.mu_eps    = 0;   % mean & std for epsilon
        model.sig_eps   = max(noise_level,1e-4);
    case 5 % sys id: convect-react-diffuse eq. th(1)=diffusivity, th(2)=convect vel, th(3:4)=react terms
        model.th_true   = [1,1/5,-1,0];
        model.mu_th     = [0, 0, 0];   % mean & std for prior
        model.sig_th    = [1, 1, 1]*0.3; % make convection velocity have smaller variance
        model.mu_eps    = 0;   % mean & std for likelihood
        model.sig_eps   = max(noise_level,1e-4);
    case 6 % sys id: 2-species react-diffuse eq. th(1:2)=diffusivities, th(3:7)=reaction params
        model.th_true   = [1/10, 40/10, 0.1, -1, 1, 0.9, -1]; % true params
        model.eqtype    = 'reactdiffuse1d2sp'; 
        model.n_spec    = 2;
        model.qoiType   = 'Downsampled';
        model.tgap      = 100;
        model.xgap      = 5;
        model.mu_th     = [1/10, 40/10, 0.1, -1, 1, 0.9, -1]; %[0, 0, 0, 0, 0, 0, 0];   % mean & std for prior
        model.sig_th    = [1, 1, 1, 1, 1, 1, 1]*0.3; % make convection velocity have smaller variance
        model.mu_eps    = 0;   % mean & std for likelihood
        model.sig_eps   = max(noise_level,1e-4);
        
    case 7 % sys id: 2-species react-diffuse eq. (less params) th(1:2)=diffusivities, th(3)=reaction param
        model.th_true   = [log(0.1), log(4), -1];
        model.th_ind    = [1,2,4];
        model.eqtype    = 'reactdiffuse1d2sp'; 
        model.n_spec    = 2;
        model.qoiType   = 'Downsampled';
        model.tgap      = 100;
        model.xgap      = 5;
        model.mu_eps    = 0;   % mean & std for likelihood
        model.sig_eps   = max(noise_level,1e-4);
    case 8 % same as case 7, but use a threshold-area statistical QoI
        model.th_true   = [log(0.1), log(4), -1];
        model.th_ind    = [1,2,4];
        model.eqtype    = 'reactdiffuse1d2sp'; 
        model.n_spec    = 2;
        model.qoiType   = 'ThresholdArea';
        model.mu_eps    = 0;   % mean & std for likelihood
        model.sig_eps   = 0.25;
        model.n_snap    = 10;
    case 9 % stats QoI: size distribution. 3 param: th(1:2)=diffusivities, th(3)=reaction
        model.th_true   = [log(0.1), log(4), -1];
        model.th_ind    = [1,2,4];
        model.eqtype    = 'reactdiffuse1d2sp'; 
        model.n_spec    = 2;
        model.qoiType   = 'SizeDistr';
        model.mu_eps    = 0;   % mean & std for likelihood
        model.sig_eps   = 0.1;
        model.threshold = 0.25;
        model.n_snap    = 1;
    case 10 % stats QoI: size distribution. 4 param: th(1:2)=diffusivities, th(3:4)=reaction param
        model.th_true   = [log(0.1), log(4), -1, 1]; 
        model.th_ind    = [1,2,4,5];
        model.eqtype    = 'reactdiffuse1d2sp'; 
        model.n_spec    = 2;
        model.qoiType   = 'SizeDistr';
        model.mu_eps    = 0;   % mean & std for likelihood
        model.sig_eps   = 0.1;
        model.threshold = 0.25;
        model.n_snap    = 1;
    case 11 % same QoI as case 10, on 1-species Cahn-Hilliard system
        model.th_true   = [log(0.5),-1,1]; 
        model.th_ind    = [1,2,4];
        model.eqtype    = 'cahnhilliard1d'; 
        model.n_spec    = 1;
        model.qoiType   = 'SizeDistr';
        model.mu_eps    = 0;   % mean & std for likelihood
        model.sig_eps   = 0.03;
        model.threshold = 0.5;
        model.n_snap    = 5;
    case 12 % matrix composition QoI
        model.th_true   = [log(0.5),-1,1]; 
        model.th_ind    = [1,2,4];
        model.eqtype    = 'cahnhilliard1d'; 
        model.n_spec    = 1;
        model.qoiType   = 'MatrixComp';
        model.mu_eps    = 0;   % mean & std for likelihood
        model.sig_eps   = 0.03;
        model.n_snap    = 5;
end

% define data and QoI
model           = buildModel(model);                % data = model.y; QoI function = model.modelFun
if ~isfield(model,'mu_th') && ~isfield(model,'sig_th')
    model.mu_th     = zeros(size(model.th_true));   % prior mean
    model.sig_th    = 3*ones(size(model.th_true));  % prior std
end

if model.annealingProposal
    model.propose   = @(th_t,iter) propose(th_t,model,iter);  % proposal algorithm with simulated annealing
else
    model.propose   = @(th_t) propose(th_t,model);  % proposal algorithm
end
model.posterior = @(th) posterior(th,model);    % posterior

% Initialize
model.th        = cell(niter,1);    % theta samples
model.p         = zeros(niter,1);   % posterior probs

th0             = model.mu_th;      % pick a starting point
model.th{1}     = th0;              % initial theta
model.p(1)      = model.posterior(th0); % initial posterior


