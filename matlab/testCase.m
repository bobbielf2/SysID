function model = testCase(test_case,noise_level)

% Initialization
model.y         = genData(test_case,noise_level); % data
model.modelFun  = @(th) myModel(th,test_case); % model

switch test_case
    case 1 % solve equation: y = sin(theta)
        model.mu_th     = 0; % mean & std for prior
        model.sig_th    = pi;
        model.mu_eps    = 0; % mean & std for likelihood
        model.sig_eps   = 0.02;
        
        niter           = 9000; % num of MCMC interations
        model.niter     = niter;
        model.th        = cell(niter,1);   % theta samples
        model.p         = zeros(niter,1);  % posterior probs
        model.th{1}     = 0;               % pick a starting point
        model.p(1)      = postProb(model.th{1},model);
    case 2  % system id: react-diffuse eq. th=prefactors in reaction terms to be identified
        model.mu_th     = [0, 0];   % mean & std for prior
        model.sig_th    = [1, 1]*0.3;
        model.mu_eps    = 0;   % mean & std for likelihood
        model.sig_eps   = 0.5;
        
        niter           = 3000;             % num of MCMC interations
        model.niter     = niter;
        model.th        = cell(niter,1);   % theta samples
        model.p         = zeros(niter,1);  % posterior probs
        model.th{1}     = [0,0];               % pick a starting point
        model.p(1)      = postProb(model.th{1},model);
    case 3  % system id: react-diffuse eq. th=prefactors in reaction terms to be identified
        model.mu_th     = [0, 0, 0];   % mean & std for prior
        model.sig_th    = [1, 1, 1]*0.3;
        model.mu_eps    = 0;   % mean & std for likelihood
        model.sig_eps   = 0.2;
        
        niter           = 5000;             % num of MCMC interations
        model.niter     = niter;
        model.th        = cell(niter,1);   % theta samples
        model.p         = zeros(niter,1);  % posterior probs
        model.th{1}     = [0,0,0];               % pick a starting point
        model.p(1)      = postProb(model.th{1},model);
    case 5 % sys id: convect-react-diffuse eq. th(1)=diffusivity, th(2)=convect vel, th(3:4)=react terms
        model.mu_th     = [0, 0, 0];   % mean & std for prior
        model.sig_th    = [1, 1/3, 1]*0.3; % make convection velocity have smaller variance
        model.mu_eps    = 0;   % mean & std for likelihood
        model.sig_eps   = 0.1;
        
        niter           = 5000;             % num of MCMC interations
        model.niter     = niter;
        model.th        = cell(niter,1);   % theta samples
        model.p         = zeros(niter,1);  % posterior probs
        model.th{1}     = [0,0,0];               % pick a starting point
        model.p(1)      = postProb(model.th{1},model);
end

model.propose    = @(th_t) genSamp(th_t,model); % proposal algorithm
model.posterior  = @(th) postProb(th,model);    % posterior
