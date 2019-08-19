function model = metropolis_hastings(model,niter_more)
% Metropolis-Hastings iterations
% Input
%   model: system to be learned
%   niter_more: resume/continue MCMC simulation for an additional
%                   niter_more number of iterations

if nargin > 1
    % continue iteration
    niter_start         = model.niter;
    niter_end           = niter_start + niter_more;
    th                  = cell(niter_end,1);
    th(1:niter_start)   = model.th;
    p                   = zeros(niter_end,1);
    p(1:niter_start)    = model.p;
    th_T                = cell(niter_end,1);
    th_T(1:niter_start) = model.th_T;           % store the Markov chain {theta_t,t=1,2,...}
    p_T                 = zeros(niter_end,1);
    p_T(1:niter_start)  = model.p_T;
    model.niter         = niter_end;
    accept              = zeros(niter_end,1);
    accept(1:niter_start) = model.accept;
else
    niter_start = 1;
    niter_end   = model.niter;
    th          = model.th;
    p           = model.p;
    th_T        = th;           % store the Markov chain {theta_t,t=1,2,...}
    p_T         = p;
    accept      = zeros(niter_end,1);
end

th_t        = th_T{niter_start};  % current theta
p_t         = p_T(niter_start);   % current p
accept_t    = accept(niter_start);% current number of acceptance


for i = niter_start+1:niter_end
    if mod(i,ceil((niter_end-niter_start)/10)) == 1 || (niter_end - i < 5)
        disp(['Iter ',num2str(i)])
    end
    th{i}   = model.propose(th_t);      % propose new sample based on previous
    p(i)    = model.posterior(th{i});   % calculate posterior probability
    if model.likelihoodType == 0
        alpha   = min([0, p(i) - p_t]);
        r       = log10(rand());
    else
    %if isnan(p(i)/p_t), alpha = 0; else
        alpha   = min([1, p(i)/p_t]);   % acceptance probability
        r       = rand();
    %end
    end
    if r <= alpha  % accept or not
        th_t = th{i};
        p_t  = p(i);
        accept_t = accept_t + 1;
    end
    th_T{i} = th_t;
    p_T(i) = p_t;
    accept(i) = accept_t;
    
    % adjust proposal
    if i >= 50
        model.sig_q = (2.4)^2/numel(th_t)*(cov(cell2mat(th_T(min(floor(i/2),500):i,:)))+0.1*eye(numel(th_t)));
    end
    if mod(i,50) == 1
        disp('adjust sig_q... ')
        disp(num2str(model.sig_q))
    end
end

model.th    = th;
model.p     = p;
model.th_T  = th_T;
model.p_T   = p_T;
model.accept = accept;

disp(['accept rate is ',num2str(accept(niter_end)),' out of ',num2str(niter_end)])
