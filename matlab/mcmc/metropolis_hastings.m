function model = metropolis_hastings(model,niter_more)
% Metropolis-Hastings iterations
    
if nargin > 1 && model.niter < niter_more
    % continue iteration
    niter_start         = model.niter;
    niter_end           = niter_more;
    th                  = cell(niter_end,1);
    th(1:niter_start)   = model.th;
    p                   = zeros(niter_end,1);
    p(1:niter_start)    = model.p;
    th_T                = cell(niter_end,1);
    th_T(1:niter_start) = model.th_T;           % store the Markov chain {theta_t,t=1,2,...}
    model.niter         = niter_end;
else
    niter_start = 1;
    niter_end   = model.niter;
    th          = model.th;
    p           = model.p;
    th_T        = th;           % store the Markov chain {theta_t,t=1,2,...}
end

th_t        = model.th{niter_start};  % current theta
p_t         = model.p(niter_start);   % current p


for i = niter_start+1:niter_end
    if mod(i,floor((niter_end-niter_start)/10)) == 1 || (niter_end - i < 5)
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
    end
    th_T{i} = th_t;
end

model.th    = th;
model.p     = p;
model.th_T  = th_T;
