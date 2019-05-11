function model = metropolis_hastings(model)
% Metropolis-Hastings iterations

niter   = model.niter;
th      = model.th;
p       = model.p;

th_t    = model.th{1};  % current theta
p_t     = model.p(1);   % current p
th_T    = th;           % store the Markov chain {theta_t,t=1,2,...}

for i = 2:niter
    if mod(i,floor(niter/10)) == 1 || (niter - i < 2)
        disp(['Iter ',num2str(i)])
    end
    th{i}   = model.propose(th_t);      % propose new sample based on previous
    p(i)    = model.posterior(th{i});   % calculate posterior probability
    %if isnan(p(i)/p_t), alpha = 0; else
        alpha   = min([1, p(i)/p_t]);   % acceptance probability
    %end
    if rand() <= alpha  % accept or not
        th_t = th{i};
        p_t  = p(i);
    end
    th_T{i} = th_t;
end

model.th    = th;
model.p     = p;
model.th_T  = th_T;
