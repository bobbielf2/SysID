function p = postProb(th,model)

y_th = model.modelFun(th); % compute y(theta) from my model
epsilon = norm(model.y - y_th); % y = y_th + eps, where eps is assumed normally distr'd

p_th = mvnpdf(th,model.mu_th,model.sig_th.^2); % prior (assume std normal distr for now)
% pd_eps = makedist('Normal','mu',model.mu_eps,'sigma',model.sig_eps); % assume std normal distr for now
% p_eps = pdf(pd_eps,epsilon);
p_eps = mvnpdf(epsilon,model.mu_eps,model.sig_eps.^2); % likelihood (assume std normal distr for now)


p = p_eps * p_th; % posterior ~ likelihood * prior
