function [p, mu_eps, sig_eps] = posterior(th,model)

y_th = model.modelFun(th); % compute y(theta) from my model

% Compute prior
if model.isSparse
    b = model.sig_th/sqrt(2);
    p_th = exp(-norm((th-model.mu_th)./b,1))/prod(2*b);
else
    p_th = mvnpdf(th,model.mu_th,model.sig_th.^2); % prior (assume std normal distr for now)
end

% Compute likelihood for eps
%               y = y_th + eps, where eps is assumed normally distr'd
switch model.likelihoodType

    case 1 % eps 2-norm squared, likelihood ~ Gamma distribution
        eps_modified = norm(model.y(:) - y_th(:))^2; 
        gam_k = numel(y_th)/2; gam_th = 2*model.sig_eps^2; % Gamma distr params
        mu_eps = model.sig_eps^2*numel(y_th); % mean of err (for plots)
        sig_eps = model.sig_eps^2*sqrt(2*numel(y_th)); % std of err (for plots)
        p_eps = gampdf(eps_modified,gam_k,gam_th); % sum of squares of normal distr N(0,sig^2)

    case 2 % mean of eps, likelihood ~ normal distribution of smaller sigma
        eps_modified = mean(model.y - y_th, 'all'); 
        mu_eps = 0; sig_eps = model.sig_eps / sqrt(numel(y_th));
        p_eps = mvnpdf(eps_modified,mu_eps,sig_eps.^2); % likelihood (assume std normal distr for now)
end

% Compute posterior
p = p_eps * p_th; % posterior ~ likelihood * prior
