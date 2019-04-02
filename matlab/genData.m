function y = genData(test_case,noise_level)

switch test_case
    case 1
        y = sin(2);
        %y = y  + randn()*noise_level; % additive noise
        y = y * (1 + randn()*noise_level); % multiplicative noise
    case {2, 3}
        U = reactDiffuse1d([1,1,0]);
        U = U + randn(size(U))*noise_level; % additive noise
        %U = U .* (1 + randn(size(U))*noise_level); % multiplicative noise
        y = U(:);
end