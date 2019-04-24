function y = genData(model)

switch model.testCase
    case 1
        th_true = 2; % the true theta
        y = model.modelFun(th_true);
        %y = y  + randn()*model.noiseLevel; % additive noise
        y = y * (1 + randn()*model.noiseLevel); % multiplicative noise
    case {2, 3, 4}
        th_true = [1,1,0];
        U = model.modelFun(th_true);
        Un = U + randn(size(U))*model.noiseLevel; % additive noise
        %Un = U .* (1 + randn(size(U))*model.noiseLevel); % multiplicative noise
        y = Un;
    case 5
        U = model.modelFun([1,1/5,-1,0]);
        Un = U + randn(size(U))*model.noiseLevel; % additive noise
        y = Un;
end

switch model.testCase
    case 2
        figure(1)
        subplot(2,2,1), surf(U), title('solution')
        subplot(2,2,2), surf(Un), %title(['add noise, $\|\epsilon^*\|_2/\sqrt{m} =$',num2str(norm(U(:)-Un(:))/sqrt(numel(U)))],'interpreter','latex')
        title(['add noise, $\|\epsilon^*\|_2^2 =$',num2str(norm(U(:)-Un(:))^2)],'interpreter','latex')
        drawnow
    case {3, 5}
        figure(1)
        subplot(2,3,1), surf(U), title('solution')
        subplot(2,3,4), surf(Un), %title(['add noise, $\|\epsilon^*\|_2/\sqrt{m} =$',num2str(norm(U(:)-Un(:))/sqrt(numel(U)))],'interpreter','latex')
        title(['add noise, $\|\epsilon^*\|_2^2 =$',num2str(norm(U(:)-Un(:))^2)],'interpreter','latex')
        drawnow
end