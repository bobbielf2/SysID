function model = buildModel(model)
% Output:
%   model.y           = data
%   model.modelFun    = QoI 

switch model.testCase
    case 1
        th_true = 2; % the true theta
    case {2, 3, 4}
        th_true = [1,1,0];
    case 5
        th_true = [1,1/5,-1,0];
    case 6
        th_true = [1/10, 40/10, 0.1, -1, 1, 0.9, -1];
        [~, ~, model.UV0] = reactDiffuse1d2sp; % get initial condition
    case 7
        th_true = [1/10*10, 40/10, -1];
        [~, ~, model.UV0] = reactDiffuse1d2sp; % get initial condition
end

model.modelFun = @(th) myModel(th,model); % model
U = model.modelFun(th_true);
Un = U + randn(size(U))*model.noiseLevel; % additive noise
%Un = U .* (1 + randn(size(U))*model.noiseLevel); % multiplicative noise
model.y = Un;

% visualization
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
    case {6, 7}
        figure(1)
        subplot(2,3,1), plot(reshape(U(:,end),[],2)), title('solution')
        subplot(2,3,4), plot(reshape(Un(:,end),[],2))
        title(['add noise, $\|\epsilon^*\|_2^2 =$',num2str(norm(U(:)-Un(:))^2)],'interpreter','latex')
        drawnow
end


function U = myModel(th,model)
% model functions

switch model.testCase
    case 1
        U = sin(th);
    case 2
        U = reactDiffuse1d([th(1), th(1), th(2)]);
    case {3, 4}
        U = reactDiffuse1d(th);
    case 5
        U = convectReactDiffuse1d([th,0]);
    case 6
        U = reactDiffuse1d2sp(th,model.UV0);
    case 7
        U = reactDiffuse1d2sp([th(1)/10,th(2),0.1, th(3), 1, 0.9, -1],model.UV0);
end