function model = buildModel(model)
% Output:
%   model.y           = data
%   model.modelFun    = QoI

% generate solution
if model.testCase >= 6
    [samps, model.C0, t, x] = Solvers(model.th_true,model);
end

% define QoIs
if model.testCase < 6
    model.modelFun = @(th) myModel(th,model); % model
    U = model.modelFun(model.th_true);
    Un = U + randn(size(U))*model.noiseLevel; % additive noise
    %Un = U .* (1 + randn(size(U))*model.noiseLevel); % multiplicative noise
    model.y = Un;
else
    if model.testCase >= 6 && model.testCase <= 8
        samps0 = samps; % samples w\o noise
        samps = samps + randn(size(samps))*model.noiseLevel; % noisy samples
    end
    model.t = t; model.x = x;
    model = qoiInit(samps,model);
    U = QoIs(samps,model.qoi);
    model.y = U;
    model.modelFun = @(th) myModel(th,model); % model
end

% visualization
figure(1)
switch model.testCase
    case 2
        subplot(2,2,1), surf(U), title('solution')
        subplot(2,2,2), surf(Un), %title(['add noise, $\|\epsilon^*\|_2/\sqrt{m} =$',num2str(norm(U(:)-Un(:))/sqrt(numel(U)))],'interpreter','latex')
        title(['add noise, $\|\epsilon^*\|_2^2 =$',num2str(norm(U(:)-Un(:))^2)],'interpreter','latex')
    case {3, 5}
        subplot(2,3,1), surf(U), title('solution')
        subplot(2,3,4), surf(Un), %title(['add noise, $\|\epsilon^*\|_2/\sqrt{m} =$',num2str(norm(U(:)-Un(:))/sqrt(numel(U)))],'interpreter','latex')
        title(['add noise, $\|\epsilon^*\|_2^2 =$',num2str(norm(U(:)-Un(:))^2)],'interpreter','latex')
    case {6, 7}
        x = x(model.qoi.xsampind);
        U0 = QoIs(samps0,model.qoi);
        Ut = U0(:,1:end/2); Vt = U0(:,1+end/2:end);
        Utn = U(:,1:end/2); Vtn = U(:,1+end/2:end);
        subplot(2,3,1), plot(x,[Ut(:,end),Vt(:,end)]), title(['solution, true \theta =',num2str(model.th_true)])
        subplot(2,3,4), plot(x,[Utn(:,end),Vtn(:,end)])
        title(['add noise, $\|\epsilon^*\|_2^2 =$',num2str(norm(U0(:)-U(:))^2)],'interpreter','latex')
    case 8
        Ut = samps0(:,1:end/2); Vt = samps0(:,1+end/2:end);
        Utn = samps(:,1:end/2); Vtn = samps(:,1+end/2:end);
        subplot(2,3,1), plot(x,[Ut(:,end),Vt(:,end)]), title(sprintf('solution, true $\\theta$ = %.3f, %.3f, %.1f',model.th_true),'interpreter','latex')
        subplot(2,3,4), plot(x,[Utn(:,end),Vtn(:,end)])
        title(['add noise, $\|\epsilon^*\|_2^2 =$',num2str(norm(samps0(:)-samps(:))^2)],'interpreter','latex')
    case {9, 10}
        Ut = samps(:,1:end/2); Vt = samps(:,1+end/2:end);
        subplot(1,3,1), plot(x,reshape([Ut(:,end),Vt(:,end)],[],2)), title(sprintf('solution, true $\\theta$ = %.3f, %.3f, %.1f, %.1f',model.th_true),'interpreter','latex')
        
        bar_centers = model.qoi.edges(2:end)-model.qoi.dx;
        subplot(1,3,2), for j = 1:size(U,2), bar(bar_centers,U(:,j)); hold on, end; hold off; xlim([0,15]);
    case 11
        Ct = samps; i_snap = model.qoi.i_snap;
        subplot(1,3,1), plot(x,Ct(:,i_snap)), title(sprintf('solution, true $\\theta$ = %.3f, %.2f, %.2f, %.2f',model.th_true),'interpreter','latex')
        bar_centers = model.qoi.edges(2:end)-model.qoi.dx;
        subplot(1,3,2), for j = 1:size(U,2), bar(bar_centers,U(:,j)); hold on; end; hold off; xlim([0,15]);
    case 12
        Ct = samps; i_snap = model.qoi.i_snap;
        subplot(1,2,1), plot(x,Ct(:,i_snap)), title(sprintf('solution, true $\\theta$ = %.3f, %.2f, %.2f, %.2f',model.th_true),'interpreter','latex')
        hold on; 
        if model.qoi.imax == 1, plot(x,U(end-1:end)+0*x,'--'); 
        else, plot(x,U(end)+0*x,'--');  end; hold off;
end
drawnow % plot it

function U = myModel(th,model)
% model functions

if model.testCase >= 6
    samps = Solvers(th,model); % generate solution
    U = QoIs(samps,model.qoi); % compute QoI
else
    switch model.testCase
        case 1
            U = sin(th);
        case 2
            U = reactDiffuse1d([th(1), th(1), th(2)]);
        case {3, 4}
            U = reactDiffuse1d(th);
        case 5
            U = convectReactDiffuse1d([th,0]);
    end
end