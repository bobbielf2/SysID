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
        [~, ~, model.UV0,x] = reactDiffuse1d2sp; % get initial condition
    case {7, 8, 9}
        th_true = [log(1/10), log(40/10), -1];
        [~, ~, model.UV0,x,t] = reactDiffuse1d2sp; % get initial condition
end


if model.testCase == 8
    [V1, V2] = reactDiffuse1d2sp([exp(th_true(1)),exp(th_true(2)),0.1, th_true(3), 1, 0.9, -1],model.UV0);
    U = [V1,V2];
    Un = U + randn(size(U))*model.noiseLevel; % additive noise
    
    % pick snapshots
    N = 10; % number of snapshots to extract for each species
    i_snap = round(linspace(1,size(V1,2),N));
    i_snap = [i_snap, i_snap+N];
    model.i_snap = i_snap;
    
    % setup sampling thresholds
    Nt = 101; thresholds = linspace(0,1,Nt);
    model.thresholds = thresholds;
    
    area = zeros(Nt, 1);
    model.y = zeros(Nt, 2*N);
    model.cmax = zeros(2*N,1); model.cmin = zeros(2*N,1);
    for i = 1:2*N
        sample = Un(:,i_snap(i));
        cmin = min(sample); model.cmin(i) = cmin;
        cmax = max(sample); model.cmax(i) = cmax;
        sample_sc = (sample - cmin)/(cmax - cmin); % normalized sample to [0,1]
        for j = 1:Nt
            thd = thresholds(j);
            area(j) = sum(sample_sc >= thd);
        end
        area = area/numel(sample_sc);
        model.y(:,i) = area;
    end
    model.modelFun = @(th) myModel(th,model); % model
elseif model.testCase == 9
    model.x = x; model.t = t;
    model.threshold = 0.25;
    model.th_true = th_true;
    model.modelFun = @(th) myModel(th,model); % model
    U = model.modelFun(th_true);
    model.y = U;
else
    model.modelFun = @(th) myModel(th,model); % model
    U = model.modelFun(th_true);
    Un = U + randn(size(U))*model.noiseLevel; % additive noise
    %Un = U .* (1 + randn(size(U))*model.noiseLevel); % multiplicative noise
    model.y = Un;
end
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
        subplot(2,3,1), plot(x,reshape(U(:,end),[],2)), title(['solution, true \theta =',num2str(th_true)])
        subplot(2,3,4), plot(x,reshape(Un(:,end),[],2))
        title(['add noise, $\|\epsilon^*\|_2^2 =$',num2str(norm(U(:)-Un(:))^2)],'interpreter','latex')
        drawnow
    case 8
        figure(1)
        subplot(2,3,1), plot(x,reshape([V1(:,end),V2(:,end)],[],2)), title(sprintf('solution, true $\\theta$ = %.3f, %.3f, %.3f',th_true),'interpreter','latex')
        subplot(2,3,4), plot(x,reshape([Un(:,end/2),Un(:,end)],[],2))
        title(['add noise, $\|\epsilon^*\|_2^2 =$',num2str(norm(U(:)-Un(:))^2)],'interpreter','latex')
        drawnow
    case 9
        figure(1)
        [Ut, Vt, ~, x, ~] = reactDiffuse1d2sp([exp(th_true(1)),exp(th_true(2)),0.1, th_true(3), 1, 0.9, -1],model.UV0);
        subplot(2,3,1), plot(x,reshape([Ut(:,end),Vt(:,end)],[],2)), title(sprintf('solution, true $\\theta$ = %.3f, %.3f, %.3f',th_true),'interpreter','latex')
        
        dx = mean(diff(x)); % spatial grid size
        edges = [0:2*dx:max(x)-min(x)].'; % bins of histogram for the sizes
        edges = edges(2:end)-dx;
        subplot(2,3,4), for j = 1:size(U,2), bar(edges,U(:,j)); hold on, end
        hold off
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
        U = reactDiffuse1d2sp([exp(th(1)),exp(th(2)),0.1, th(3), 1, 0.9, -1],model.UV0);
    case 8
        if model.fixInit
            UV0 = model.UV0;
        else
            UV0 = [];
        end
        [Ut, Vt] = reactDiffuse1d2sp([exp(th(1)),exp(th(2)),0.1, th(3), 1, 0.9, -1],UV0);
        UV = [Ut,Vt];
        
        % snapshots and threholds
        i_snap = model.i_snap; N = numel(i_snap);
        thresholds = model.thresholds;  Nt = numel(thresholds);
        
        U = zeros(Nt, N);
        area = zeros(Nt, 1);
        for i = 1:N
            sample = UV(:,i_snap(i));
            cmin = model.cmin(i); cmax = model.cmax(i);
            sample_sc = (sample - cmin)/(cmax - cmin); % normalized sample to [0,1]
            for j = 1:Nt
                thd = thresholds(j);
                area(j) = sum(sample_sc >= thd);
            end
            area = area/numel(sample_sc);
            U(:,i) = area;
        end
    case 9
        if model.fixInit % better fix initial condition!
            UV0 = model.UV0;
        else
            UV0 = [];
        end
        [Ut, Vt, ~, x, ~] = reactDiffuse1d2sp([exp(th(1)),exp(th(2)),0.1, th(3), 1, 0.9, -1],UV0);
        dx = mean(diff(x)); % spatial grid size
        UV = [Ut(:,end),Vt(:,end)]; % use last snapshot of each species, should correspond to time T >= 15
        
        % initialize QoI: size distribution
        N = size(UV,2); % num of snapshots
        edges = [0:2*dx:max(x)-min(x)]; % bins of histogram for the sizes
        U = zeros(numel(edges)-1,N);
        for j = 1:2
            % normalized sample to [0,1]
            sample = UV(:,j);
            cmin = min(sample); cmax = max(sample);
            sample_sc = (sample - cmin)/(cmax - cmin);
            
            % find concentration sizes at a threshold
            threshold = 0.25;
            A = (sample_sc >= threshold); % find continuous components
            A = diff(A); % starting & end points of components (starting labeled 1, ending labeled -1)
            startpt = find(A == 1); endpt = find(A == -1);
            if isempty(startpt) || isempty(endpt) || (numel(startpt)==1 && numel(endpt)==1 && endpt<startpt)
                U(end,j) = 1;
            else
                try
                    if endpt(1) < startpt(1), endpt(1) = []; end % remove incomplete component at left bdry
                    if startpt(end) > endpt(end), startpt(end) = []; end % remove incomplete component  at right bdry
                catch
                    keyboard
                end
                S = (endpt - startpt)*dx; % get sizes
                
                % size distribution
                vals = histcounts(S,edges);
                U(:,j) = vals(:)/sum(vals);
            end
        end
end