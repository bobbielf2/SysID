%% Experiment with Statistical QoI
% generate data
[~, C0, ~, ~] = cahnhilliard1d([],[0.5,-1,0,1]); % get initial condition
%%
figure(1);clf; set(gcf,'position',[0 0 1000 1000]);
figure(2);clf; set(gcf,'position',[500 0 1000 1000]);
splt = 1;
for tha = -1.1:0.05:-0.9
    for thb = 0.8:0.1:1.2
        [U, ~, t, x] = cahnhilliard1d(C0,[0.5,tha,0,thb]); % Cahn-Hilliard 1 species
        dx = mean(diff(x));
        sample = U(:,end);
        
        %%%
        figure(1), subplot(5,5,splt)
        plot(x,sample); hold on
        if splt <= 5, title(['\theta_4 = ',num2str(thb)]), end
        if mod(splt, 5)==1, ylabel(['\theta_2 = ',num2str(tha)]), end
        %xticklabels([])
        
        % normalize concentration
        
        cmin = min(sample); cmax = max(sample);
        sample_sc = (sample - cmin)/(cmax - cmin); % normalized sample to [0,1]
        plot(x,cmin+0*x,'--k')
        plot(x,cmax+0*x,'--k')
        
        % size distribution QoI
        
        threshold = 0.5;
        plot(x,(1-threshold)*cmin+threshold*cmax+0*x,'k')
        
        A = (sample_sc >= threshold); % find continuous components
        A = diff(A); % starting & end points of components (starting labeled 1, ending labeled -1)
        startpt = find(A == 1); endpt = find(A == -1);
        %%%%%% ignore incomplete components %%%%%%%%%%%%%%
        if endpt(1) < startpt(1), endpt(1) = []; end % remove incomplete component at left bdry
        if startpt(end) > endpt(end), startpt(end) = []; end % remove incomplete component  at right bdry
        %%%%%% alternatively, include the incomplete components %%%%%%%%%%%
        %if endpt(1) < startpt(1), startpt = [x(1); startpt]; end % include incomplete component at left bdry
        %if startpt(end) > endpt(end), endpt = [endpt;x(end)]; end % include incomplete component at right bdry
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        S = (endpt - startpt)*dx; % sizes
        
        edges = [0:2*dx:max(x)-min(x)]; % bins of histogram
        figure(2), subplot(5,5,splt)
        histogram(S,edges)
        if splt <= 5, title(['\theta_4 = ',num2str(thb)]), end
        if mod(splt, 5)==1, ylabel(['\theta_2 = ',num2str(tha)]), end
        xlim([0,10])
        [vals, edges] = histcounts(S,edges);
        splt = splt + 1;
    end
end
%% sensitivity of size distribution QoI vs Cahn-Hilliard parameters
load sysid_statsQoI_sizDis_testcase6_cahnhilliard
ths = {model.th_true
    [log(0.5)+0.427,-2,2]
    [log(0.5)+0.67572,-3,3]
    [log(0.5)+0.994,-5,5]
    [log(0.5)+1.768,-10,10]};
figure(1)
for i = 1:5
thth = ths{i};
U = model.modelFun(thth);
pp = model.posterior(thth);
edges = model.edges; % bins of histogram for the sizes
subplot(2,5,i),bar(edges,U); xlim([0,10]), title(sprintf('$\\theta$ =[%.4f,%d,%d], p = %.2f',thth,pp),'interpreter','latex')

thth(1) = exp(thth(1));
[C,~,t,x] = cahnhilliard1d(model.C0,thth,model.th_ind);
subplot(2,5,i+5),for j = size(C,2), plot(x,C(:,j)), drawnow; end
end