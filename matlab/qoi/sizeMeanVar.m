function U = sizeMeanVar(samps, qoi)
% compute size distribution of samples

% samps = solution field to be sampled, ordered as "timestep-fast, species-slow" manner
%    e.g. say have n=2 species s1, s2, and k=3 snapshots i1, i2, i3
%    then order is samps = [s1(i1), s1(i2), s1(i3), s2(i1), s2(i2), s2(i3)]
% qoi   = struct containing parameters needed to compute the QoI. (See qoiInit.m)

i_snap = qoi.i_snap; N = numel(i_snap);  % total number of samples

% compute size distribution
dx = qoi.dx; % grid spacing
U = zeros(2,N);
for i = 1:N
    % normalized sample to [0,1]
    sample = samps(:,i_snap(i));
    cmin = min(sample); cmax = max(sample);
    sample_sc = (sample - cmin)/(cmax - cmin);
    
    % find concentration sizes at a threshold
    threshold = qoi.threshold;
    A = (sample_sc >= threshold); % find continuous components
    A = diff(A); % starting & end points of components (starting labeled 1, ending labeled -1)
    startpt = find(A == 1); endpt = find(A == -1);
    if isempty(startpt) || isempty(endpt) || (numel(startpt)==1 && numel(endpt)==1 && endpt<startpt)
        S = dx*numel(sample);
        vals = [S;0]; % whole domain as one bump
        U(:,i) = vals;
    else
        try
            if endpt(1) < startpt(1), endpt(1) = []; end % remove incomplete component at left bdry
            if startpt(end) > endpt(end), startpt(end) = []; end % remove incomplete component  at right bdry
        catch
            keyboard
        end
        S = (endpt - startpt)*dx; % get sizes
        
        % size distribution
        vals = [mean(S);var(S)];
        U(:,i) = vals;
    end
end

end