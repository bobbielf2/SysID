function U = matrixComposition(samps, qoi)
% compute matrix composition (mean of local minima) of samples

% samps = solution field to be sampled, ordered as "timestep-fast, species-slow" manner
%    e.g. say have n=2 species s1, s2, and k=3 snapshots i1, i2, i3
%    then order is samps = [s1(i1), s1(i2), s1(i3), s2(i1), s2(i2), s2(i3)]
% qoi   = struct containing parameters needed to compute the QoI. (See qoiInit.m)

i_snap = qoi.i_snap; N = numel(i_snap);  % total number of samples

% compute matrix composition (mean of local minima)
U = zeros(2,N); 
for i = 1:N
    sample = samps(:,i_snap(i));
    % find local max & min
    variation = sign(diff(sample));
    extrema = sign(diff(variation));
    imax = find(extrema == -1) + 1;
    imin = find(extrema == 1) + 1;
    
    if ~isempty(imin)
        U(1,i) = mean(sample(imin)); % mean of loc mins
    else
        U(1,i) = min(sample,[],'all');
    end
    if ~isempty(imax)
        U(2,i) = mean(sample(imax)); % mean of loc mins
    else
        U(2,i) = max(sample,[],'all');
    end
end

if qoi.imax
    U = U(:); % also use mean(local max) as QoI
else
    U = U(1,:).';
end

end