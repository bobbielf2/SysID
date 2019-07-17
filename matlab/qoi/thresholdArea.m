function U = thresholdArea(samps, qoi)
% compute the functional QoI: total x-area above a given threshold of samples

% samps = solution field to be sampled, ordered as "timestep-fast, species-slow" manner
%    e.g. say have n=2 species s1, s2, and k=3 snapshots i1, i2, i3
%    then order is samps = [s1(i1), s1(i2), s1(i3), s2(i1), s2(i2), s2(i3)]
% qoi   = struct containing parameters needed to compute the QoI. (See qoiInit.m)

i_snap = qoi.i_snap; N = numel(i_snap);  % total number of samples

% compute the threshold-area QoI
thresholds = qoi.thresholds;  nt = numel(thresholds);
U = zeros(nt, N);
area = zeros(nt, 1);
for i = 1:N
    % normalize sample to [0,1]
    sample = samps(:,i_snap(i));
    cmin = qoi.cmin(i); cmax = qoi.cmax(i);
    sample_sc = (sample - cmin)/(cmax - cmin);
    
    % area of field above threshold
    for j = 1:nt
        thd = thresholds(j);
        area(j) = sum(sample_sc >= thd);
    end
    area = area/numel(sample_sc);
    U(:,i) = area;
end
end