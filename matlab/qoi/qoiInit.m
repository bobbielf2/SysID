function model = qoiInit(type,samps,model)
% Initialize QoI computation
% Input
%   type:       type of QoI
%   samps:      solution field to be sampled, ordered as "timestep-fast, species-slow" manner
%               e.g. say have n=2 species s1, s2, and k=3 snapshots i1, i2, i3
%               then order is samps = [s1(i1), s1(i2), s1(i3), s2(i1), s2(i2), s2(i3)]
%   n_spec:     number of species
% Output
%   model.qoi   struct containing parameters needed to compute the QoI.

if isfield(model,'n_spec')
    n_spec = model.n_spec;
else
    n_spec = 1;
end

Nt      = size(samps,2)/n_spec;                 % number of timesteps

if isfield(model,'n_snap')
    n_snap  = model.n_snap;                     % number of snapshots
    i_snap  = ceil((1:n_snap).'/n_snap * Nt);   % picking sampled snapshots
    i_snap  = i_snap + (0:n_spec-1)*Nt; i_snap = i_snap(:); % concat snapshots for multi species
    qoi.i_snap = i_snap;        % save snapshots numbers
    qoi.n_snap = model.n_snap;  % number of snapshots to pick
end

type = model.qoiType;       % type of QoI
qoi.type = type;     
switch type
    case 'SizeDistr'
        x = model.x; dx = mean(diff(x)); % spatial grid size
        qoi.edges = (0:2*dx:max(x)-min(x)).'; % bins of histogram for the sizes
        qoi.threshold = model.threshold;
        qoi.dx = dx; % grid spacing
    case 'MatrixComp'
        qoi.imax = 1; % also use mean(local max) as QoI?
    case 'ThresholdArea'
        N = numel(i_snap);
        qoi.cmin = zeros(N,1); qoi.cmax = zeros(N,1);
        for i = 1:N
            sample = samps(:,i_snap(i));
            cmin = min(sample); qoi.cmin(i) = cmin;
            cmax = max(sample); qoi.cmax(i) = cmax;
        end
        nt = 101; thresholds = linspace(0,1,nt);
        qoi.thresholds = thresholds;
    case 'Downsampled'
        tgap = model.tgap;
        xgap = model.xgap;
        tsampind1 = ceil(1:tgap:Nt).'; % t sampling indices of 1 species
        tsamps = tsampind1 + (0:n_spec-1)*Nt; % t sampling in all species
        qoi.tsampind = tsamps(:); % t sampling indices of all species
        qoi.xsampind = ceil(1:xgap:size(samps,1)); % x sampling indices
end
model.qoi = qoi;
end