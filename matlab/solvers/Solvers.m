function [C,C0,t,x] = Solvers(th,model)

% fixed initial condition?
if model.fixInit && isfield(model,'C0')
    C0 = model.C0;
else
    C0 = [];
end

eqtype = model.eqtype; % type of equation (solver)

if ~isfield(model,'th_ind'), model.th_ind = []; end

switch eqtype
    case 'reactdiffuse1d'
    case 'reactdiffuse1d2sp'
        th(model.th_ind <= 2) = exp(th(model.th_ind <= 2));
        [Ut,Vt,C0,t,x] = reactDiffuse1d2sp(C0,th,model.th_ind);
        C = [Ut, Vt];
    case 'cahnhilliard1d'
        th(model.th_ind == 1) = exp(th(model.th_ind == 1));
        [C,C0,t,x] = cahnhilliard1d(C0,th,model.th_ind);
end