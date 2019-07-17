function setup(v)

if nargin == 0
    mfilepath=fileparts(mfilename('fullpath'));
    addpath([mfilepath, '/mcmc'])
    addpath([mfilepath, '/solvers'])
    addpath([mfilepath, '/qoi'])
else
    mfilepath=fileparts(mfilename('fullpath'));
    rmpath([mfilepath, '/mcmc'])
    rmpath([mfilepath, '/solvers'])
    rmpath([mfilepath, '/qoi'])
end


