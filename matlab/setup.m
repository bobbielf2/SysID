function setup(v)

if nargin == 0
    mfilepath=fileparts(mfilename('fullpath'));
    addpath([mfilepath, '/mcmc'])
    addpath([mfilepath, '/solvers'])
else
    mfilepath=fileparts(mfilename('fullpath'));
    rmpath([mfilepath, '/mcmc'])
    rmpath([mfilepath, '/solvers'])
end


