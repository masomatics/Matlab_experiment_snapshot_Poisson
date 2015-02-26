parameters = load('good_parameters.mat');
snapshots = load('good_snapshots.mat');

parameters = parameters.s

theta= parameters.theta
timesample = parameters.timesample
tend = parameters.tend
sigW = parameters.sigW
deltat = parameters.deltat
N = parameters.N
init = parameters.init

snapshots = snapshots.sss;
