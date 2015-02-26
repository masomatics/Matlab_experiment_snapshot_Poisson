parameters = load('good_parameters2.mat');
snapshots = load('good_snapshots2.mat');

parameters = parameters.s

theta= parameters.theta
timesample = parameters.timesample
tend = parameters.tend
sigW = parameters.sigW
deltat = parameters.deltat
N = parameters.N
init = parameters.init

snapshots = snapshots.sss;
