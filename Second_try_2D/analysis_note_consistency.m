
%%
init = [0;0];
theta = [2 10 0.6, 0.2];
%timesample = [5, 10, 12, 15];
timesample=[4, 5,7,10];
tend = timesample(end);
%% First two part checks if the snapshot are indeed appropriate. 
tic,
N = 2000;
max_num_jumps = 5000;
rnsource1 = rand(N,max_num_jumps);
rnsource2 = rand(N,max_num_jumps);
snapshots = zeros(length(init), N, length(timesample));
%Something is wrong with the slice_index loop. As annoying as it is,
%I will do what the Biologists actually do. Run a new batch of simulation
%for each timepoint. 
%snapshots = debug_analysis_data_generation_Gillespie(init, theta, tend, ...
%    timesample, rnsource1, rnsource2, N);

Gillespie_mean = zeros(length(init),length(timesample));
for(tt = 1:length(timesample))

    snapshot = analysis_data_generation_Gillespie(init, theta, tend, ...
    timesample(tt), rnsource1, rnsource2, N);

    Gillespie_mean(:,tt) =squeeze(mean(snapshot,2));

    snapshots(:,:,tt) = snapshot;
    display(['snapshot at t=', num2str(timesample(tt)), 'complete']);
end 
toc

figure(33)
for plotframe = 1:length(timesample)
    subplot(2,ceil(length(timesample)/2), plotframe)
    scatter(snapshots(1,:,plotframe), snapshots(2,:,plotframe))
end

%%
close all;
tic,
N = 10000;
deltat = 0.1;
[euler_mean, deriv] = analysis_derivative_tauleap(init, theta, tend, ...
    deltat, timesample, snapshots, N);

toc


Gillespie_mean =squeeze(mean(snapshots,2))
euler_mean

%% Try taking derivative.  sigW = [a,b,c] must be chosen appropriately.
N = 5000;
sigW = [0.5;2] ;
tic,
[deriv, energy] =analysis_snap_deriv_tauleap(init, theta, tend, ...
    deltat, sigW, timesample, snapshots, N)
toc


