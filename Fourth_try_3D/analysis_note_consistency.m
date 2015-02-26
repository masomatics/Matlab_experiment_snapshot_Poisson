
%%
init = [0;0;0];
theta = [2 10 0.6, 0.2,0,0.1,1];
%timesample = [5, 10, 12, 15];
num_parameters = length(theta);
%timesample=[3, 5,7,10];
timesample = [5];
tend = timesample(end);
%% First two part checks if the snapshot are indeed appropriate. 
tic,
N = 2000;
max_num_jumps = 10000;
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
    scatter3(snapshots(1,:,plotframe), snapshots(2,:,plotframe),snapshots(3,:,plotframe))
    xlabel('M')
    ylabel('P')
    zlabel('D')
end

%%
%close all;
tic,
N = 10000;
deltat = 0.1;
[euler_mean, deriv] = analysis_derivative_tauleap(init, theta, tend, ...
    deltat, timesample, snapshots, N);

toc


Gillespie_mean =squeeze(mean(snapshots,2))
euler_mean

%% Try taking derivative.  sigW = [a,b,c] must be chosen appropriately.
N = 3;
sigW = [0.5;2; 1] ;
tic,
[deriv, energy] =analysis_snap_deriv_tauleap(init, theta, tend, ...
    deltat, sigW, timesample, snapshots, N)
toc

theta_now = rand(1,num_parameters)

tic,
[deriv_now, energy_now] =analysis_snap_deriv_tauleap(init, theta_now, tend, ...
    deltat, sigW, timesample, snapshots, N)
toc

[euler_mean_now, deriv] = analysis_derivative_tauleap(init, theta_now, tend, ...
    deltat, timesample, snapshots, N);

%% KL 
%[deriv_now_KL, energy_now_KL] =analysis_snap_deriv_tauleapKL(init, theta, tend, ...
%    deltat, sigW, timesample, snapshots, N); 

