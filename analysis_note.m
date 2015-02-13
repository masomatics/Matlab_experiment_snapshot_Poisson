%%
init = [1;0;0];
theta = [2 10 1/4, 0.2,0.1];
timesample = [1,4,10,15];
tend = timesample(end);
%%
tic,
N = 10000;
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

%%
tic,
N = 10000;
[euler_mean, deriv] = analysis_derivative_tauleap(init, theta, tend, ...
    0.1, timesample, snapshots, N);

toc

Gillespie_mean =squeeze(mean(snapshots,2))
euler_mean