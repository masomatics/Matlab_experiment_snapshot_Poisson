%% Generate data. This worked.
init = [0;0];
theta = [2, 10, 1/4,1];
timesample = [1,3,5, 10];
tend = timesample(end);
%%
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

    snapshot = debug_analysis_data_generation_Gillespie(init, theta, tend, ...
    timesample(tt), rnsource1, rnsource2, N);

    Gillespie_mean(:,tt) =squeeze(mean(snapshot,2));

    snapshots(:,:,tt) = snapshot;
    display(['snapshot at t=', num2str(timesample(tt)), 'complete']);
end 





%%
tic,
N = 10000;
[meandat, deriv] = debug_analysis_derivative_tauleap(init, theta, tend, ...
    0.1, timesample, snapshots, N);

toc

%% 
Xt_real = repmat(init, [1, length(timesample)]);
for r = 1:length(timesample)
    Xt_real(:,r) = debug_sanity_check(init,theta, timesample(r));
end 
%Xt_real

%% 
Gillespie_mean =squeeze(mean(snapshots,2))
Xt_real
meandat