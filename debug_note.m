%% Generate data. This worked.
init = [0;0];
theta = [2, 10, 1/4 ,1];
timesample = [5, 11, 15, 30];

%%
N = 2000;
max_num_jumps = 5000;
rnsource1 = rand(N,max_num_jumps);
rnsource2 = rand(N,max_num_jumps);
tend = 31;


snapshots = analysis_data_generation_Gillespie(init, theta, tend, ...
    timesample, rnsource1, rnsource2, N);

mean(snapshots,2)


%%
tic,
N = 10000;
deriv = debug_analysis_derivative_tauleap(init, theta, 30, ...
    0.01, timesample, snapshots, N)

toc