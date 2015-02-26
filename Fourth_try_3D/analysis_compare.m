tic,
theta_approx = [1.5759      9.8103       0.461     0.25197    0.012318    0.084743     0.84547]
N = 2000;
max_num_jumps = 10000;
rnsource1_exhibit = rand(N,max_num_jumps);
rnsource2_exhibit = rand(N,max_num_jumps);
snapshots_real = zeros(length(init), N, length(timesample));
%Something is wrong with the slice_index loop. As annoying as it is,
%I will do what the Biologists actually do. Run a new batch of simulation
%for each timepoint. 
%snapshots = debug_analysis_data_generation_Gillespie(init, theta, tend, ...
%    timesample, rnsource1, rnsource2, N);

Gillespie_mean_real = zeros(length(init),length(timesample));
for(tt = 1:length(timesample))

    snapshot = analysis_data_generation_Gillespie(init, theta, tend, ...
    timesample(tt), rnsource1_exhibit, rnsource2_exhibit, N);

    Gillespie_mean_real(:,tt) =squeeze(mean(snapshot,2));

    snapshots_real(:,:,tt) = snapshot;
    display(['snapshot at t=', num2str(timesample(tt)), 'complete']);
end 
toc

snapshots_approx = zeros(length(init), N, length(timesample));
Gillespie_mean_approx = zeros(length(init),length(timesample));
for(tt = 1:length(timesample))

    snapshot = analysis_data_generation_Gillespie(init, theta_approx, tend, ...
    timesample(tt), rnsource1_exhibit, rnsource2_exhibit, N);

    Gillespie_mean_approx(:,tt) =squeeze(mean(snapshot,2));

    snapshots_approx(:,:,tt) = snapshot;
    display(['snapshot at t=', num2str(timesample(tt)), 'complete']);
end 
toc


figure(1)
scatter3(snapshots_real(1,:,plotframe), snapshots_real(2,:,plotframe),snapshots_real(3,:,plotframe))
hold on;
scatter3(snapshots_approx(1,:,plotframe), snapshots_approx(2,:,plotframe),snapshots_approx(3,:,plotframe))
    xlabel('M')
    ylabel('P')
    zlabel('D')
figure(2)
subplot(2,2,1) 
scatter(snapshots_real(1,:,plotframe), snapshots_real(2,:,plotframe))
hold on;
scatter(snapshots_approx(1,:,plotframe), snapshots_approx(2,:,plotframe))
    xlabel('M')
    ylabel('P')

subplot(2,2,2) 
scatter(snapshots_real(2,:,plotframe), snapshots_real(3,:,plotframe))
hold on;
scatter(snapshots_approx(2,:,plotframe), snapshots_approx(3,:,plotframe))
    xlabel('P')
    ylabel('D')

subplot(2,2,3) 
scatter(snapshots_real(1,:,plotframe), snapshots_real(3,:,plotframe))
hold on;
scatter(snapshots_approx(1,:,plotframe), snapshots_approx(3,:,plotframe))
    xlabel('M')
    ylabel('D')