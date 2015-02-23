%% TRUE
N = 2000;
max_num_jumps = 5000;

%rnsource1_test = rand(N,max_num_jumps);
%rnsource2_test = rand(N,max_num_jumps);

for(tt = 1:length(timesample))

    snapshot_true = analysis_data_generation_Gillespie(init, theta, tend, ...
    timesample(tt), rnsource1, rnsource2, N);


    snapshots_true(:,:,tt) = snapshot_true;
    display(['snapshot at t=', num2str(timesample(tt)), 'complete']);
end 





%% APPROX1
theta_init = thetahistory(1,:);

for(tt = 1:length(timesample))

    snapshot_init = analysis_data_generation_Gillespie(init, theta_init, tend, ...
    timesample(tt), rnsource1, rnsource2, N);

    snapshots_init(:,:,tt) = snapshot_init;
    display(['snapshot at t=', num2str(timesample(tt)), 'complete']);
end 


figure(98)
for plotframe = 1:length(timesample)
    subplot(2,ceil(length(timesample)/2), plotframe)
    scatter(snapshots(1,:,plotframe), snapshots(2,:,plotframe), 'b')
    hold on;
    scatter(snapshots_init(1,:,plotframe), snapshots_init(2,:,plotframe), 'r', 'x')
    legend('true parameter', 'approximate')
    title(['t = ', num2str(timesample(plotframe))])
end

%%
%% APPROX2
theta_init = thetahistory(10292,:);

for(tt = 1:length(timesample))

    snapshot_init = analysis_data_generation_Gillespie(init, theta_init, tend, ...
    timesample(tt), rnsource1, rnsource2, N);

    snapshots_init(:,:,tt) = snapshot_init;
    display(['snapshot at t=', num2str(timesample(tt)), 'complete']);
end 


figure(99)
for plotframe = 1:length(timesample)
    subplot(2,ceil(length(timesample)/2), plotframe)
    scatter(snapshots(1,:,plotframe), snapshots(2,:,plotframe), 'b')
    hold on;
    scatter(snapshots_init(1,:,plotframe), snapshots_init(2,:,plotframe), 'r', 'x')
    legend('true parameter', 'approximate')
    title(['t = ', num2str(timesample(plotframe))])
end

hold off;

