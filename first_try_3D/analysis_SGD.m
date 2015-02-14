analysis_note_consistency
%%
theta_now = [2,10,abs(rand(1)), abs(rand(1)), abs(rand(1))];
num_iter = 2000;
energyhistory = energy

eta = 0.001;
for k = 1: num_iter
    tic,
    [deriv, energy] = analysis_snap_deriv_tauleap(init, theta_now, tend, ...
        deltat, sigW, timesample, snapshots, N);
    theta_now = theta_now - eta*deriv';
    theta_now(1) =2;
    theta_now(2) =10;
    theta_now = max(theta_now, 0);
    display(['iteration', num2str(num_iter), ' complete', ' theta = ', num2str(theta_now)]);
    display(['energy = ', num2str(sum(energy))]);
    
    if(sum(energy) < 4)
        eta = 0.0005;
    end
    
    toc
end