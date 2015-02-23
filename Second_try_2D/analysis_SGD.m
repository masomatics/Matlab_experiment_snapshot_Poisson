analysis_note_consistency
%%
N = 5000
theta_now = rand(1,4);
%%
%theta_now(2) = 6;
%theta_now = [1.0257      6.1349           0           0]
%theta_now = [1.2538    6.5528         0.0001         0.0001];
%theta_now = [2, 10, 0.0001,0.0001];

num_iter = 11000;
energyhistory = zeros(1, num_iter);
thetahistory = zeros(num_iter, length(theta_now));
thetahistory(1,:) = theta_now;

display([' REal theta = ', num2str(theta)])
eta = 0.001;
for k = 1: num_iter
    tic,
    [deriv, energy] = analysis_snap_deriv_tauleap(init, theta_now, tend, ...
        deltat, sigW, timesample, snapshots, N);
    theta_now = theta_now - eta*deriv';
 %   theta_now(1) =2; 
 deriv
    theta_now(2) =10;
    theta_now = max(theta_now, 0.0001);
    display(['iteration', num2str(k), ' complete', ' theta = ', num2str(theta_now)]);
    display(['energy = ', num2str(sum(energy))]);
    
    energyhistory(k) = sum(energy) ;
    thetahistory(k,:) = theta_now;
    if(sum(energy) < 2)
        eta = 0.0005;
    end
    
%    if(sum(energy) < 1)
            if(sum(energy) < 0.33)

        eta = 0.0001;
    end
    
    toc
end