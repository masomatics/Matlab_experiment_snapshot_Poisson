theta_now_try = theta_now
eta = 0.001;
for k = 1: num_iter
    tic,
    [deriv, energy] = analysis_snap_deriv_tauleap(init, theta_now_try, tend, ...
        deltat, sigW, timesample, snapshots, N);
    theta_now_try = theta_now_try - eta*deriv';
 %   theta_now(1) =2; 
 deriv
    theta_now_try = max(theta_now_try, 0.0001);
    display(['iteration', num2str(k), ' complete', ' theta = ', num2str(theta_now_try)]);
    display(['energy = ', num2str(sum(energy))]);
    
    energyhistory_try(k) = sum(energy) ;
    thetahistory_try(k,:) = theta_now_try;
    
%    if(sum(energy) < 1) 0.33 for good parameter set2 
            if(sum(energy) < 0.38)

        eta = 0.0001;
    end
    
    toc
end