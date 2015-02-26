analysis_note_consistency
%%
N = 5000
theta_now = rand(1,num_parameters);
%%
%theta_now(2) = 6;
%theta_now = [1.0257      6.1349           0           0]
%theta_now = [1.2538    6.5528         0.0001         0.0001];
%theta_now = [2, 10, 0.0001,0.0001];

num_iter = 50000;
energyhistory_M = zeros(1, num_iter);
thetahistory_M = zeros(num_iter, length(theta_now));
thetahistory_M(1,:) = theta_now;
theta_now(2) =10

display([' REal theta = ', num2str(theta)])
eta = 0.00001;
tic,
for k = 1: num_iter
    
    [deriv, energy] = analysis_snap_deriv_tauleap_mean(init, theta_now, tend, ...
        deltat, sigW, timesample, snapshots, N);
    theta_now = theta_now - eta*deriv';
 %   theta_now(1) =2; 
% deriv
 %   theta_now(2) =10;
    theta_now = max(theta_now, 0.0001);
    if(mod(k, 100) == 0)
    display(['iteration', num2str(k), ' complete', ' theta = ', num2str(theta_now)]);
        display(['energy = ', num2str(sum(energy))]);
        toc
        tic,
    end
    energyhistory_M(k) = sum(energy) ;
    thetahistory_M(k,:) = theta_now;
    if(sum(energy) < 30)
        eta = 0.0005;
    end
    
%    if(sum(energy) < 1) 0.33 for good parameter set2 

    
    
end
toc

%eee = energyhistory(1 : k);
%ttt = thetahistory(1:k, :);
%sss = snapshots;

% value1 = theta;
% value2 = timesample; 
% value3 = tend; 
% value4 = sigW;
% value5 = deltat;
% value6 = N; 
% value7 = init;
% field1 = 'theta';
% field2 = 'timesample'; 
% field3 = 'tend';
% field4 = 'sigW';
% field5 = 'deltat'; 
% field6 = 'N';
% field7 = 'init'; 
% 
% s = struct(field1, value1,field2, value2,field3, value3,field4, value4, ...
% field5, value5, field6, value6, field7, value7);  
% 
% save('good_parameters2.mat', 's')

