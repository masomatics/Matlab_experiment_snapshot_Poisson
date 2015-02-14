initx = [1;0;0];
theta = [2 10 1/4, 0.2,0.1];
tend = 20;
sigV = 2;
sigW = 2;

%num_timepts = 2500;
num_timepts = 2500;
%Ntry = 10000;
Ntry = 4000
tic,
rnsource = randn([3, Ntry, num_timepts]);
toc
tic,
%[timepts,datapts] = datagen_mass2(initx, tend, theta, sigV, num_timepts, rnsource, Ntry);

[timepts,datapts] = and_CFD_datagen_mass_diffusion(initx, tend, theta, sigV, num_timepts, rnsource, Ntry);

%%

timesample = [5:5:20];
timeindex = find(ismember(timepts, timesample));

rnsource2 = randn([3, Ntry, length(timesample)]);
snapshots = datapts(:, :, timeindex) + rnsource2;

%%
Ntry = 4000;
%tic,
rnsource = randn([3, Ntry, num_timepts]);
%toc
theta0 = [1,5,1,1,1];
theta0 = theta
tic,
[tilde_pys, deriv] = and_CFD_datagen_mass_derivStat...
                        (initx, tend, theta0, sigV, sigW, num_timepts, rnsource, snapshots, timesample, Ntry);
toc
%%
figure(99)

for frame_index = 1:length(timesample)
    subplot(2,ceil((length(timesample)-1)/2),frame_index)
 scatter3(snapshots(1,:,frame_index),snapshots(2,:,frame_index), snapshots(3,:,frame_index), 8, tilde_pys(frame_index,:))
 title(['time at ', num2str(timesample(frame_index))]);
end
%tic,
%deriv = and_CFD_datagen_mass_derivStat2...
%                        (initx, tend, theta0, sigV, sigW, num_timepts, rnsource, snapshots, timesample, Ntry)
%toc



%%
close all;
Ntry = 4000;
%tic,
rnsource = randn([3, Ntry, num_timepts]);
%toc

eta = 0.01; 
num_iter = 1500
%theta0 = [1,5,1,1,1];


thetahistory = zeros(1+num_iter,5);
thetahistory(1,:) = theta0;

writerObj = VideoWriter('snapshots_sequence.avi');
writerObj.FrameRate= 1;
set(gcf,'Position',get(0,'ScreenSize'))
open(writerObj)
for iter = 1:num_iter
    display(['initiating the ', num2str(iter), 'th iterations...']);
    tic,
    [tilde_pys, deriv] = and_CFD_datagen_mass_derivStat...
                        (initx, tend, theta0, sigV, sigW, num_timepts, rnsource, snapshots, timesample, Ntry);
    display([num2str(iter), 'th iteration : ' , num2str(deriv')]);
    
    for frame_index = 1:length(timesample)
    subplot(2,ceil((length(timesample)-1)/2),frame_index)
    scatter3(snapshots(1,:,frame_index),snapshots(2,:,frame_index), snapshots(3,:,frame_index),8, tilde_pys(frame_index,:));
    title(['t =', num2str(timesample(frame_index)), ' Iteration ', num2str(iter), ' theta0=', num2str(theta0)]);
    end

    M(iter) = getframe(gcf);
    writeVideo(writerObj,M(iter));
    %I am trying to increase the energy, so I will take positive gradient
    theta0 = max(theta0 + eta*deriv',0);
    theta0(5) = min(theta0(5),1);
        display([num2str(iter), 'th iteration complete:  new theta0 is =' , num2str(theta0)]);
        thetahistory(iter+1,:) = theta0;
        
    display(['Completing', num2str(iter), 'th iterations...']);
    toc
end 
close(writerObj);
save('thetahistory1.mat', 'thetahistory')



