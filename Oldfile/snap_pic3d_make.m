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
timesample=[30, 500, 1000, 2000];
snapshots = datapts(:,:,timesample);
%%
figure(99)

for frame_index = 1:length(timesample)
    subplot(2,ceil((length(timesample)-1)/2),frame_index)
 scatter3(snapshots(1,:,frame_index),snapshots(2,:,frame_index), snapshots(3,:,frame_index))
 xlabel('x_1', 'FontSize', 15)
 ylabel('x_2', 'FontSize', 15)
 zlabel('x_3', 'FontSize', 15)
 title(['t= ', num2str(timesample(frame_index))], 'FontSize', 15);
end
%tic,


%%
figure(2)
for frame_index = 1:length(timesample)
    subplot(2,ceil((length(timesample)-1)/2),frame_index)
 scatter(snapshots(1,:,frame_index),snapshots(2,:,frame_index))
 %xlabel('x_1', 'FontSize', 15)
 %ylabel('x_2', 'FontSize', 15)
 title(['t= ', num2str(timesample(frame_index))], 'FontSize', 15);
 
 
 set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
 
end