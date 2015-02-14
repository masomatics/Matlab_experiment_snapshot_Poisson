initx = [1 ; 0; 0];
theta = [2, 10, 1/4, 0.01, 5];
tend = 5;
sigV = 2;
sigW = 3;
num_timepts = 2500;
Ntry_sample = 2000;
tic,
rnsource_sample_sys = randn([2, Ntry_sample, num_timepts]);
rnsource_sample_obs = randn([2, Ntry_sample, num_timepts]);

toc
tic,

timesample= sort(tend*(rand([1,5])))
%Create Snapshots Data
tic,
[timemat, datmatx] =and_CFD_datagen_mass_diffusion(initx, tend, theta, sigV,...
                                num_timepts, rnsource_sample_sys, Ntry_sample);
toc 
datmaty = max(datmatx + sigW * rnsource_sample_obs,0);



snapshot = datmaty(:, :, timesample);  

%num_timepts = 10000;
%Ntry_sim = 5000;
%tic,
%rnsource_sim = randn([2, Ntry, num_timepts]);
%toc
%tic,