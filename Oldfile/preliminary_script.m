initx = [1;0;0];
theta = [2 10 1/4, 0.2,0.1];
theta = [1 5 1 1 2]
theta =  [3.0377    7.1908         0         0    0.2078]
tend = 20;
sigV = 2;
sigW = 2;


num_timepts = 2500;
%Ntry = 10;
Ntry = 400
tic,
rnsource = randn([3, Ntry, num_timepts]);
toc
tic,

[timepts,datapts] = and_CFD_datagen_mass_diffusion(initx, tend, theta, sigV, num_timepts, rnsource, Ntry);

display('mean is') 
[mean(datapts(1,:)) ; mean(datapts(2,:)) ; mean(datapts(3,:)) ]

