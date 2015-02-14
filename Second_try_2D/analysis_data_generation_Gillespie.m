%analysis_data_generation_Gillespie
%
% I must admit that this is a very slow implementation.
%
% This is a Gillespie implmenation for the model 
% rxn1 0 -> M 
% rxn2 M -> M + P 
% rxn3 M -> 0
% rxn4 M+ P -> M
%
%Input
%
% int init[num_species,1] :initial value
% double theta[num_parameters] : parameters
% rnsource1[NN]:   uniform random variable, number of expected jumps  * N  
% rnsource2[NN]:   uniform random variable, number of expected jumps  * N  
% timeample    :array of sampling time, listed in increasing order.
% double tend : terminal time
% int N : number of samples
%
%Output
%
% snapshots
%
%




function snapshots = analysis_data_generation_Gillespie(init, theta, tend, ...
    timesample, rnsource1, rnsource2, N)

    %Preset variables
    num_slices = length(timesample);  
    num_species = length(init);
    
    %Initialization 
    t_sub0 = 0;
    init_now = repmat(init, [1,N]);
    %Target variables
    snapshots = zeros(num_species, N, num_slices);

    %loop of filling in the data
    for(slice_index = 1 : num_slices) 

        t_subtend = timesample(slice_index);

        for(j = 1:N)
            snapshots(:,j,slice_index) = singlepath(init_now(:,j), theta,...
                                (t_subtend- t_sub0),rnsource1(j,:), rnsource2(j,:));
        end 
        t_sub0 = t_subtend;
        init_now = squeeze(snapshots(:,:,slice_index));  
    end 

end


function data_t_subend = singlepath(init, theta, t_subend, rnsource1, rnsource2) 

%Preset variables
num_parameters = length(theta);
num_species = length(init); 
rxn_matrix = [1     0       -1      0       ; 
                  0     1       0       -1      ];
max_num_jumps = 10000;
%Initialization  
dat_now = init;
time_now = 0; 
rand_index = 1; 

%Temporary variables
rxn_choice = 0;
deltat = 0; 
rate = zeros(1, num_parameters);
cumrate = zeros(1, num_parameters-1);
%Debug 

%Gillespie Loop
for(t = 1:max_num_jumps)
    
    %Compute rate and generate the necessary random variables
    [cumrate, rate] = compute_rate(dat_now, theta);
    rxn_choice = 1 + sum(cumrate < (rnsource1(rand_index)*sum(rate))); 
    deltat = -log(rnsource2(rand_index))/sum(rate);
    time_now = deltat + time_now;
    if(time_now >t_subend)
        break;
    end 
    %Update the  datmat, random index, and time
    dat_now = dat_now + rxn_matrix(:,rxn_choice);
    rand_index = rand_index +1;
       

end %end of the time loop

data_t_subend = dat_now;
end 

%takes in the current state x and returns the rate at the point 
function [cumrate, rate] = compute_rate(x, theta)
    rate(1) = theta(1);
    rate(2) = theta(2)*x(1);
    rate(3) = theta(3)*x(1);
    rate(4) = theta(4)*x(2)*x(1);

    cumrate(1) = rate(1);
    cumrate(2) = cumrate(1) +rate(2);
    cumrate(3) = cumrate(2) +rate(3);
end 


