%analysis_snap_deriv_tauleap
%
% This code computes the sensitivity of the model Anderson's Dimerization 
% Model
%
% rxn1 0 -> M 
% rxn2 M -> M + P 
% rxn3 M -> 0
% rxn4 M+ P -> M
%using  Girsanov's method.  This is the implementation of SnapNote 6.8. The
%three chosen pairs of dimensions are [1,2], [2,3], [1,3].  
%
%Input
%
% int init[num_species,1] :initial value
% double theta[num_parameters] : parameters
% rnsource1[NN]:   uniform random variable, number of expected jumps  * N  
% rnsource2[NN]:   uniform random variable, number of expected jumps  * N  
% timeample    :   array of sampling time, listed in increasing order.
% deltat       :   deltat for Tauleap. 
% sigW         :   THREE dimensional COLUMN vector, independent noises. 
%
%Output
%
% debug_derivative : deriv X_2(timesample)   num_parameters * num_slices
%                                            matrix
%
%
%

function [derivative, energy] = analysis_snap_deriv_tauleap_mean(init, theta, tend, ...
    deltat, sigW, timesample, snapshots, N)

   % include= [];

    %Preset variables
    [num_species, num_particles, num_slices] = size(snapshots);
    num_parameters = length(theta);
    num_species = length(init); 
    rxn_matrix = [1     0       -1      0       ; 
                  0     1       0       -1      ];
    [num_species, num_rxns] = size(rxn_matrix);
    compress_snap_wgts = 1/num_particles * ones(num_particles, num_slices);    
    num_leaps = ceil(tend/deltat);     
        
    %Target variables
    derivative = zeros(num_parameters, 1);    
        
    %Temporary variables    
    rxn_rate = zeros(num_parameters, N);
    rxn_cnt = zeros(1,num_rxns); %number of rxns per channel that occurs in deltat
    
    %Initialization 
    dat_now = repmat(init,[1,N]);
    time_now = 0; 
    slice_index = 1;
    deriv_loglike = zeros(num_parameters,N);
    derivative_alltime = zeros(num_parameters, num_slices); 
    mean_snaps = mean(snapshots, 2);
    pre_derivative = zeros(num_species, num_parameters); 

    %Debug 
    energy = zeros(1, num_slices);  
    deriv_mean = NaN; 

    %Monte Carlo Simulation 
    for(h = 1:(num_leaps)) 
        
        %Update rate, generate Poisson
        rxn_rate = compute_rate(dat_now, theta); 
        rxn_cnt = poissrnd(rxn_rate*deltat);
        dat_now = dat_now + rxn_matrix * rxn_cnt;
        dat_now = max(dat_now, 0);
               
        %compute DP,  N * num_parameters
        deriv_loglike = deriv_loglike + ...
            (rxn_cnt - deltat * rxn_rate).*repmat(1./theta',[1,N]);
        time_now = time_now + deltat;

        %Compute derivative
        if(abs(time_now - timesample(slice_index))< deltat/10)
            
            derivative_now = zeros(1,num_parameters);       
            datmean = mean(dat_now, 2);
            pre_derivative = 1/N*dat_now * deriv_loglike';
            derivative_now = sign(datmean - mean(snapshots(:,:,slice_index),2))'...
                *pre_derivative;
            derivative_alltime(:,slice_index) = derivative_now; %[1,5] signed to [5,1] is OK in matlab                   

            energy(slice_index) = sum(abs(datmean - mean(snapshots(:,:,slice_index),2))); 
            slice_index = slice_index  + 1;

        end %end of Frame 
        
    end % end of Tauleap loop 
 derivative = sum(derivative_alltime,2);
end %End of the main function

%takes in the current state x and returns the rate at the point 
function rate = compute_rate(x, theta)
    
    num_parameters = length(theta);
    [num_species , N ] = size(x);
    rate= zeros(num_parameters, N);
    rate(1,:) = theta(1);
    rate(2,:) = theta(2)*x(1,:);
    rate(3,:) = theta(3)*x(1,:);
    rate(4,:) = theta(4).*x(2,:).*x(1,:);
end 



