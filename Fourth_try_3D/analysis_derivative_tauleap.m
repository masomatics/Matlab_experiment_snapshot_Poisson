%analysis_data_generation_Gillespie
%
% This code computes the sensitivity of the model Anderson's Dimerization 
% Model
%
% rxn1 0 -> M 
% rxn2 M -> M + P 
% rxn3 M -> 0
% rxn4 M+ P -> M
% rxn5 P  -> 0 
% rxn6 2P -> D
% rxn7 D  -> 0
%using  Girsanov's method.  This is a Debug version and only compute the
% sensitivity with respect to the mean of the second coordinate. 
%
%Input
%
% int init[num_species,1] :initial value
% double theta[num_parameters] : parameters
% rnsource1[NN]:   uniform random variable, number of expected jumps  * N  
% rnsource2[NN]:   uniform random variable, number of expected jumps  * N  
% timeample    :array of sampling time, listed in increasing order.
%
%Output
%
% debug_derivative : deriv X_2(timesample)   num_parameters * num_slices
%                                            matrix
%
%
%

function [meandat, debug_derivative] = analysis_derivative_tauleap(init, theta, tend, ...
    deltat, timesample, snapshots, N)

    include= [];

    %Preset variables
    [num_species, num_particles, num_slices] = size(snapshots);
    num_parameters = length(theta);
    num_species = length(init); 
    rxn_matrix = [1     0       -1      0       0      0    0; 
                  0     1       0       -1    -1  -2    0;
                  0     0       0        0     0   1    -1];
    [num_species, num_rxns] = size(rxn_matrix);
    compress_snap_wgts = 1/num_particles * ones(num_particles, num_slices);    
    num_leaps = ceil(tend/deltat);
 
    %Target variables
    derivative = zeros(num_parameters, 1);
    tilde_p_ymk = zeros(num_particles);
    p_ymkj = zeros(num_particles, N);
    dp_jr = zeros(N,num_particles);
    dEP_kr= zeros(N, num_particles);     
    
    %Temporary variables    
    ycopy = zeros(num_species, num_particles, N); 
    xcopy = zeros(num_species, num_particles, N);    
    scale_p_ymkj = zeros(1,N) ;
    scale_p_ymkj_mat = zeros(num_particles,N);
    rxn_rate = zeros(num_parameters, N);
    rxn_cnt = zeros(1,num_rxns); %number of rxns per channel that occurs in deltat
    
    %Initialization 
    dat_now = repmat(init,[1,N]);
    time_now = 0; 
    slice_index = 1;
    deriv_loglike = zeros(num_parameters,N);
    derivative_alltime = zeros(num_parameters, num_slices); 

    %Debug 
    energy = zeros(1, num_slices);  
    deriv_mean = NaN; 

%    debug_x2_copy = zeros(N, num_parameters);
%    debug_dEP_2 = zeros(N, num_parameters); 
%    debug_derivative = zeros(num_parameters, num_slices);
    
    %Monte Carlo Simulation 
    for(h = 1:(num_leaps)) 
        rxn_rate = compute_rate(dat_now, theta); 
        rxn_cnt = poissrnd(rxn_rate*deltat);
        dat_now = dat_now + rxn_matrix * rxn_cnt;
        dat_now = max(dat_now, 0);
        
%         deriv_loglike(1,:)= deriv_loglike(1,:) + ...
%             (rxn_cnt(1,:)- deltat *rxn_rate(1,:))/theta(1);
%         deriv_loglike(2,:)= deriv_loglike(2,:) + ...
%             (rxn_cnt(2,:)- deltat *rxn_rate(2,:))/theta(2);
%         deriv_loglike(3,:)= deriv_loglike(3,:) + ...
%             (rxn_cnt(3,:)- deltat *rxn_rate(3,:))/theta(3);
%         deriv_loglike(4,:)= deriv_loglike(4,:) + ...
%             (rxn_cnt(4,:)- deltat *rxn_rate(4,:))/theta(4);       
        
        deriv_loglike = deriv_loglike + ...
            (rxn_cnt - deltat * rxn_rate).*repmat(1./theta',[1,N]);
        time_now = time_now + deltat;
        if(abs(time_now - timesample(slice_index))< deltat/10)
             time_now
             debug_dEP_2 = dat_now(2,:) * deriv_loglike'/N;
             debug_derivative(:,slice_index) = debug_dEP_2';
             
             meandat(:,slice_index) = mean(dat_now, 2);             
             slice_index = slice_index  + 1;

        end 
        

    end 

end

%takes in the current state x and returns the rate at the point 
function rate = compute_rate(x, theta)
    
    num_parameters = length(theta);
    [num_species , N ] = size(x);
    rate= zeros(num_parameters, N);
    rate(1,:) = theta(1);
    rate(2,:) = theta(2)*x(1,:);
    rate(3,:) = theta(3)*x(1,:);
    rate(4,:) = theta(4).*x(2,:).*x(1,:);
    rate(5,:) = theta(5)*x(2,:);
    rate(6,:) = theta(6).*x(2,:).*x(2,:);
    rate(7,:) = theta(7)*x(3,:);
end 



