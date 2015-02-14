%analysis_data_generation_Gillespie
%
% This code computes the sensitivity of the model Anderson's Dimerization 
% Model
%
% rxn1 0 -> M 
% rxn2 M -> M + P 
% rxn3 M -> 0
% rxn4 M+ P -> M
% rxn5 P -> D
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

function [derivative, energy] = analysis_snap_deriv_tauleap(init, theta, tend, ...
    deltat, sigW, timesample, snapshots, N)

    include= [];

    %Preset variables
    [num_species, num_particles, num_slices] = size(snapshots);
    num_parameters = length(theta);
    num_species = length(init); 
    rxn_matrix = [1     0       -1      0       0; 
                  0     1       0       -1      0;
                  0     0       0       0       1];  
    [num_species, num_rxns] = size(rxn_matrix);
    compress_snap_wgts = 1/num_particles * ones(num_particles, num_slices);    
    num_leaps = ceil(tend/deltat);
    partial_dim = [1,2 ; 2,3; 1,3];  
     
    [num_partial_snaps,num_partial_species] = size(partial_dim); % one[1,2], two[2,3], three[1,3].
    
    %Target variables: This will be used three times. 
    derivative = zeros(num_parameters, 1);
    tilde_p_ymk = zeros(num_particles);
    p_ymkj = zeros(num_particles, N);
    dp_jr = zeros(N,num_particles);
    dEP_kr= zeros(N, num_particles);     
    
    
    %Temporary variables    
    ycopy = zeros(num_partial_species, num_particles, N); 
    xcopy = zeros(num_partial_species, num_particles, N);    
    scale_p_ymkj = zeros(1,N) ;
    scale_p_ymkj_mat = zeros(num_particles,N);
    rxn_rate = zeros(num_parameters, N);
    rxn_cnt = zeros(1,num_rxns); %number of rxns per channel that occurs in deltat
    sigW_copy = zeros(num_partial_species, num_particles, N);
    
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
        
        %Update rate, generate Poisson
        rxn_rate = compute_rate(dat_now, theta); 
        rxn_cnt = poissrnd(rxn_rate*deltat);
        dat_now = dat_now + rxn_matrix * rxn_cnt;
        dat_now = max(dat_now, 0);
               
        %compute DP,  N * num_parameters
        deriv_loglike = deriv_loglike + ...
            (rxn_cnt - deltat * rxn_rate).*repmat(1./theta',[1,N]);
        time_now = time_now + deltat;
        if(abs(time_now - timesample(slice_index))< deltat/10)
            
            includeon = 1;
            %derivative to be stored for this snapshot. 
            derivative_now = zeros(1,num_parameters);  
            
            for(alpha = 1: num_partial_snaps)
             ycopy = repmat(snapshots(partial_dim(alpha,:),:,slice_index), [1,1,N]);
             xcopy = permute(repmat(dat_now(partial_dim(alpha,:),:),[1,1,num_particles])...
                ,[1,3,2]);
             %Each dimension shall have different observation noise. 
             sigW_copy = repmat(sigW(partial_dim(alpha,:)),[1, num_particles, N]);
             distance_pair = min(squeeze(sum((xcopy - ycopy).^2.*(1./sigW_copy.^2),1)), 300);
             p_ymkj = exp(-(distance_pair));
             scale_p_ymkj = 1./sum(p_ymkj,1);
             scale_p_ymkj_mat = repmat(scale_p_ymkj,[num_particles,1]) ;
             p_ymkj = p_ymkj .* scale_p_ymkj_mat;
             tilde_p_ymk = mean(p_ymkj,2);
             dEP_kr = 1/N *p_ymkj * deriv_loglike' ;
             derivative_now = derivative_now + sign(tilde_p_ymk - ...
                compress_snap_wgts(:,slice_index))' * (dEP_kr); 
            energy(slice_index) = energy(slice_index) + sum(abs(tilde_p_ymk -compress_snap_wgts(:,slice_index)));

            %Check if the gradient is any meaningful (tilde p is not uniform) 
            includeon = includeon *(max(p_ymkj(:)) > min(p_ymkj(:)));
                         
            end   
            derivative_alltime(:,slice_index) = derivative_now; %[1,5] signed to [5,1] is OK in matlab
            
            
            if(includeon)
                include = [include, slice_index];
            end            
            slice_index = slice_index  + 1;
       
            
        end %end of Frame 
        
        derivative = sum(derivative_alltime(:,include),2);
    end % end of Tauleap loop 

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
    rate(5,:) = theta(5)*x(2,:);
end 



