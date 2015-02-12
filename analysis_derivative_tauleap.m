%analysis_data_generation_Gillespie
%
% This is meant to be run for three dimensional system with 
% snapshots being the 2D scatter plot of the three pairs of species. 
% The algorithm follows the SnapNote 6.8.  To change the system, one needs
% only to change the rate subroutine and the rxn_matrix.
%
%Input
%
% int init[num_species,1] :initial value
% double theta[num_parameters] : parameters
% rnsource1[NN]:   uniform random variable, number of expected jumps  * N  
% rnsource2[NN]:   uniform random variable, number of expected jumps  * N  
% timeample    :array of sampling time, listed in increasing order.
% snapshots    :num_species  * num_particles * num_slices 
%
%Output
%
% debug_derivative : deriv X_2(timesample)   num_parameters * num_slices
%                                            matrix
%
%
%

function debug_derivative = analysis_derivative_tauleap(init, theta, tend, ...
    deltat, timesample, snapshots, N)

    include= [];

    %Preset variables
    [num_species, num_particles, num_slices] = size(snapshots);
    snapshots1 = snapshots([1,2],:,:);
    snapshots2 = snapshots([2,3],:,:);
    snapshots3 = snapshots([1,3],:,:);
    num_parameters = length(theta);
    num_species = length(init); 
    rxn_matrix = [1 0 -1 0; 0 1 0 -1]; 
    [num_species, num_rxns] = size(rxn_matrix);
    compress_snap_wgts1 = 1/num_particles * ones(num_particles, num_slices);    
    compress_snap_wgts2 = 1/num_particles * ones(num_particles, num_slices);    
    compress_snap_wgts3 = 1/num_particles * ones(num_particles, num_slices);    

    num_leaps = ceil(tend/deltat);
    obs_dim = 2;  %Simultanesouly observed species
 
    %Target variables
    derivative = zeros(num_parameters, 1);    
    tilde_p_ymk1 = zeros(num_particles);
    tilde_p_ymk2 = zeros(num_particles);
    tilde_p_ymk3 = zeros(num_particles);
    p_ymkj1 = zeros(num_particles, N);
    p_ymkj2 = zeros(num_particles, N);
    p_ymkj3 = zeros(num_particles, N);
    dp_jr1 = zeros(N,num_particles);
    dp_jr2 = zeros(N,num_particles);
    dp_jr3 = zeros(N,num_particles);    
    dEP_kr1= zeros(N, num_particles);     
    dEP_kr2= zeros(N, num_particles);     
    dEP_kr3= zeros(N, num_particles);     
    
    %Temporary variables    
    ycopy = zeros(obs_dim, num_particles, N); 
    xcopy = zeros(obs_dim, num_particles, N);        
    scale_p_ymkj = zeros(1,N) ;
    scale_p_ymkj_mat = zeros(num_particles,N);
    rxn_rate = zeros(num_parameters, N);
    rxn_cnt = zeros(1,num_rxns); %number of rxns per channel that occurs in deltat
    derivative_now = zeros(num_parameters,1);
    
    %Initialization 
    data_now = repmat(init,[1,N]);
    time_now = 0; 
    slice_index = 1;
    deriv_loglike = zeros(num_parameters,N);
    derivative_alltime = zeros(num_parameters, num_slices); 

    %Debug 
    energy = zeros(1, num_slices);  
    deriv_mean = NaN; 

    %Monte Carlo Simulation 
    for(h = 1:(num_leaps+1)) 
        rxn_rate = compute_rate(data_now, theta); 
        rxn_cnt = poissrnd(rxn_rate*deltat);
        data_now = data_now + rxn_matrix * rxn_cnt;
        data_now = max(data_now, 0);
                     
        deriv_loglike = deriv_loglike + ...
            (rxn_cnt - deltat * rxn_rate).*repmat(1./theta',[1,N]);
        
        if(abs(time_now - timesample(slice_index))< deltat/10)
            
             derivative_now = zeros(num_parameters,1);

             %P1
             ycopy = repmat(snapshots1(:,:,snaptime_now), [1,1,N]);
             xcopy = permute(repmat(data_now([1,2],:),[1,1,num_particles])...
                ,[1,3,2]);
             distance_pair = min(squeeze(sum((xcopy - ycopy).^2,1)/sigW^2), 300);
             p_ymkj1 = exp(-(distance_pair));
             scale_p_ymkj = 1./sum(p_ymkj1,1);
             scale_p_ymkj_mat = repmat(scale_p_ymkj,[num_particles,1]) ;
             p_ymkj1 = p_ymkj1 .* scale_p_ymkj_mat;
             tilde_p_ymk1 = mean(p_ymkj1,2);
             dEP_kr1 = 1/N *p_ymkj1 * deriv_loglike' ;
             derivative_now = derivative_now + sign(tilde_p_ymk1 - ...
                compress_snap_wgts1(:,snaptime_now))' * (dEP_kr1);
             energy(snaptime_now) = energy(snaptime_now) + ...
                 sum(abs(tilde_p_ymk1 -compress_snap_wgts1(:,snaptime_now)));

             %P2
             ycopy = repmat(snapshots2(:,:,snaptime_now), [1,1,N]);
             xcopy = permute(repmat(data_now([2,3],:),[1,1,num_particles])...
                ,[1,3,2]);
             distance_pair = min(squeeze(sum((xcopy - ycopy).^2,1)/sigW^2), 300);
             p_ymkj2 = exp(-(distance_pair));
             scale_p_ymkj = 1./sum(p_ymkj2,1);
             scale_p_ymkj_mat = repmat(scale_p_ymkj,[num_particles,1]) ;
             p_ymkj2 = p_ymkj2 .* scale_p_ymkj_mat;
             tilde_p_ymk2 = mean(p_ymkj2,2);
             dEP_kr2 = 1/N *p_ymkj2 * deriv_loglike' ;
             derivative_now = derivative_now + sign(tilde_p_ymk2 - ...
                compress_snap_wgts2(:,snaptime_now))' * (dEP_kr2);             
             energy(snaptime_now) = energy(snaptime_now) + ...
                 sum(abs(tilde_p_ymk2 -compress_snap_wgts2(:,snaptime_now)));             
             %P3
             ycopy = repmat(snapshots3(:,:,snaptime_now), [1,1,N]);
             xcopy = permute(repmat(data_now([1,3],:),[1,1,num_particles])...
                ,[1,3,2]);
             distance_pair = min(squeeze(sum((xcopy - ycopy).^2,1)/sigW^2), 300);
             p_ymkj3 = exp(-(distance_pair));
             scale_p_ymkj = 1./sum(p_ymkj3,1);
             scale_p_ymkj_mat = repmat(scale_p_ymkj,[num_particles,1]) ;
             p_ymkj3 = p_ymkj3 .* scale_p_ymkj_mat;
             tilde_p_ymk3 = mean(p_ymkj3,2);
             dEP_kr3 = 1/N *p_ymkj3 * deriv_loglike' ;
             derivative_now = derivative_now + sign(tilde_p_ymk3 - ...
                compress_snap_wgts3(:,snaptime_now))' * (dEP_kr3);             
             energy(snaptime_now) = energy(snaptime_now) + ...
                 sum(abs(tilde_p_ymk3 -compress_snap_wgts3(:,snaptime_now)));  
             %Store
             derivative_alltime(:,snaptime_now) = derivative_now;
             slice_index = slice_index  + 1;
             %mean(data_now, 2)
        end %end of frame 
        
        time_now = time_now + deltat;

    end %end of time loop 
    derivative = mean(derivative_alltime(:,(include-1)),2);

end 

%takes in the current state x and returns the rate at the point 
function rate = compute_rate(x, theta)
    
    num_parameters = length(theta);
    [num_species , N ] = size(x);
    rate= zeros(num_parameters, N);
    rate(1,:) = theta(1);
    rate(2,:) = theta(2)*x(1,:);
    rate(3,:) = theta(3)*x(1,:);
    rate(4,:) = theta(4)*x(2,:);
end 


