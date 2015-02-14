
%%
%
% rxn1 0 -> M 
% rxn2 M -> M + P 
% rxn3 M -> 0
% rxn4 M+ P -> M
% rxn5 P -> D
% Default =  
%
%% 

function [timemat, datmat] = and_CFD_datagen_mass_diffusion(init, tend, theta, sigV, num_timepts, rnsource, N)
    h = waitbar(0,'Adabra Catabra'); 
    delta = tend / num_timepts;
    rxn_mat =    [1     0       -1      0       0; 
                  0     1       0       -1      0;
                  0     0       0       0       1]; 
    rxn_rate = zeros(5,N);
    
    
        
    timemat = delta*(0:num_timepts);
    
    datmat = zeros(3, N, num_timepts+1);
    datmathat = zeros(3,N);
    
    datmat(:,:,1) = repmat(init, 1,N);    
    matgrowth = repmat([1;0],1,N);
    for(k = 1 : num_timepts)
        
        rxn_rate(1,:) = theta(1); 
        rxn_rate(2,:) = theta(2) * datmat(1,:,k);
        rxn_rate(3,:) = theta(3) * datmat(1,:,k);
        rxn_rate(4,:) = theta(4) * datmat(2,:,k).*datmat(1,:,k);
        rxn_rate(5,:) = theta(5) * datmat(2,:,k);
        
        
        
        waitbar(k/num_timepts);
        datmathat = datmat(:, :, k) + rxn_mat*rxn_rate*delta;
        datmat(:, :, k+1)  =  max(datmathat + sigV* sqrt(delta)* rnsource(:,:,k),0)  ;
        %datmat(:, :, k+1)= datmathat + sigV* sqrt(delta)* rnsource(:,:,k);
    end 
    close(h)
end 