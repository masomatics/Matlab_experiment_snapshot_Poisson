
%%
% This is the diffusion version of Dave's CFD paper. 
%
% dXt = A Xt dt + dW.  
% 
% rnsource  2 by  sample number by numtimepoints 
% dataset   2 by  sample number by numtimepoints 
% init      2 by 1
%% 

function [timemat, datmat, derivdat] = and_CFD_datagen_mass_deriv(init, tend, theta, sigV, num_timepts, rnsource, N)
    h = waitbar(0,'Adabra Catabra'); 
    delta = tend / num_timepts;
    

    rxn_mat =    [1     0       -1      0       0; 
                  0     1       0       -2      2;
                  0     0       0       1       -1]; 
    rxn_rate = zeros(5,N);
    
    
        
    timemat = delta*(0:num_timepts);
    
    datmat = zeros(3, N, num_timepts+1);
    derivdat = NaN(5,N);
    deriv_loglike = zeros(5,N);
    
    datmat(:,:,1) = repmat(init, 1,N);    
    matgrowth = repmat([1;0],1,N);
    for(k = 1 : num_timepts)
        
        rxn_rate(1,:) = theta(1); 
        rxn_rate(2,:) = theta(2) * datmat(1,:,k);
        rxn_rate(3,:) = theta(3) * datmat(1,:,k);
        rxn_rate(4,:) = theta(4) * datmat(2,:,k).^2;
        rxn_rate(5,:) = theta(5) * datmat(3,:,k);
        
        
        
        waitbar(k/num_timepts);
        datmathat = datmat(:, :, k) + rxn_mat*rxn_rate*delta;
        datmat(:, :, k+1)  =  max(datmathat + sigV* sqrt(delta)* rnsource(:,:,k),0)  ;
        
        %deriv_loglike(1,:) = deriv_loglike(1,:) +   sum( (datmat(:, :, k+1)- datmathat).*([-1, 0 ;0,0] * datmat(:, :, k)), 1) /(sigV^2) ;
        
        deriv_loglike(1,:) = deriv_loglike(1,:) +   (rxn_mat(:,1)' *datmat(:,:,k+1)).*datmat(1, :, k) /(sigV^2) ;
        deriv_loglike(2,:) = deriv_loglike(2,:) +   (rxn_mat(:,2)' *datmat(:,:,k+1)).*datmat(1, :, k)/(sigV^2) ;
        deriv_loglike(3,:) = deriv_loglike(3,:) +   (rxn_mat(:,3)' *datmat(:,:,k+1)).*datmat(1, :, k)/(sigV^2) ;
        deriv_loglike(4,:) = deriv_loglike(4,:) +   (rxn_mat(:,4)' *datmat(:,:,k+1)).*datmat(2, :, k).^2/(sigV^2) ;
        deriv_loglike(5,:) = deriv_loglike(5,:) +   (rxn_mat(:,5)' *datmat(:,:,k+1)).*datmat(3, :, k)/(sigV^2) ;
        
    end 
    derivdat =datmat(:,:,num_timepts+1).*repmat(deriv_loglike(3,:), [3,1]) ;
end 