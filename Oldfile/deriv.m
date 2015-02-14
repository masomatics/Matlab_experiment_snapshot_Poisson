function derivdata = deriv(datmat, timepts, theta, sigV)

    [N, num_timepts] = size(datmat);
    derivdata = NaN(1,N); 
    for k = 1:N
        derivdata(k) = 0;
        deriv_loglike = 0; 
        for j = 2:num_timepts
            
            delta = timepts(j) -  timepts(j-1);
            xdat_estimate = datmat(k, j-1) + theta*datmat(k, j-1)*delta; 
%            loglike = loglike +  datmat(k, j-1)*delta ...
%                        *(datmat(k, j) -  xdat_estimate)/(sigV^2*delta); 

            deriv_loglike = deriv_loglike +  datmat(k, j-1)* ...
                        (datmat(k, j) -  xdat_estimate)/(sigV^2); 


            
        end  
        derivdata(k) = datmat(k,num_timepts)* deriv_loglike;
    end
        
end 