
%
% Solution to ther Anderson CFD. Use this to confirm the formulas for the
% first and the second derivative. 
%
%

theta(1) = 2 ; 
theta(2) = 10;
theta(3) = 1/4;
theta(4) = 1;

%%

t = 3; 
X0 = [0; 0]

A = [-theta(3), 0 ; theta(2) , -theta(4)]; 
b = [theta(1); 0];

P   = [theta(4) - theta(3) , 0; theta(2), 1];
Pin = [1, 0 ; -theta(2), theta(4) - theta(3)]./(theta(4) - theta(3)) ;




D = diag([-theta(3) , - theta(4)]);
expDt = diag([ exp(-theta(3)*t) , exp(- theta(4)*t)])
intD_t = diag([(exp(theta(3)*t)-1)/theta(3), (exp(theta(4)*t)-1)/theta(4)]) 

expAt = P * expDt * Pin;

int_expA_t_s =   P*(expDt * intD_t)*Pin

Xt = expAt * X0 + int_expA_t_s * b

%Full Formula
Xt = P* diag([exp(-theta(3)*t), exp(-theta(4)*t)]) * Pin  * X0+ ... 
P*diag([1/theta(3)*(1- exp(-theta(3)*t)) ,1/theta(4)*(1- exp(-theta(4)*t))])*Pin...
*[theta(1);0]


DD = diag([1/theta(3)*(1- exp(-theta(3)*t)) ,1/theta(4)*(1- exp(-theta(4)*t))]);
rho_3P = [-1, 0; 0,0];
rho_3Pin = [0, 0 ; 0, -1] ./(theta(4) - theta(3)) + [1, 0 ; -theta(2), theta(4) - theta(3)]./((theta(4) - theta(3))^2)
rho_3DD = diag([-1/(theta(3)^2)*(1- exp(-theta(3)*t)) + 1/theta(3) * t* exp(-theta(3)*t),0]);


rho_3Xt = (P* diag([-t* exp(-theta(3)*t), 0 ]) * Pin + ... 
rho_3P*     diag([exp(-theta(3)*t), 0 ]) *Pin + ...
P*     diag([exp(-theta(3)*t), 0 ]) *rho_3Pin) * X0+ ... 
(P*diag([-1/(theta(3)^2)*(1- exp(-theta(3)*t)) + 1/theta(3) * t* exp(-theta(3)*t),0])*Pin +...
rho_3P*diag([1/theta(3)*(1- exp(-theta(3)*t)) ,1/theta(4)*(1- exp(-theta(4)*t))])*Pin + ...
P*diag([1/theta(3)*(1- exp(-theta(3)*t)) ,1/theta(4)*(1- exp(-theta(4)*t))])*rho_3Pin)...
*[theta(1);0]

rho_3Xt = (P*rho_3DD*Pin + rho_3P*DD*Pin + P * DD * rho_3Pin)* [theta(1);0]


%I will ignore the first term since X0 in our example is the zero vector. 
rho_33P = [0, 0; 0,0];
rho_33Pin = [0, 0 ; 0, -1] ./(theta(4) - theta(3))^2 + ...
    [0, 0 ; 0, -1]./((theta(4) - theta(3))^2) + ...
    [1, 0 ; -theta(2), theta(4) - theta(3)]./((theta(4) - theta(3))^3)*2 ;   

rho_33DD = diag([ ...
2/(theta(3)^3)*(1- exp(-theta(3)*t))-...    
2/(theta(3)^2)*t*exp(- theta(3)*t)+...
1/(theta(3))* (-t^2)*exp(- theta(3)*t) ...
,0]);

rho_33Xt= ( rho_3P*rho_3DD*Pin + P*rho_33DD*Pin + P*rho_3DD*rho_3Pin + ...
  rho_33P*DD*Pin + rho_3P*rho_3DD*Pin  + rho_3P*DD*rho_3Pin + ...
  rho_3P*DD*rho_3Pin + P*rho_3DD*rho_3Pin + P*DD*rho_33Pin)* [theta(1);0] 




