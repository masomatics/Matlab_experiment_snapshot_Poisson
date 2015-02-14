% Analytic Solution to the system of Anderson CFD.
%
%
%
function Xt = debug_sanity_check(X0, theta, t)

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
end 