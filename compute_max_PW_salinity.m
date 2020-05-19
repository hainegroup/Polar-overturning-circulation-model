% Function to compute the maximum polar water salinity using OC3D formula.
% twnh Aug '17, Jan, May '20

function S_PW = compute_max_PW_salinity(T_1,S_1,geophysical_parameters,static_parameters)

% Extract parameters:
pressure = geophysical_parameters.pressure ;
L        = static_parameters.L ;
c_p      = static_parameters.c_p ;
c_i      = static_parameters.c_i ;
T_i      = static_parameters.T_i ;
S_i      = static_parameters.S_i ;

% Compute:
T_f1  = gsw_CT_freezing(S_1    ,pressure) ;      % Use freezing temperature at temperature T_1.
alpha = gsw_alpha(      S_1,T_1,pressure) ;
beta  = gsw_beta(       S_1,T_1,pressure) ;
L1    = L  + c_i.*(T_f1 - T_i) ;
L2    = L1 + c_p.*(T_1 - T_f1) ;
S_PW  = (beta.*(S_1 - S_i).*S_1.*L1 + alpha.*(T_1 - T_f1).*S_i.*L2)./(beta.*(S_1 - S_i).*L1 + alpha.*(T_1 - T_f1).*L2) ;

end