% Function to progressively compute solution to the polar overturning
% conceptual model.
% The progmat matrix records the progress reached through the progressive
% solve:
% progmat(ss,pp,ff) Stage passed
%   1               U_2 negative
%   2               U_3 negative
%   3               U_i negative
%   4               Good Phi fit
%   5               F_2/Q_2 inequalities
%   6               u_1/u_i inequalities
%   7               OW/AW density stability
% twnh Sep, Nov '19, Jan, Mar '20

function [Phis,S_ss,progmat,...
    U_is,U_2s,U_3s,u_ss,T_3s,S_3s,rho_3s,rho_ss,T_fss,...
    Q_2_los,Q_2_his,F_2_los,F_2_his,u_1_los,u_1_his,u_i_los,u_i_his,...
    S_s0,DPhis,Lfacs,...
    S_2,T_a,S_a,rho_2,T_f2,mat1_cd] = solve_POC_model(dynamic_parameters,this_index,static_parameters,geophysical_parameters)

%% Extract variables from arguments
Q          = dynamic_parameters.Q(this_index)  
F          = dynamic_parameters.F(this_index) ;
U_1        = dynamic_parameters.U_1(this_index) ;
T_1        = dynamic_parameters.T_1(this_index) ;
S_1        = dynamic_parameters.S_1(this_index) ;
rh1        = dynamic_parameters.rh1(this_index) ;
DS2        = dynamic_parameters.DS2(this_index) ;
phi        = dynamic_parameters.phi(this_index) ;
rho_i      = static_parameters.rhi ;
T_i        = static_parameters.T_i ;
S_i        = static_parameters.S_i ;
c_p        = static_parameters.c_p ;
c_i        = static_parameters.c_i ;
gamma      = static_parameters.gam ;
L          = static_parameters.L ;
DPhi_limit = static_parameters.DPhi_limit ;
pressure   = geophysical_parameters.pressure ;

% Compute PW and ambient water properties, plus AW density. These variables
% are also passed out.
S_2   = compute_max_PW_salinity(T_1, S_1, pressure, L, c_p, c_i, T_i, S_i) + DS2 ;
T_f2  = gsw_CT_freezing(S_2,pressure) ;  % Freezing temperature of Polar    Water, degC
rho_2 = gsw_rho(S_2, T_f2, pressure) ;   % Density              of Polar    Water, kg/m^3
T_a   = phi*T_f2 + (1 - phi)*T_1 ;       % Temperature          of ambient  Water, degC
S_a   = phi*S_2  + (1 - phi)*S_1 ;       % Salinity             of ambient  Water, g/kg
rho_a = gsw_rho(S_a, T_a, pressure) ;    % Density              of ambient  Water, kg/m^3

%% Compute progressively
Phis_vec = linspace(0,1,static_parameters.N_Phis_vec)' ;
N_Phis   = static_parameters.N_Phis_vec ;
S_s0     = gsw_SA_from_rho(rh1,T_f2,pressure) ;        % Not exactly the right freezing temperature, but ensures that S_s0 is a slight overestimate.
S_ss_vec = linspace(S_s0,static_parameters.S_ssmax,static_parameters.N_S_ss_vec)' ;
N_S_ss   = static_parameters.N_S_ss_vec ;

% Specify static matrices and vectors. NaNs get replaced in the loop.
rhs1    = [
    F -     rh1*U_1     ; ...
    -     rh1*U_1*S_1 ; ...
    Q - c_p*rh1*U_1*T_1 ] ;

matrix1 = [
    rho_2      NaN  rho_i     ; ...
    rho_2*S_2  NaN  rho_i*S_i ; ...
    c_p*rho_2*T_f2 NaN  NaN       ] ;

matrix3 = [
    rh1      rho_i     ; ...
    rh1*S_1  rho_i*S_i ; ...
    c_p*rh1*T_1  NaN       ] ;

% Allocate output matrices
progmat  = NaN(N_S_ss,N_Phis,7) ;
Q_2_los  = NaN(N_S_ss,N_Phis) ;
Q_2_his  = Q_2_los ;
F_2_los  = Q_2_los ;
F_2_his  = Q_2_los ;
u_1_los  = Q_2_los ;
u_1_his  = Q_2_los ;
u_i_los  = Q_2_los ;
u_i_his  = Q_2_los ;
DPhis    = Q_2_los ;
U_is     = Q_2_los ;
U_2s     = Q_2_los ;
U_3s     = Q_2_los ;
u_ss     = Q_2_los ;
% Condition number diagnostics for debug:
mat1_cd  = Q_2_los ;

% Precompute as much as possible
[S_ss, Phis] = ndgrid(S_ss_vec,Phis_vec) ;
T_fss        = gsw_CT_freezing(S_ss,pressure) ;         % Freezing temperature of Shelf Water, oC.
rho_ss       = gsw_rho(S_ss, T_fss,pressure) ;          % Density              of Shelf Water, kg/m^3
T_3s         = T_fss - (T_fss - T_a).*Phis ;            % Temperature          of Overflow Water, degC
S_3s         = S_ss  - (S_ss  - S_a).*Phis ;            % Salinity             of Overflow Water, g/kg
rho_3s       = gsw_rho(S_3s, T_3s, pressure) ;          % Density              of Overflow Water, kg/m^3
Lfacs        = L - c_p.*T_fss + c_i.*(T_fss - T_i) ;

% Loop over S_s and Phi
for ss = 1:N_S_ss
    Lfac_s       = Lfacs(ss,1) ;
    matrix1(3,3) = -rho_i*Lfac_s ;
    matrix3(3,2) = -rho_i*Lfac_s ;
    soln_mat     = (matrix3'*matrix3)\matrix3' ;
    
    for pp = 1:N_Phis
        Phi   = Phis(  ss,pp) ;
        T_3   = T_3s(  ss,pp) ;
        S_3   = S_3s(  ss,pp) ;
        rho_3 = rho_3s(ss,pp) ;
        S_s   = S_ss(  ss,pp) ;
        T_fs  = T_fss( ss,pp) ;
        rho_s = rho_ss(ss,pp) ;
        
        % Test for OW static stability
        if(rho_3 < rh1)
            progmat(ss,pp,7) = 0 ;
        else
            progmat(ss,pp,7) = 1 ;
        end % if
        
        if( progmat(ss,pp,7) == 1)
            % Define and solve the problem for [U_2, U_3, U_i], the basin+shelf system:
            matrix1(:,2) = [rho_3; rho_3*S_3; c_p*rho_3*T_3] ;
            soln1 = matrix1\rhs1 ;
            U_2   = soln1(1) ;
            U_3   = soln1(2) ;
            U_i   = soln1(3) ;
            
            % Check satisfaction of inequality constraints at stage 1.
            if(U_2 <= 0)
                progmat(ss,pp,1) = 1 ;
            else
                progmat(ss,pp,1) = 0 ;
            end % if
            if(U_3 <= 0)
                progmat(ss,pp,2) = 1 ;
            else
                progmat(ss,pp,2) = 0 ;
            end % if
            if(U_i <= 0)
                progmat(ss,pp,3) = 1 ;
            else
                progmat(ss,pp,3) = 0 ;
            end % if
            
            if(progmat(ss,pp,1) == 1 && progmat(ss,pp,2) == 1 && progmat(ss,pp,3) == 1)
                % Compute u_s
                u_s = U_3*(1 - Phi) ;
                
                % Check entrainment Phi. Compute difference and compare to threshold:
                DPhi = Phi - (1 - (gamma*(abs(u_s)^(1/3))*(rho_s - rho_a)^(-2/3))) ;
                if(abs(DPhi) <= DPhi_limit && Phi < 1)          % Avoid Phi = 1.
                    progmat(ss,pp,4) = 1 ;
                    % Solve overdetermined problem for (u_1, u_i) using this u_s
                    % and the inequality constraints.
                    rhs3       = [
                        -     rho_s*u_s      ; ...
                        -     rho_s*u_s*S_s  ; ...
                        - c_p*rho_s*u_s*T_fs ] ;
                    [U,~,~]    = svd(matrix3) ;
                    range_vec3 = U(:,3) ;
                    
                    % Form constraint equation for F_2 in terms of Q_2
                    % F_2 = const1 + slope1*Q_2
                    const1 = -(range_vec3'*rhs3)/range_vec3(1) ;
                    slope1 = - range_vec3(3)    /range_vec3(1) ;
                    F_2_lo = const1 ;                          % Value of F_2 for minimum Q_2, which is zero.
                    F_2_hi = const1 + slope1*Q ;               % Value of F_2 for maximum Q_2, which is Q.
                    
                    % Test these F_2 limits for consistency with the F_2 inequalities
                    if(~(F_2_lo > 0 && F_2_hi > 0) || (F_2_lo < F && F_2_hi < F))    % Consistent
                        progmat(ss,pp,5) = 1 ;
                        % Now find the consistent (F_2,Q_2) value with smaller Q_2 by testing F_2_lo.
                        F_2_lo = min([F_2_lo         0]) ;
                        F_2_lo = max([F_2_lo F]) ;
                        Q_2_lo = (F_2_lo - const1)/slope1 ;
                        % Now find the consistent (F_2,Q_2) value with larger Q_2 by testing F_2_hi.
                        F_2_hi = min([F_2_hi         0]) ;
                        F_2_hi = max([F_2_hi F]) ;
                        Q_2_hi = (F_2_hi - const1)/slope1 ;
                        
                        % With the line segment of (Q_2_lo,F_2_lo) and
                        % (Q_2_hi,F_2_hi) values, we can now compute (u_1, u_i)
                        % solutions and test them for satisfaction of
                        % the inequalities on (u_1, u_i).
                        y_vec_lo = [F_2_lo; 0; Q_2_lo] + rhs3 ;
                        y_vec_hi = [F_2_hi; 0; Q_2_hi] + rhs3 ;
                        soln2    = soln_mat*[y_vec_lo y_vec_hi] ;
                        u_1_lo   = soln2(1,1) ;
                        u_i_lo   = soln2(2,1) ;
                        u_1_hi   = soln2(1,2) ;
                        u_i_hi   = soln2(2,2) ;
                        
                        % Test the (u_1, u_i) inequalities.
                        if((u_1_lo > 0 && u_i_lo < U_i) || (u_1_hi > 0 && u_i_hi < U_i))    % At least one point is consistent
                            progmat(ss,pp,6) = 1 ;
                            
                            % Adjust one of the (u_1, u_i) points if necessary.
                            % Get equation for the (u_1, u_i) line:
                            % u_i = const2 + slope2*u_1
                            slope2 = (u_i_hi - u_i_lo)/(u_1_hi - u_1_lo) ;
                            const2 = u_i_lo - slope2*u_1_lo ;
                            
                            % Test each coordinate and adjust if necessary
                            if(u_1_lo < 0)
                                u_1_lo = 0 ;
                                u_i_lo = const2 + slope2*u_1_lo ;
                            end % if
                            if(u_1_hi < 0)
                                u_1_hi = 0 ;
                                u_i_hi = const2 + slope2*u_1_hi ;
                            end % if
                            if(u_i_lo > U_i)
                                u_i_lo = U_i ;
                                u_1_lo = (u_i_lo - const2)/slope2 ;
                            end % if
                            if(u_i_hi > U_i)
                                u_i_hi = U_i ;
                                u_1_hi = (u_i_hi - const2)/slope2 ;
                            end % if
                            
                            soln3 = [
                                u_1_lo, u_1_hi ;
                                u_i_lo, u_i_hi ] ;
                            
                            % Adjust one of the (F_2, Q_2) points if necessary.
                            soln4 = matrix3*soln3 - rhs3 ;
                            F_2_lo = soln4(1,1) ;
                            Q_2_lo = soln4(3,1) ;
                            F_2_hi = soln4(1,2) ;
                            Q_2_hi = soln4(3,2) ;
                            
                            % Build output data
                            F_2_los(ss,pp) = F_2_lo ;
                            F_2_his(ss,pp) = F_2_hi ;
                            Q_2_los(ss,pp) = Q_2_lo ;
                            Q_2_his(ss,pp) = Q_2_hi ;
                            u_1_los(ss,pp) = u_1_lo ;
                            u_1_his(ss,pp) = u_1_hi ;
                            u_i_los(ss,pp) = u_i_lo ;
                            u_i_his(ss,pp) = u_i_hi ;
                            U_is(   ss,pp) = U_i ;
                            U_2s(   ss,pp) = U_2 ;
                            U_3s(   ss,pp) = U_3 ;
                            u_ss(   ss,pp) = u_s ;
                            DPhis(  ss,pp) = DPhi ;
                            mat1_cd(ss,pp) = cond(matrix1) ;
                                                        
                        end % if
                    end % if
                end % if
            end % if
        end % if
    end % ss
end % pp

end

% Function to compute the maximum polar water salinity using OC3D formula.
% twnh Aug '17, Jan '20

function S_PW = compute_max_PW_salinity(T_1, S_1, pressure, L, c_p, c_i, T_i, S_i)

T_f1  = gsw_CT_freezing(S_1    ,pressure) ;      % Use freezing temperature at temperature T_1.
alpha = gsw_alpha(      S_1,T_1,pressure) ;
beta  = gsw_beta(       S_1,T_1,pressure) ;
L1    = L  + c_i*(T_f1 - T_i) ;
L2    = L1 + c_p*(T_1 - T_f1) ;
S_PW  = (beta*(S_1 - S_i)*S_1*L1 + alpha*(T_1 - T_f1)*S_i*L2)/(beta*(S_1 - S_i)*L1 + alpha*(T_1 - T_f1)*L2) ;

end