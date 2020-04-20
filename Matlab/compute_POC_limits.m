% Function to compute the Q+, Q-, Qstar limits for U_2 and U_3 vanishing.
% twnh Mar '20
function param = compute_POC_limits(this_case,inds,flag)

% Define intermediate variables. They are used to assemble the output by
% computing a non-dimensional compound forcing variable.
tmp1  = 1 - this_case.static_parameters.S_i./this_case.dynamic_parameters.S_1(inds) ;
tmp2  = this_case.sens.Lpr(inds).*this_case.sens.F(inds) ;
tmp3  = this_case.dynamic_parameters.rh1(inds).*this_case.static_parameters.c_p.*this_case.dynamic_parameters.T_1(inds).*this_case.dynamic_parameters.U_1(inds) ;
tmp4  = this_case.static_parameters.S_i./this_case.dynamic_parameters.S_1(inds) ;

switch flag
    case 'salt'
        % Salt crisis: zero entrainment
        [this_S,this_T] = compute_SW_properties_zero_Phi(this_case,inds) ;
        Q_theory        = compute_limit(this_case,inds,this_T,this_S) ;
        param           = (tmp1.*Q_theory + tmp2)./tmp3 + tmp4 ;
        
    case 'emergency'
        % OW emergency: U_3 = 0
        Q_theory = compute_limit(this_case,inds,this_case.sens.Tf2(inds),this_case.sens.S_2(inds)) ;
        param    = (tmp1.*Q_theory + tmp2)./tmp3 + tmp4 ;
        
    case 'heat'
        % Heat crisis: This is tricky and there are a few theoretical options:
        
        % 0. Entrainment equals one. This underestimates the forcing variable because for phi > 0 these
        % T and S values give OW denser than AW. Included for completeness.
        Q_theory0  = compute_limit(this_case,inds,this_case.sens.T_a(inds),this_case.sens.S_a(inds)) ;
        param.mid0 = (tmp1.*Q_theory0 + tmp2)./tmp3 + tmp4 ;
        
        % 1. Entrainment is maximized (which is less than one), but this is tricky to get
        % accurate. It also depends on the maximum shelf salinity (which is
        % fragile). And even then, it's an approximation.
        [this_S,this_T, this_Phi, this_Ss] = compute_OW_properties_max_Phi(this_case,inds) ;
        Q_theory1   = compute_limit(this_case,inds,this_T,this_S) ;
        param.mid1  = (tmp1.*Q_theory1 + tmp2)./tmp3 + tmp4 ;
        param.T_3_1 = this_T ;
        param.S_3_1 = this_S ;
        param.Phi_1 = this_Phi ;
        param.S_s_1 = this_Ss ;
        
        % 2. Use limits of AW properties to define limits of possible AW. This is a
        % robust limit (static stability of OW and AW), but gives wide bounds.
        % 2a. Limiting case for phi = 1.
        [Sf,Tf]       = compute_freezing_properties_AW_density(this_case,inds) ;
        Q_theory2.hi  = compute_limit(this_case,inds,Tf,Sf) ;
        param.hi      = (tmp1.*Q_theory2.hi  + tmp2)./tmp3 + tmp4 ;
        % 2b. Limiting case for this phi and S_s = infinity.
        this_T        = this_case.sens.T_a(inds) ;
        this_S        = gsw_SA_from_rho(this_case.dynamic_parameters.rh1(inds),this_T,this_case.geophysical_parameters.pressure) ;
        Q_theory2.mid = compute_limit(this_case,inds,this_T,this_S) ;
        param.mid2    = (tmp1.*Q_theory2.mid + tmp2)./tmp3 + tmp4 ;
        param.T_3_2 = this_T ;
        param.S_3_2 = this_S ;
        % 2c. Limiting case for phi = 0.
        Q_theory2.lo  = compute_limit(this_case,inds,this_case.dynamic_parameters.T_1(inds),this_case.dynamic_parameters.S_1(inds)) ;
        param.lo      = (tmp1.*Q_theory2.lo  + tmp2)./tmp3 + tmp4 ;
        
        % 3. Assume we know T3 and S3. This is an unreasonable
        % assumption, but it tests if the theory and code works!
        Q_theory3   = compute_limit(this_case,inds,this_case.sens.T_3can(inds),this_case.sens.S_3can(inds)) ;
        param.cheat = (tmp1.*Q_theory3 + tmp2)./tmp3 + tmp4 ;
        
        % 4. Compute T3 and S3 by assuming U3 = -U1 (not exactly true), we know the ambient
        % water mixing fraction, and we know rho(T3,S3) = rho(T1,S1) (not always true).
        [this_S, this_T, this_Phi, this_Ss] = compute_OW_properties(this_case,inds) ;
        Q_theory4   = compute_limit(this_case,inds,this_T,this_S) ;
        param.mid4  = (tmp1.*Q_theory4 + tmp2)./tmp3 + tmp4 ;
        param.T_3_4 = this_T ;
        param.S_3_4 = this_S ;
        param.Phi_4 = this_Phi ;
        param.S_s_4 = this_Ss ;
        
end % switch

end

function Qlim = compute_limit(this_case,inds,T,S)
% Implements formula for the Q value that gives U2 or U3 = 0.

F    = this_case.sens.F(inds) ;
U1   = this_case.dynamic_parameters.U_1(inds) ;
T1   = this_case.dynamic_parameters.T_1(inds) ;
S1   = this_case.dynamic_parameters.S_1(inds) ;
rho1 = this_case.dynamic_parameters.rh1(inds) ;
Sice = this_case.static_parameters.S_i ;
c_p  = this_case.static_parameters.c_p ;
Lp   = this_case.sens.Lpr(inds) ;

fac1 = Lp.*S + c_p.*T.*Sice ;
fac2 = Lp    + c_p.*T ;
fac3 = S - Sice ;

Qlim = (fac1.*(rho1.*U1 - F) - rho1.*fac2.*S1.*U1)./fac3 + c_p*rho1.*U1.*T1 ;

end

function [Ss,Tf] = compute_SW_properties_zero_Phi(this_case,inds)
% Finds shelf water salinity and (freezing) temperature that gives zero
% entrainment.

U1       = this_case.dynamic_parameters.U_1(inds) ;
Ta       = this_case.sens.T_a(inds) ;
Sa       = this_case.sens.S_a(inds) ;
gamma    = this_case.static_parameters.gam ;
pressure = this_case.geophysical_parameters.pressure ;

rha = gsw_rho(Sa,Ta,pressure) ;    % Density of ambient  Water, kg/m^3
rhs = gamma^(3/2).*abs(U1).^(1/2) + rha ;

Ss = NaN(numel(Sa),1) ;
Tf = NaN(numel(Sa),1) ;
init_Tf = gsw_CT_freezing(nanmean(Sa),pressure) ;     % Initial guess of freezing temperature.

% Find shelf salinity and freezing temperature by searching for each zero
% crossing.
for ii = 1:length(rhs)
    if(~isnan(rhs(ii)))
        % This function crosses zero for the freezing temperature that
        % gives density rhs(ii).
        minimizer = @(Tf) gsw_rho(gsw_SA_freezing_from_CT(Tf,pressure),Tf,pressure) - rhs(ii) ;
        Tf(ii) = fzero(minimizer,init_Tf) ;
        Ss(ii) = gsw_SA_freezing_from_CT(Tf(ii),pressure) ;
    end % if
end % ii

end

function [Sf,Tf] = compute_freezing_properties_AW_density(this_case,inds)
% Finds OW salinity and temperature assuming that OW density equals the AW
% density (to maintain marginal static stability) and the OW is at the 
% freezing temperature.

S1       = this_case.dynamic_parameters.S_1(inds) ;
rh1      = this_case.dynamic_parameters.rh1(inds) ;
N        = numel(inds) ;
pressure = this_case.geophysical_parameters.pressure ;

% This function crosses zero for the freezing temperature that
% gives AW density.
init_Tfs = gsw_CT_freezing(S1,pressure) ;
Tf       = NaN(N,1) ;
Sf       = NaN(N,1) ;
for ii = 1:N
    minimizer = @(Tf) gsw_rho(gsw_SA_freezing_from_CT(Tf,pressure),Tf,pressure) - rh1(ii) ;
    Tf(ii)    = fzero(minimizer,init_Tfs(ii)) ;
    Sf(ii)    = gsw_SA_freezing_from_CT(Tf(ii),pressure) ;
end % ii

end

function [S3,T3,Phi_max,Ss] = compute_OW_properties_max_Phi(this_case,inds)
% Finds OW salinity and temperature that gives maximum entrainment for
% upper limit of Smax on SW salinity. The formula assumes that U_3 = - U_1,
% which is not exactly right because it neglects U_i>0. The formula also
% assumes knowledge of Smax and phi, which are themselves, fragile.

U1       = this_case.dynamic_parameters.U_1(inds) ;
Ta       = this_case.sens.T_a(inds) ;
Sa       = this_case.sens.S_a(inds) ;
gamma    = this_case.static_parameters.gam ;
pressure = this_case.geophysical_parameters.pressure ;
Smax     = this_case.static_parameters.S_ssmax ;
DPhi     = this_case.static_parameters.Dphi_limit ;

rha      = gsw_rho(Sa,Ta,pressure) ;              % Density of ambient  Water, kg/m^3
Tf_max   = gsw_CT_freezing(Smax,pressure) ;       % Maximum possible SW density
rhs_max  = gsw_rho(Smax,Tf_max,pressure) ;        % Maximum possible SW density
Phi_max  = 1 - (gamma^(3/2)).*(abs(U1).^(1/2))./(rhs_max - rha) + DPhi ; % Add DPhi because that's the largest acceptable Phi.

T3       = Tf_max - (Tf_max - Ta).*Phi_max ;      % Temperature of Overflow Water, degC
S3       = Smax   - (Smax   - Sa).*Phi_max ;      % Salinity    of Overflow Water, g/kg
Ss       = Smax.*ones(size(S3)) ;

end

function [S3,T3,best_Phi,Ss] = compute_OW_properties(this_case,inds)
% Finds OW salinity and temperature that solves for entrainment and Ss
% assuming U_2 = 0, U_3 = -U_1 (not exactly right because it neglects U_i),
% and we know the ambient water mixing fraction, phi (hence Ta, Sa).

% Numerically THIS GIVES SOLUTIONS WITH Ss ~ Inf and Phi ~ 1. These aren't
% very REALISTIC!!!

T1       =  this_case.dynamic_parameters.T_1(inds) ;
S1       =  this_case.dynamic_parameters.S_1(inds) ;
U1       =  this_case.dynamic_parameters.U_1(inds) ;
Ta       =  this_case.sens.T_a(inds) ;
Sa       =  this_case.sens.S_a(inds) ;
gam_fac  = (this_case.static_parameters.gam)^(3/2) ;
pressure =  this_case.geophysical_parameters.pressure ;

rha      = gsw_rho(Sa,Ta,pressure) ;               % Density of ambient  Water, kg/m^3
rh1      = gsw_rho(S1,T1,pressure) ;               % Density of Atlantic Water, kg/m^3
Tf       = gsw_CT_freezing(Sa,pressure) ;

% %% Testing!
% Smax     = this_case.static_parameters.S_ssmax ;
% DPhi     = this_case.static_parameters.Dphi_limit ;
% 
% Tf_max   = gsw_CT_freezing(Smax,pressure) ;        % Maximum possible SW density
% rhs_max  = gsw_rho(Smax,Tf_max,pressure) ;        % Maximum possible SW density
% Phi_max  = 1 - gam_fac.*(abs(U1).^(1/2))./(rhs_max - rha) + DPhi ; % Add DPhi because that's the largest acceptable Phi.
% %%% WARNING!! TESTING!!
% U1       =  abs(this_case.sens.U_3can(inds)) ;
% % End testing!

% Loop over the cases of interest in the inds variable and compute T3 and S3.
N        = numel(inds) ;
T3       = NaN(N,1) ;
S3       = NaN(N,1) ;
best_Phi = NaN(N,1) ;
Ss       = NaN(N,1) ;
for ii = 1:N
    T3_fn        = @(Phi) Tf(ii) - (Tf(ii) - Ta(ii))*Phi ;
    S3_fn        = @(Phi) gsw_SA_from_rho(rh1(ii).*ones(size(Phi)),T3_fn(Phi),pressure);
    Ss_fn        = @(Phi) (S3_fn(Phi) - Sa(ii).*Phi)./(1 - Phi) ;
    minimizer    = @(Phi) (Phi - 1 + (sqrt(U1(ii))*gam_fac)./(gsw_rho(Ss_fn(Phi),Tf(ii).*ones(size(Phi)),pressure) - rha(ii))).^2 ;
    best_Phi(ii) = fminbnd(minimizer,0,1) ;
    %best_Phi(ii) = fminbnd(minimizer,0,Phi_max(ii)) ;

    T3(ii)       = T3_fn(best_Phi(ii)) ;
    S3(ii)       = S3_fn(best_Phi(ii)) ;
    Ss(ii)       = Ss_fn(best_Phi(ii)) ;
end % ii

end