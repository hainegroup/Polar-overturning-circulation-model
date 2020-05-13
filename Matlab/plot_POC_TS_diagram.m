% Make a plot of the TS diagram
% twnh May '20

function  objects = plot_POC_TS_diagram(this_case)

% TS diagram
theta      = -2.5:0.05:5 ;
salts      = 33  :0.02:37.0 ;
dens_labs  = 18  :0.25:32 ;
dens_labs2 = 18  :0.50:32 ;
[Y,X]      = meshgrid(theta,salts) ;
sigma0     = gsw_sigma0(X,Y) ;
iceT       = gsw_CT_freezing(salts,0) ;
% Cut out densities below the freezing point.
for ss = 1:length(salts)
    inds2            = find(Y(ss,:)<iceT(ss)) ;
    sigma0(ss,inds2) = NaN(length(inds2),1) ;
end % ss

% Plot TS diagram
[cs,h] = contour(X,Y,sigma0,dens_labs,'k') ;
hold on
plot(salts,iceT,'k-','linewidth',2)
clabel(cs,h,dens_labs2,'labelspacing',0.5*288,'fontsize',7)
axis([min(salts) max(salts) min(theta) max(theta)]) ;

% Plot results from theory
if(this_case.plot_options.diag_level > 0)
    p_theory  = plot(this_case.theory.S3s,this_case.theory.denom0.vals    ,'color',this_case.theory.denom0.color) ;
    p_theory2 = plot(this_case.theory.S3s,this_case.theory.numer0.U_2.vals,'color',this_case.theory.numer0.U_2.color) ;
    p_theoryi = plot(this_case.theory.S3s,this_case.theory.numer0.U_i.vals,'color',this_case.theory.numer0.U_i.color) ;
    if(this_case.plot_options.diag_level > 1)
        p_theory_S_star = plot(this_case.theory.numer0.U_2.S_star,this_case.theory.numer0.U_2.Tf2,'color',this_case.theory.numer0.U_2.color,'Marker','diamond') ;
    end
end % if
if(isfield(this_case.theory.patch_data,'S1'))
    p1 = patch(this_case.theory.patch_data.S1,this_case.theory.patch_data.T1,this_case.plot_options.OW_col) ;
    set(p1,'edgecolor','none','FaceAlpha',this_case.plot_options.transp_factor) ;
    uistack(p1,'bottom') ;
end % if
if(isfield(this_case.theory.patch_data,'S2'))
    p2 = patch(this_case.theory.patch_data.S2,this_case.theory.patch_data.T2,this_case.plot_options.OW_col) ;
    set(p2,'edgecolor','none','FaceAlpha',this_case.plot_options.transp_factor) ;
    uistack(p2,'bottom') ;
end % if
if(isfield(this_case.theory.patch_data,'S3'))
    p3 = patch(this_case.theory.patch_data.S3,this_case.theory.patch_data.T3,this_case.plot_options.aW_col) ;
    set(p3,'edgecolor','none','FaceAlpha',this_case.plot_options.transp_factor) ;
    uistack(p3,'bottom') ;
end % if

% Plot water masses
p_AW  = scatter(this_case.base_case.S_1,this_case.base_case.T_1 ,80,'o','filled','Markerfacecolor',this_case.plot_options.AW_col) ;
p_PW  = scatter(this_case.dynamic_parameters.S_2,this_case.dynamic_parameters.Tf2,80,'o','filled','Markerfacecolor',this_case.plot_options.PW_col) ;
if(isfield(this_case,'vals'))
    p_OW_vals = scatter(cell2mat(this_case.vals.S_3),cell2mat(this_case.vals.T_3),20,'o','filled','Markerfacecolor',this_case.plot_options.OW_col) ;
    p_SW_vals = scatter(cell2mat(this_case.vals.S_s),cell2mat(this_case.vals.T_f),20,'o','filled','Markerfacecolor',this_case.plot_options.SW_col) ;
end % if
p_SW  = scatter(this_case.sens.S_scan,this_case.sens.Tfs,80,'o','filled','Markerfacecolor',this_case.plot_options.SW_col,'MarkerEdgecolor','w') ;
p_aW  = scatter(this_case.dynamic_parameters.S_a,this_case.dynamic_parameters.T_a,80,'x','Markerfacecolor',this_case.plot_options.aW_col,'MarkerEdgecolor',this_case.plot_options.aW_col) ;
p_OW  = scatter(this_case.sens.S_3can,this_case.sens.T_3can ,80,'o','filled','Markerfacecolor',this_case.plot_options.OW_col,'MarkerEdgecolor','w') ;

xlabel('Salinity [g/kg]')
ylabel('Temperature [^{\rm o}C]')
grid on
if(this_case.plot_options.diag_level > 1)
    legend([p_AW,p_PW,p_OW,p_SW,p_aW,p_theory,p_theory2,p_theoryi],...
        'AW','PW','OW','SW','aW',this_case.theory.denom0.legend,this_case.theory.numer0.U_2.legend,this_case.theory.numer0.U_i.legend,'location','eastoutside','interpreter','latex')
else
    legend([p_AW,p_PW,p_OW,p_SW,p_aW],...
        'AW','PW','OW','SW','aW','location','eastoutside','interpreter','latex')
end % if
set(gca,'Box','on') ;
title([this_case.base_case.exp_str,' ',this_case.base_case.tit_str],'Interpreter','latex','FontSize',14)

pos = get(gca,'Position') ;
set(gca,'Position',[0.3 pos(2) 0.35 pos(4)]) ;

% Output structure for interactive use
objects = struct('p_AW',p_AW,'p_OW',p_OW,'p_PW',p_PW,'p_aW',p_aW,'p_SW',p_SW,'p_OW_vals',p_OW_vals,'p_SW_vals',p_SW_vals,'p_theory_S_star',p_theory_S_star) ;
if(this_case.plot_options.diag_level > 0)
    objects.p_theory  = p_theory ;
    objects.p_theory2 = p_theory2 ;
    objects.p_theoryi = p_theoryi ;
    objects.p1        = p1 ;
    objects.p2        = p2 ;
    objects.p3        = p3 ;
end % if

end