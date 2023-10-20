% Make a plot of the TS diagram
% twnh May, Sep '20, Oct '23

function  objects = plot_POC_TS_diagram(this_case)

% TS diagram
persistent theta salts sigma0 iceT X Y dens_labs dens_labs2
if(isempty(theta) || this_case.base_case.app.replot_TS_diagram)
    theta      = -2.5:0.05:12 ;
    salts      = 30  :0.02:40.0 ;
    dens_labs  = 18  :0.25:32 ;
    dens_labs2 = 18  :0.50:32 ;
    [Y,X]      = meshgrid(theta,salts) ;
    sigma0     = this_case.static_fns.rho(X,Y,0) - 1000 ;
    iceT       = this_case.static_fns.T_freezing(salts,0) ;
    % Cut out densities below the freezing point.
    for ss = 1:length(salts)
        inds2            = find(Y(ss,:)<iceT(ss)) ;
        sigma0(ss,inds2) = NaN(length(inds2),1) ;
    end % ss
    this_case.base_case.app.replot_TS_diagram = false ;         % This is necessary if the freezing temperature or the EOS changes in the app.
end % if

% Plot TS diagram
[cs,h] = contour(X,Y,sigma0,dens_labs,'k') ;
hold on
plot(salts,iceT,'k-','linewidth',2)
clabel(cs,h,dens_labs2,'labelspacing',0.5*288,'fontsize',7)
axis([min(salts) max(salts) min(theta) max(theta)]) ;

% Plot results from theory
if(this_case.plot_options.diag_level > 0)
    objects.p_theory2  = plot(this_case.theory.S3s,this_case.theory.U_2_max.vals,'color',this_case.theory.U_2_max.color) ;
    objects.p_theoryi  = plot(this_case.theory.S3s,this_case.theory.U_i_max.vals,'color',this_case.theory.U_i_max.color) ;
    if(this_case.plot_options.diag_level > 1)
        objects.p_theoryib = plot(this_case.theory.S3s,this_case.theory.U_i_min.vals,'color',this_case.theory.U_i_min.color) ;
        objects.p_theory2b = plot(this_case.theory.S3s,this_case.theory.U_2_min.vals,'color',this_case.theory.U_2_min.color) ;
        objects.p_theory   = plot(this_case.theory.S3s,this_case.theory.denom0.vals ,'color',this_case.theory.denom0.color) ;
        objects.p_theory3b = plot(this_case.theory.S3s,this_case.theory.U_3_min.vals,'color',this_case.theory.U_3_min.color) ;
    end % if
end % if

if(this_case.plot_options.diag_level > 0 || this_case.base_case.app_mode)
    
    if(isfield(this_case.theory.patch_data,'S1'))
        objects.p1 = patch(this_case.theory.patch_data.S1,this_case.theory.patch_data.T1,this_case.plot_options.OW_col) ;
        set(objects.p1,'edgecolor','none','FaceAlpha',this_case.plot_options.transp_factor) ;
        uistack(objects.p1,'bottom') ;
    end % if
    if(isfield(this_case.theory.patch_data,'S2'))
        objects.p2 = patch(this_case.theory.patch_data.S2,this_case.theory.patch_data.T2,this_case.plot_options.OW_col) ;
        set(objects.p2,'edgecolor','none','FaceAlpha',this_case.plot_options.transp_factor) ;
        uistack(objects.p2,'bottom') ;
    end % if
    if(isfield(this_case.theory.patch_data,'S3'))
        objects.p3 = patch(this_case.theory.patch_data.S3,this_case.theory.patch_data.T3,this_case.plot_options.aW_col) ;
        set(objects.p3,'edgecolor','none','FaceAlpha',this_case.plot_options.transp_factor) ;
        uistack(objects.p3,'bottom') ;
    end % if
    
end % if

% Plot water masses
objects.p_AW  = scatter(this_case.base_case.S_1,this_case.base_case.T_1 ,80,'o','filled','Markerfacecolor',this_case.plot_options.AW_col) ;
objects.p_PW  = scatter(this_case.dynamic_parameters.S_2,this_case.dynamic_parameters.Tf2,80,'o','filled','Markerfacecolor',this_case.plot_options.PW_col) ;
if(isfield(this_case,'vals'))
    objects.p_OW_vals = scatter(cell2mat(this_case.vals.S_3),cell2mat(this_case.vals.T_3),20,'o','filled','Markerfacecolor',this_case.plot_options.OW_col) ;
    objects.p_SW_vals = scatter(cell2mat(this_case.vals.S_s),cell2mat(this_case.vals.T_f),20,'o','filled','Markerfacecolor',this_case.plot_options.SW_col) ;
end % if
objects.p_SW  = scatter(this_case.sens.S_scan,this_case.sens.Tfs,80,'o','filled','Markerfacecolor',this_case.plot_options.SW_col,'MarkerEdgecolor','w') ;
objects.p_aW  = scatter(this_case.dynamic_parameters.S_a,this_case.dynamic_parameters.T_a,80,'x','Markerfacecolor',this_case.plot_options.aW_col,'MarkerEdgecolor',this_case.plot_options.aW_col) ;
objects.p_OW  = scatter(this_case.sens.S_3can,this_case.sens.T_3can ,80,'o','filled','Markerfacecolor',this_case.plot_options.OW_col,'MarkerEdgecolor','w') ;

xlabel('Salinity [g/kg]')
ylabel('Temperature [^{\rm o}C]')
set(gca,'XLim',this_case.base_case.app.S1_lims+[-1 1]) ;
set(gca,'YLim',this_case.base_case.app.T1_lims+[-0.5 0.5]) ;
grid on

switch(this_case.plot_options.diag_level)
    case 0      % Production
        l = legend([objects.p_AW,objects.p_PW,objects.p_OW,objects.p_SW,objects.p_aW],...
            'AW','PW','OW','SW','aW') ;
    case 1      % Partial diagnostics
        l = legend([objects.p_AW,objects.p_PW,objects.p_OW,objects.p_SW,objects.p_aW,objects.p_theory2,objects.p_theoryi],...
            'AW','PW','OW','SW','aW',this_case.theory.U_2_max.legend,this_case.theory.U_i_max.legend) ;
    case 2      % Full diagnostics
        l = legend([objects.p_AW,objects.p_PW,objects.p_OW,objects.p_SW,objects.p_aW,objects.p_theory,objects.p_theory2,objects.p_theory2b,objects.p_theoryi,objects.p_theoryib,objects.p_theory3b],...
            'AW','PW','OW','SW','aW',this_case.theory.denom0.legend,this_case.theory.U_2_max.legend,this_case.theory.U_2_min.legend,this_case.theory.U_i_max.legend,this_case.theory.U_i_min.legend,this_case.theory.U_3_min.legend) ;
end % switch
set(l,'location','eastoutside','interpreter','latex') ;
set(gca,'Box','on') ;
title([this_case.base_case.exp_str,' ',this_case.base_case.tit_str],'Interpreter','latex','FontSize',14) ;

if(~this_case.base_case.app_mode)
    pos = get(gca,'Position') ;
    set(gca,'Position',[0.3 pos(2) 0.35 pos(4)]) ;
else
    annotation('textbox',[0.725 0.12 0.25 0.14],'String',{'Patches show OW theory constraints:','Brown=budget equations','Green=entrainment model'},'EdgeColor','k','backgroundcolor','w') ;
end % if

end