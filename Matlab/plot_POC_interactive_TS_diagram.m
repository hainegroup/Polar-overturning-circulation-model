% Interactive TS diagram.
% twnh May '20

function plot_POC_interactive_TS_diagram(this_case)


dragging = [];
orPos = [];
this_case.plot_options.diag_level  = 2 ;

%% Make Interactive TS diagram
close all
figure
objects = plot_POC_TS_diagram(this_case) ;
this_case.dynamic_parameters.original_S_2 = this_case.dynamic_parameters.S_2 ;
objects.p_AW.ButtonDownFcn = @dragObject ;
objects.p_PW.ButtonDownFcn = @dragObject ;
objects.p_aW.ButtonDownFcn = @dragObject ;
delete(objects.p_theory_S_star) ;
delete(objects.p_OW) ;
delete(objects.p_SW) ;
delete(objects.p_OW_vals) ;
delete(objects.p_SW_vals) ;
set(get(gca,'Title'),'String','Interactive Theory Constraints') ;
set(gcf,'WindowButtonUpFcn',@dropObject,'units','normalized','WindowButtonMotionFcn',@moveObject,'Visible','on');

%% Interactive functions
% These are NESTED functions, so they can see the calling functions
% variables, which is important!

    function dragObject(hObject,~)
        dragging = hObject;
        orPos = get(gca,'CurrentPoint');
    end

    function dropObject(~,~)
        if ~isempty(dragging)
            newPos = get(gca,'CurrentPoint');
            posDiff = newPos - orPos;
            set(dragging,'XData',get(dragging,'Xdata') + posDiff(1,1));
            set(dragging,'YData',get(dragging,'Ydata') + posDiff(1,2));
            set(gca,'selected','on') ;
            dragging = [];
        end
    end

    function moveObject(hObject,~)
        if ~isempty(dragging)
            newPos = get(gca,'CurrentPoint');
            posDiff = newPos - orPos;
            orPos = newPos;
            new_S = get(dragging,'Xdata') + posDiff(1,1) ;
            new_T = get(dragging,'Ydata') + posDiff(1,2) ;
            
            old_phi = this_case.dynamic_parameters.phi ;
            new_phi = old_phi ;
            old_AWS = objects.p_AW.XData ;
            old_AWT = objects.p_AW.YData ;
            old_PWS = objects.p_PW.XData ;
            old_PWT = objects.p_PW.YData ;
            
            new_AWS = old_AWS ;
            new_AWT = old_AWT ;
            new_PWS = old_PWS ;
            new_PWT = old_PWT ;
            new_aWS = objects.p_aW.XData ;
            new_aWT = objects.p_aW.YData ;
            
            switch(hObject.CurrentObject.DisplayName)
                case 'AW'
                    new_AWS = new_S ;
                    new_AWT = max(new_T,gsw_CT_freezing(new_AWS,this_case.geophysical_parameters.pressure)) ;
                    
                    % Also update PW if necessary
                    new_PWS = min(new_PWS,compute_max_PW_salinity(new_AWT,new_AWS,this_case.geophysical_parameters,this_case.static_parameters)) ;
                    new_PWT = max(gsw_CT_freezing(new_PWS,this_case.geophysical_parameters.pressure)) ;
                    
                    % Also update aW
                    new_aWS = new_PWS*new_phi + new_AWS*(1-new_phi) ;
                    new_aWT = new_PWT*new_phi + new_AWT*(1-new_phi) ;
                    
                case 'PW'
                    new_PWS = min(new_S,this_case.dynamic_parameters.original_S_2) ;
                    new_PWT = max(gsw_CT_freezing(new_S,this_case.geophysical_parameters.pressure)) ;
                    
                    % Also update aW
                    new_aWS = new_PWS*new_phi + this_case.dynamic_parameters.S_1*(1-new_phi) ;
                    new_aWT = new_PWT*new_phi + this_case.dynamic_parameters.T_1*(1-new_phi) ;
                    
                case 'aW'
                    new_aWT = min(new_T,this_case.dynamic_parameters.T_1) ;
                    new_aWT = max(new_aWT,gsw_CT_freezing(this_case.dynamic_parameters.S_2,this_case.geophysical_parameters.pressure)) ;
                    new_aWS = min(new_S,gsw_SA_from_rho(this_case.dynamic_parameters.rh1,new_aWT,this_case.geophysical_parameters.pressure)) ;
                    new_aWT = max(new_aWT,gsw_CT_freezing(new_aWS,this_case.geophysical_parameters.pressure)) ;
                    
                    % Also update PW
                    S_2     = @(phi) (new_aWS + (phi - 1)*this_case.dynamic_parameters.S_1)/phi ;
                    Tf2     = @(phi) (new_aWT + (phi - 1)*this_case.dynamic_parameters.T_1)/phi ;
                    bnd_fn  = @(phi) (Tf2(phi) - gsw_CT_freezing(S_2(phi),this_case.geophysical_parameters.pressure))^2 ;
                    new_phi = fminbnd(bnd_fn,0,1) ;
                    
                    new_PWS = S_2(new_phi) ;
                    new_PWT = Tf2(new_phi) ;
                    
                    % Ensure that PW salinity can't exceed the maximum
                    % value.
                    new_PWS = min(new_PWS,this_case.dynamic_parameters.original_S_2) ;
                    
                    % Also update aW
                    new_aWS = new_PWS*new_phi + this_case.dynamic_parameters.S_1*(1-new_phi) ;
                    new_aWT = new_PWT*new_phi + this_case.dynamic_parameters.T_1*(1-new_phi) ;
                    
            end % switch
            
            % Update graphics objects and data structure
            this_case.dynamic_parameters.phi = new_phi ;
            this_case.dynamic_parameters.S_1 = new_AWS ;
            this_case.dynamic_parameters.T_1 = new_AWT ;
            this_case.dynamic_parameters.rh1 = gsw_rho(new_AWS,new_AWT,this_case.geophysical_parameters.pressure) ;
            this_case.dynamic_parameters.S_2 = new_PWS ;
            this_case.dynamic_parameters.Tf2 = new_PWT ;
            this_case.dynamic_parameters.S_a = new_aWS ;
            this_case.dynamic_parameters.T_a = new_aWT ;
            
            % Plot
            set(objects.p_AW,'XData',new_AWS) ;
            set(objects.p_AW,'YData',new_AWT) ;
            set(objects.p_PW,'XData',new_PWS) ;
            set(objects.p_PW,'YData',new_PWT) ;
            set(objects.p_aW,'XData',new_aWS) ;
            set(objects.p_aW,'YData',new_aWT) ;
            
            % Update lines and patches from theory:
            this_case.theory = compute_POC_theory(this_case.dynamic_parameters,this_case.static_parameters,this_case.geophysical_parameters) ;
            set(objects.p_theory ,'XData',this_case.theory.S3s) ;
            set(objects.p_theory ,'YData',this_case.theory.denom0.vals) ;
            set(objects.p_theory2,'XData',this_case.theory.S3s) ;
            set(objects.p_theory2,'YData',this_case.theory.numer0.U_2.vals) ;
            set(objects.p_theoryi,'XData',this_case.theory.S3s) ;
            set(objects.p_theoryi,'YData',this_case.theory.numer0.U_i.vals) ;
            
            if(isfield(this_case.theory.patch_data,'S1'))
                set(objects.p1       ,'XData',this_case.theory.patch_data.S1) ;
                set(objects.p1       ,'YData',this_case.theory.patch_data.T1) ;
            else
                set(objects.p1       ,'XData',[]) ;
                set(objects.p1       ,'YData',[]) ;
            end % if
            if(isfield(this_case.theory.patch_data,'S2'))
                set(objects.p2       ,'XData',this_case.theory.patch_data.S2) ;
                set(objects.p2       ,'YData',this_case.theory.patch_data.T2) ;
            else
                set(objects.p2       ,'XData',[]) ;
                set(objects.p2       ,'YData',[]) ;
            end % if
            if(isfield(this_case.theory.patch_data,'S3'))
                set(objects.p3       ,'XData',this_case.theory.patch_data.S3) ;
                set(objects.p3       ,'YData',this_case.theory.patch_data.T3) ;
            else
                set(objects.p3       ,'XData',[]) ;
                set(objects.p3       ,'YData',[]) ;
            end % if
        end
    end

end