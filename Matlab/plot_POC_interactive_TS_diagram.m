function plot_POC_interactive_TS_diagram(this_case)

dragging = [];
orPos = [];
this_case.plot_options.diag_level = 2 ;

%% Make Interactive TS diagram
figure
objects = plot_POC_TS_diagram(this_case) ;
objects.p_OW.ButtonDownFcn = @dragObject ;
objects.p_PW.ButtonDownFcn = @dragObject ;
objects.p_aW.ButtonDownFcn = @dragObject ;

set(gcf,'WindowButtonUpFcn',@dropObject,'units','normalized','WindowButtonMotionFcn',{@moveObject,this_case,objects});
keyboard

%% Interactive functions

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
            dragging = [];
        end
    end
    function moveObject(hObject,~,this_case,objects)
        if ~isempty(dragging)
            newPos = get(gca,'CurrentPoint');
            posDiff = newPos - orPos;
            orPos = newPos;
            new_S = get(dragging,'Xdata') + posDiff(1,1) ;
            new_T = get(dragging,'Ydata') + posDiff(1,2) ;
            
            switch(hObject.CurrentObject.DisplayName)
                case 'OW'
                    aWS   = objects.p_aW.XData ;
                    aWT   = objects.p_aW.YData ;
                    
                    new_rho = gsw_rho(new_S,new_T,this_case.geophysical_parameters.pressure) ;
                    if(new_rho < this_case.dynamic_parameters.rh1)      % OW denisty exceeds AW density
                        new_S = gsw_SA_from_rho(this_case.dynamic_parameters.rh1,new_T,this_case.geophysical_parameters.pressure) ;
                    end
                    
                    new_T = max(new_T,gsw_CT_freezing(new_S,this_case.geophysical_parameters.pressure)) ;
                    new_T = min(new_T,aWT) ;
                    
                    S_s = @(Phi) (new_S - Phi*aWS)/(1-Phi) ;
                    Tfs = @(Phi) (new_T - Phi*aWT)/(1-Phi) ;
                    bnd_fn = @(Phi) (Tfs(Phi) - gsw_CT_freezing(S_s(Phi),this_case.geophysical_parameters.pressure))^2 ;
                    new_Phi = fminbnd(bnd_fn,0,1) 
                    
                    new_S_s = S_s(new_Phi) ;
                    new_Tfs = Tfs(new_Phi) ;
                    set(objects.p_SW,'XData',new_S_s) ;
                    set(objects.p_SW,'YData',new_Tfs) ;
                    
                    % Test consistency with entrainment formula.
                    aW_rho = gsw_rho(aWS    ,aWT    ,this_case.geophysical_parameters.pressure) ;
                    SW_rho = gsw_rho(new_S_s,new_Tfs,this_case.geophysical_parameters.pressure) ;
                    
                    Phi_theory = 1 - this_case.static_parameters.gam*((this_case.dynamic_parameters.U_1)^(1/3))/((SW_rho - aW_rho)^(2/3)) 
                    if(abs(Phi_theory - new_Phi) < this_case.static_parameters.DPhi_limit)
                        plot(new_S,new_T,'r*')
                    end % if
                    
                case 'PW'
                    new_S = min(new_S,this_case.sens.S_2) ;
                    new_T = gsw_CT_freezing(new_S,this_case.geophysical_parameters.pressure) ;
                    
                    phi = (objects.p_aW.XData - this_case.dynamic_parameters.S_1)/(objects.p_PW.XData - this_case.dynamic_parameters.S_1) ;
                    new_aS = new_S*phi + this_case.dynamic_parameters.S_1*(1-phi) ;
                    new_aT = new_T*phi + this_case.dynamic_parameters.T_1*(1-phi) ;
                    set(objects.p_aW,'XData',new_aS) ;
                    set(objects.p_aW,'YData',new_aT) ;
                    
                case 'aW'
                    new_T = min(new_T,this_case.dynamic_parameters.T_1) ;
                    new_T = max(new_T,gsw_CT_freezing(this_case.sens.S_2,this_case.geophysical_parameters.pressure)) ;
                    new_S = min(new_S,gsw_SA_from_rho(this_case.dynamic_parameters.rh1,new_T,this_case.geophysical_parameters.pressure)) ;
                    new_T = max(new_T,gsw_CT_freezing(new_S,this_case.geophysical_parameters.pressure)) ;
                    
                    S_2 = @(phi) (new_S + (phi - 1)*this_case.dynamic_parameters.S_1)/phi ;
                    Tf2 = @(phi) (new_T + (phi - 1)*this_case.dynamic_parameters.T_1)/phi ;
                    bnd_fn = @(phi) (Tf2(phi) - gsw_CT_freezing(S_2(phi),this_case.geophysical_parameters.pressure))^2 ;
                    new_phi = fminbnd(bnd_fn,0,1) ;
                    
                    new_S_2 = S_2(new_phi) ;
                    new_Tf2 = Tf2(new_phi) ;
                    set(objects.p_PW,'XData',new_S_2) ;
                    set(objects.p_PW,'YData',new_Tf2) ;
                    
            end % switch
            set(dragging,'XData',new_S);
            set(dragging,'YData',new_T);
        end
    end

end