%% Polar Overturning Conceptual Model
%% Live script to run experiments in paper and produce the figures.
%% twnh Aug, Sep '19, Jan?May '20, Sep '20
% Must cd to the working directory to run the notebook. The notebook writes 
% pdf figures to the figures directory. Crop them to make publication-ready figures.
% 
% Housekeeping:

clear
clear plot_POC_TS_diagram  % Clears persistent variables from this function.
close all
%parpool('MJSProfile1',84) 
if(isempty(gcp))
    %parpool('MJSProfile1',84)
    parpool('local',16) 
end
%% 
% Exp 1: Northern hemisphere Fram Strait/Barents Sea Opening.

model1 = build_POC_model('NH_FSBSO_can') ;
plot_POC_solution(model1)
%% 
% Exp 2: Modified Northern hemisphere Fram Strait/Barents Sea Opening (higher 
% Q).

model2 = build_POC_model('NH_FSBSO_mod') ;
plot_POC_solution(model2) ;
%% 
% Make figure showing entrainment/shelf circulation tradeoff

plot_POC_tradeoff(model1,model2)
%% 
% Exp 3: Northern hemisphere Fram Strait/Barents Sea Opening with many Q 
% values.

model3 = build_POC_model('NH_FSBSO_range') ;
plot_POC_solutions1(model3,model1,model2) ;
%% 
% Exp 4: Northern hemisphere Fram Strait/Barents Sea Opening with wide range 
% of all parameters.

model4 = build_POC_model('NH_FSBSO_sens') ;
plot_POC_solutions2(model4,model3,model1,model2)
%% 
% Compare numerical results to theory

plot_POC_solutions_and_theory(model4)
%% 
% Exp 5: Northern hemisphere Fram Strait/Barents Sea Opening with various 
% S_2 values.

model5 = build_POC_model('NH_FSBSO_S2_range') ;
plot_POC_solutions3(model5,model1)
%% 
% Exp 6: Southern hemisphere.

model6 = build_POC_model('SH_can') ;
plot_POC_solution(model6)
%% 
% Exp V: Southern hemisphere with warm AW following reviewer 2.

modelV = build_POC_model('SH_mod') ;
plot_POC_solution(modelV)
%% 
% Exp W: Eldevik & Nilsen (2013) parameters

modelW = build_POC_model('NH_FSBSO_EN13') ;
plot_POC_solution(modelW)
%% 
% Exp X: Modified Northern hemisphere Fram Strait/Barents Sea Opening (intermediate 
% Q).

modelX = build_POC_model('NH_FSBSO_int') ;
plot_POC_solution(modelX)
%% 
% Exp Y: Modified Northern hemisphere Fram Strait/Barents Sea Opening (small 
% S2).

modelY = build_POC_model('NH_FSBSO_small_S2') ;
plot_POC_solution(modelY)
%% 
% Exp Z: Southern hemisphere with many Q values

modelZ = build_POC_model('SH_range') ;
plot_POC_solutions1(modelZ,model6,model6)