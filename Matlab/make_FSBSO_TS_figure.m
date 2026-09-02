% Read and Process Fram Strait Data
% From OC3D original code make_Fram_Strait_figure.m
% twnh Dec '15, May '20

% Housekeeping
close all
clear 
clear global
more off
fprintf(1,'\n\n make_FSBSO_TS_figure.m \n Script to plot figure of Fram Strait and BSO hydrography.\n twnh Dec ''15, May ''16\n\n') ;
base_FS = 10 ;

% Read raw FS data
filename = '../data/processed/ARK-XVIII_1_phys_oce.tab' ;
fprintf(1,' Reading Fram Strait section [%s]...',filename) ;
[Station, DateTime, lats, lons, press, temps, salts] = read_Fram_Strait_section(filename) ;
fprintf(1,'done.\n') ;

% Eliminate stations off the primary section.
target_lat     = 78.83 ;
lat_thresh     = 0.02 ;
inds           = find(abs(lats - target_lat) > lat_thresh) ;
DateTime(inds) = [] ;
lats(inds)     = [] ;
lons(inds)     = [] ;
press(inds)    = [] ;
temps(inds)    = [] ;
salts(inds)    = [] ;

% Extract individual stations
fprintf(1,' Extracting stations...') ;
start_inds = [1; 1+find(diff(press) < 0)] ;                   % Find the start (surface) of each profile.
  end_inds = [start_inds(2:end)-1; length(press)]  ;          % Find the end (deepest point) of each profile.

% Reorder data so that section has monotonic increasing longitude.
[~,inds]   = sortrows(lons(start_inds)) ;
inds       = flipud(inds) ;
start_inds = start_inds(inds) ;
  end_inds =   end_inds(inds) ;
  
% Cut stations that mess up plot
start_inds(8) = [] ;
end_inds(8)   = [] ;
no_stations   = length(start_inds) ;  
fprintf(1,'done. Found [%d] stations.\n',no_stations) ;

% Convert to absolute salts
SA    = gsw_SA_from_SP(salts, press, lons, lats) ;

% Convert to conservative temperature
Theta = gsw_CT_from_t(SA,temps,press) ;

% Make a Theta-S_A diagram.
fprintf(1,'\n Making Theta-S_A plot.\n') ;
figure
set(gcf,'PaperPosition', [0 0 7 5]);
set(gcf,'Papersize',[7 5]) ;
axes('position',[0.1 0.1 0.8 0.8]) ;

lgray = 0.8.*[1 1 1] ;
dgray = 0.6.*[1 1 1] ;
bgray = 0.4.*[1 1 1] ;
hold on
dens_labs  =  18.0:0.25:32 ;
dens_labs2 =  18.0:0.50:32 ;
tmp_theta  =  -2.5:0.05:8 ;
tmp_salts  =  32.0:0.02:35.8 ;
[Y,X]      = meshgrid(tmp_theta,tmp_salts) ;
sigma0     = gsw_rho(X,Y,0)   - 1000 ;
iceT       = gsw_CT_freezing(tmp_salts,0,0) ;

% Cut out densities below the freezing point.
for ss = 1:length(tmp_salts)
   inds   = find(Y(ss,:)<iceT(ss)) ;
   sigma0(ss,inds) = NaN(length(inds),1) ;
end % ss
[cs,h] = contour(X,Y,sigma0,dens_labs,'color',dgray) ;
plot(tmp_salts,iceT,'-','linewidth',1,'color',dgray)
clabel(cs,h,dens_labs2,'labelspacing',0.5*288,'fontsize',base_FS-2,'color',dgray)
axis([min(tmp_salts) max(tmp_salts) min(tmp_theta) max(tmp_theta)]) ;
axis square
xlabel('Salinity [g/kg]','interpreter','tex','fontsize',base_FS+5) ;
ylabel('Temperature [^{\rm o}C]','interpreter','tex','fontsize',base_FS+5) ;
grid on
orient tall
set(gca,'box','on') ;
wysiwyg
plot(SA,Theta,'.','color',lgray) ;

% Annotations
text(32.90,-0.5,'Polar Water'   ,'Fontsize',base_FS+2,'backgroundcolor','w') ;
text(35.45, 2.8,'Atlantic Water','Fontsize',base_FS+2,'rotation',90,'backgroundcolor','w') ;
text(35.30,-0.5,'Overflow W.','Fontsize',base_FS+2,'rotation',80,'backgroundcolor','w') ;

%% Read BSO data
% JCR cruise in August 2017. Not ideal because the stations aren't very close and
% the data isn't very dense in depth. Struggled to find a better section at
% NODC.
% https://www.nodc.noaa.gov/OC5/SELECT/allcruises/GB013327.html
files = {'../data/processed/ocldb1590767725.30864_CTD.nc'} ;

for tt = 1:size(files,1)
   this_file = files{tt} ;
   fprintf(1,' Reading [%s].\n',this_file) ;
   Lat    = ncread(this_file,'lat');
   Lon    = ncread(this_file,'lon');
   temps  = ncread(this_file,'Temperature');
   salts  = ncread(this_file,'Salinity');
   depths = ncread(this_file,'z');
   wh     = find(Lon < 0);
   Lon(wh) = (360+Lon(wh));
   fprintf(1,'done.\n') ;
end % tt

% Cut deep data
inds = find(depths < 600) ;

press = gsw_p_from_z(-depths(inds),mean(Lat)) ;
% Convert to absolute salts
SA    = gsw_SA_from_SP(salts(inds), press, mean(Lon), mean(Lat)) ;

% Convert to conservative temperature
Theta = gsw_CT_from_t(SA,temps(inds),press) ;

plot(SA,Theta,'.','color',bgray) ;

% Wrap up
print -dpdf FSBSO_data2.pdf