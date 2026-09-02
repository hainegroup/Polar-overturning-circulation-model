function wysiwyg
%WYSIWYG -- this function is called with no args and merely
%       changes the size of the figure on the screen to equal
%       the size of the figure that would be printed, 
%       according to the papersize attribute.  Use this function
%       to give a more accurate picture of what will be 
%       printed.
%       Dan(K) Braithwaite, Dept. of Hydrology U.of.A  11/93
 
unis = get(gcf,'units');
ppos = get(gcf,'paperposition');
set(gcf,'units',get(gcf,'paperunits'));
pos = get(gcf,'position');
pos(3:4) = ppos(3:4);
set(gcf,'position',pos);
set(gcf,'units',unis);
drawnow

% Move lower left hand corner so resized figure can be seen. Different screen resolutions will
% require different positions.
pos = get(gcf,'Position') ;
pos(1) = 380 ;
pos(2) = 10 ;
set(gcf,'Position',pos) ;
pos = get(gcf,'Position') ;
drawnow
