clc;
clear('all');

x = linspace(-1,1,5E3);
y = @(x) (x+1/2).*(x+1/4).*(x-1/4).*(x-1/2);
t = linspace(0,2*pi,500);

% Lines
plot(x,x,'k',0,0.106,'bo','MarkerFaceColor','b','LineWidth',2,'MarkerSize',14);
axis([-1,1,-2,2]);
MakePNG('NoEquilibirum');
CropOuterWhite(WhereWeAre(),'NoEquilibirum');
%
plot(x,x*0,'k',0,0.095,'bo','MarkerFaceColor','b','LineWidth',2,'MarkerSize',14);
axis([-1,1,-2,2]);
MakePNG('NeutralEquilibirum');
CropOuterWhite(WhereWeAre(),'NeutralEquilibirum');


% Ellipses
plot(cos(t)/1.4,6*sin(t)-2,'k',0,4.29,'bo','MarkerFaceColor','b','LineWidth',2,'MarkerSize',14);
axis([-1,1,-2,10]);
MakePNG('Unstable');
CropOuterWhite(WhereWeAre(),'Unstable');
%
plot(cos(t)/1.4,6*sin(t)-2,'k',0,-7.78,'bo','MarkerFaceColor','b','LineWidth',2,'MarkerSize',14);
axis([-1,1,-10,-1]);
MakePNG('Stable');
CropOuterWhite(WhereWeAre(),'Stable');

% Quartics
plot(x,100*y(x),'k',0,1.86,'bo','MarkerFaceColor','b','LineWidth',2,'MarkerSize',14);
axis([-1,1,-2,10]);
MakePNG('LinearlyUnstableNonlinearlyStable');
CropOuterWhite(WhereWeAre(),'LinearlyUnstableNonlinearlyStable');
%
plot(x,-100*y(x),'k',0,-1.265,'bo','MarkerFaceColor','b','LineWidth',2,'MarkerSize',14);
axis([-1,1,-10,2]);
MakePNG('LinearlyStableNonlinearlyUnstable');
CropOuterWhite(WhereWeAre(),'LinearlyStableNonlinearlyUnstable');


