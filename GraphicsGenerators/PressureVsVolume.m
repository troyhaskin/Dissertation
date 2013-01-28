clc;
clear('all');

T   = [linspace(TriplePointTemperature,600,200),linspace(601,CriticalTemperature,200)]';
rho = GeometricSpace(0.05,990,1.001,10);


[Psat,rhol,rhog] = SaturationStateGivenTsat(T);
rhoSat = [rhog;rhol(end:-1:1)];
Psat   = [Psat;Psat(end:-1:1)];
P300   = Pressure(rhoSat,300);
P500   = Pressure(rhoSat,500);
P647   = Pressure(rhoSat,T(end));
P700   = Pressure(rhoSat,700);

FontName = 'CMU Serif';

loglog(1./rhoSat,Psat/1E6,'Color',0.5*[1,1,1],'LineWidth',2)
hold('on');
h(1) = loglog(1./rhoSat,P300/1E6,'Color',GetColor('red'),'LineWidth',1.5);
h(2) = loglog(1./rhoSat,P500/1E6,'Color',GetColor('blue'),'LineWidth',1.5);
h(3) = loglog(1./rhoSat,P647/1E6,'Color',GetColor('green'),'LineWidth',1.5);
h(4) = loglog(1./rhoSat,P700/1E6,'Color',GetColor('purple'),'LineWidth',1.5);
xLabelHandle = xlabel('Specific Volume [m^3/kg]','FontName',FontName);
ylabel('Pressure [MPa]','FontName',FontName);
grid('on');
box('on');
axis([7E-4,2.5E2,1E-3,1E3]);
set(gca,'FontName',FontName,'FontSize',8);
g = legend(h,'T = 300 K','T = 500 K','T = 647 K','T = 700 K','Location',[0.63,0.75,0.30,0.1]);
% set(g,'FontSize',get(xLabelHandle,'FontSize'));
hold('off');
MakePNG('PressureVsVolume');
CropOuterWhite(WhereWeAre(),'PressureVsVolume');

[rho,T] = meshgrid(rho,T);
