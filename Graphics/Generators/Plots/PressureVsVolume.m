clc;
clear();
clf();

T   = [linspace(TriplePointTemperature,600,200),linspace(601,CriticalTemperature,200)]';


[Psat,rhol,rhog] = SaturationStateGivenTemperature(T);
rhoSat = [rhog;rhol(end:-1:1)];
Psat   = [Psat;Psat(end:-1:1)];
P300   = Pressure(rhoSat,300);
P500   = Pressure(rhoSat,500);
P647   = Pressure(rhoSat,T(end));
P700   = Pressure(rhoSat,700);


colors = get(gca,'ColorOrder');

loglog(1./rhoSat,Psat/1E6,'Color',0.5*[1,1,1],'LineWidth',2)
hold('on');
h(1) = loglog(1./rhoSat,P300/1E6,'Color',colors(1,:),'LineWidth',1.5);
h(2) = loglog(1./rhoSat,P500/1E6,'Color',colors(2,:),'LineWidth',1.5);
h(3) = loglog(1./rhoSat,P647/1E6,'Color',colors(3,:),'LineWidth',1.5);
h(4) = loglog(1./rhoSat,P700/1E6,'Color',colors(4,:),'LineWidth',1.5);
xLabelHandle = xlabel('Specific Volume [\si{\meter\cubed\per\kg}]','Interpreter','none');
ylabel('Pressure [\si{\mega\pascal}]','Interpreter','none');
Tstr = strcat({'$T = \SI{'},num2str([300;500;647;700]),'}{\kelvin}$');
g = legend(h,Tstr,'Location','NorthEast','Interpreter','none');
axis([7E-4,2.5E2,1E-3,1E3]);
grid('on');
box('on');
hold('off');

matlab2tikz2('..\..\PressureVsVolume.tikz','width','5in');
