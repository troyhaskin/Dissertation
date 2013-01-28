% clc;
% clear('all');
% 
% Nden = 2E3;
% 
% % Triple line generation
% Tt = TriplePointTemperature();
% [rhotl,rhotg] = TriplePointDensities();
% rhot = GeometricSpace(1.001*rhotg,0.999*rhotl,1.02,Nden);
% it = InternalEnergy(rhot,Tt);
% % Pt = Pressure(rhot,Tt);
% 
% % Saturation line generation
% Tsat = linspace(TriplePointTemperature,CriticalTemperature,500)';
% [Psat,rhol,rhog] = SaturationStateGivenTsat(Tsat);
% rhoSat = [rhog;rhol(end:-1:1)];
% iSat = InternalEnergy(rhoSat,[Tsat;Tsat(end:-1:1)]);
% 

close('all')
Ichor = 210;
MarkerSize = 9;

FontName = 'CMU Serif';
h = newplot();
set(h,'NextPlot','add','XGrid','on','YGrid','on','Box','on','XScale','log');

% Plot the saturation line
h(1) = semilogx(rhoSat,iSat/1E6,'Color',0.3*[1,1,1],'LineWidth',1.0);

% Plot the triple line
h(2) = semilogx(rhot,it/1E6,'Color',0.3*[1,1,1],'LineStyle','--','LineWidth',1.0);

% Plot an arbitrary isochore line
Color = GetColor('purple');
h(3) = semilogx(rhoSat([Ichor,Ichor]),[0,10],'Color',Color,'LineWidth',1.0);

% Plot a o at the saturation state on the isochore
Color = GetColor('purple');
h(4) = semilogx(rhoSat(Ichor),iSat(Ichor)/1E6,'s','Color',Color,'MarkerFaceColor',Color,'MarkerSize',MarkerSize);

% Plot a o inside the vapor dome on the isochore
Color = GetColor('blue');
h(5) = semilogx(rhoSat(Ichor),0.5*iSat(Ichor)/1E6,'o','Color',Color,'MarkerFaceColor',Color,'MarkerSize',MarkerSize);

% Plot a o outside the vapor dome on the isochore
Color = GetColor('green');
h(6) = semilogx(rhoSat(Ichor),1.2*iSat(Ichor)/1E6,'o','Color',Color,'MarkerFaceColor',Color,'MarkerSize',MarkerSize);

% Plot Stuff
axis([3E-3,1.5E3,-0.1,3.7]);
xLabelHandle = xlabel('Density [kg/m^3]','FontName',FontName);
ylabel('Internal Energy [MJ/kg]','FontName',FontName);
set(gca,'FontName',FontName,'FontSize',8);

g = legend(h,'Saturation line'  ,...
             'Triple line'      ,...
             'Isochore'         ,...
             'Saturation state' ,...
             'Two-Phase state'  ,...
             'One-Phase state'  ,...
             'Location',[0.62,0.75,0.31,0.1]);
         
% set(g,'FontSize',get(xLabelHandle,'FontSize'));
MakePNG('InteEnergyVsDensity');
CropOuterWhite(WhereWeAre(),'InteEnergyVsDensity');



