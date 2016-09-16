clc();
clear();
clf();

Nden = 1000;

% Triple line generation
Tt = TriplePointTemperature();
[rhotl,rhotg] = TriplePointDensities();
rhot = GeometricSpace(1.001*rhotg,0.999*rhotl,1.07,200);
it = InternalEnergy(rhot,Tt);
% Pt = Pressure(rhot,Tt);

% Saturation line generation
Tsat = linspace(TriplePointTemperature,CriticalTemperature,100)';
[Psat,rhol,rhog] = SaturationStateGivenTemperature(Tsat);
rhoSat = [rhog;rhol(end:-1:1)];
iSat = InternalEnergy(rhoSat,[Tsat;Tsat(end:-1:1)]);
%
[~,isoChor] = min(abs(rhoSat - 4));
MarkerSize = 10;


h = newplot();
colors = get(gca,'ColorOrder');
set(h,'NextPlot','add','XGrid','on','YGrid','on','Box','on','XScale','log');

% Plot the saturation line
h(1) = semilogx(rhoSat,iSat/1E6,'Color',0.3*[1,1,1],'LineWidth',0.9);

% Plot the triple line
h(2) = semilogx(rhot,it/1E6,'Color',0.3*[1,1,1],'LineStyle','--','LineWidth',0.9);

% Plot an arbitrary isochore line
Color = colors(1,:);
h(3) = semilogx(rhoSat([isoChor,isoChor]),[0,10],'Color',Color,'LineWidth',1);

% Plot a o outside the vapor dome on the isochore
Color = colors(3,:);
h(4) = semilogx(rhoSat(isoChor),1.2*iSat(isoChor)/1E6,'s','Color',Color,'MarkerFaceColor',Color,'MarkerSize',MarkerSize);

% Plot a o at the saturation state on the isochore
Color = colors(4,:);
h(5) = semilogx(rhoSat(isoChor),iSat(isoChor)/1E6,'s','Color',Color,'MarkerFaceColor',Color,'MarkerSize',MarkerSize);

% Plot a o inside the vapor dome on the isochore
Color = colors(2,:);
h(6) = semilogx(rhoSat(isoChor),0.5*iSat(isoChor)/1E6,'s','Color',Color,'MarkerFaceColor',Color,'MarkerSize',MarkerSize);


% Plot Stuff
axis([3E-3,1.5E3,-0.1,4]);
xLabelHandle = xlabel('Density [\si{\kg\per\meter\cubed}]','Interpreter','none');
ylabel('Internal Energy [\si{\mega\joule\per\kg}]','Interpreter','none');
g = legend(h,'Saturation line'  ,...
             'Triple line'      ,...
             'Isochore'         ,...
             'One-Phase state'  ,...
             'Saturated'        ,...
             'Two-Phase state'  ,...
             'Location','NorthWest');
         

matlab2tikz2('..\..\IntEnergyVsTemperature.tikz','width','5in');

