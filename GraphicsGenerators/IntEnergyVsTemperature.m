clc();
clear();

isoChor = 210;
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
h(3) = semilogx(rhoSat([isoChor,isoChor]),[0,10],'Color',Color,'LineWidth',1.0);

% Plot a o at the saturation state on the isochore
Color = GetColor('purple');
h(4) = semilogx(rhoSat(isoChor),iSat(isoChor)/1E6,'s','Color',Color,'MarkerFaceColor',Color,'MarkerSize',MarkerSize);

% Plot a o inside the vapor dome on the isochore
Color = GetColor('blue');
h(5) = semilogx(rhoSat(isoChor),0.5*iSat(isoChor)/1E6,'o','Color',Color,'MarkerFaceColor',Color,'MarkerSize',MarkerSize);

% Plot a o outside the vapor dome on the isochore
Color = GetColor('green');
h(6) = semilogx(rhoSat(isoChor),1.2*iSat(isoChor)/1E6,'o','Color',Color,'MarkerFaceColor',Color,'MarkerSize',MarkerSize);

% Plot Stuff
axis([3E-3,1.5E3,-0.1,3.7]);
xLabelHandle = xlabel('Density [\si{\kg\per\meter\cubed}]','Interpreter','none');
ylabel('Internal Energy [\si{\mega\joule\per\kg}]','Interpreter','none');
g = legend(h,'Saturation line'  ,...
             'Triple line'      ,...
             'Isochore'         ,...
             'Saturation state' ,...
             'Two-Phase state'  ,...
             'One-Phase state'  ,...
             'Location','NorthWest');
         
% set(g,'FontSize',get(xLabelHandle,'FontSize'));
MakePNG('InteEnergyVsDensity');
CropOuterWhite(WhereWeAre(),'InteEnergyVsDensity');



