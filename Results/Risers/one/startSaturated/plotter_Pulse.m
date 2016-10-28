clc();
clear();
clf();


data = load('.\qSat_372K-01kW-PointPulseInsertions.mat');

figure(1);
hold('on');
cellfun(@(t,s) plot(t-100,max(s.x,[],1)*100,'LineWidth',1.5),data.t,data.state);
hl = legend(strcat('\SI{',num2str([1,2,5,10,15].'),'}{\kW}'));
hl.Interpreter = 'none';
xlim([0,10]);
xlabel('$t$ [\si{\second}]','Interpreter','none');
ylabel('$x$ [$-$]','Interpreter','none');
box('on');
grid('on');
hold('off');

figure(2);
hold('on');
cellfun(@(t,s) plot(t-100,VoidFractionFromDensity(s.rho(end,:),s.rhoL(end,:),s.rhoG(end,:),1)*100,'LineWidth',1.5),data.t,data.state);
hl = legend(strcat('\SI{',num2str([1,2,5,10,15].'),'}{\kW}'));
hl.Interpreter = 'none';
xlim([0,10]);
xlabel('$t$ [\si{\second}]','Interpreter','none');
ylabel('$\alpha$ [$-$]','Interpreter','none');
box('on');
grid('on');
hold('off');

