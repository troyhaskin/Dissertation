clc();
% clear();


data = load('.\qSat_390K-Transient.mat');
q0   = [2:4:14,16]*1E3;

[xmax,imax] = cellfun(@(c) max(c.x,[],1),data.state,'UniformOutput',false);
voidmax = cellfun(@(c) VoidFractionFromDensity(c.rho(14,:),c.rhoL(14,:),c.rhoG(14,:),1),data.state,'UniformOutput',false);




figure(1);
clf();
ax = axes();
hold('on');
h = cellfun(@(c1,c2) plot(ax,c1,100*c2,'LineWidth',2),t,xmax,'UniformOutput',false);
ax.XAxis.Label.Interpreter = 'none';
ax.YAxis.Label.Interpreter = 'none';
leg = legend(strcat('\SI{',num2str(q0(:)/1E3),'}{\kW}'),'Location','NorthEast');
leg.Interpreter = 'none';
xlabel('$t$ \si{\second}');
ylabel('$x$ \si{\percent}');
ylim([0,0.094]);
grid('on');
box('on');
hold('off');

figure(2);
clf();
ax = axes();
hold('on');
h = cellfun(@(c1,c2) plot(ax,c1,100*c2,'LineWidth',2),t,voidmax,'UniformOutput',false);
ax.XAxis.Label.Interpreter = 'none';
ax.YAxis.Label.Interpreter = 'none';
leg = legend(strcat('\SI{',num2str(q0(:)/1E3),'}{\kW}'),'Location','NorthEast');
leg.Interpreter = 'none';
xlabel('$t$ \si{\second}');
ylabel('$\alpha$ \si{\percent}');
ylim([0,47]);
grid('on');
box('on');
hold('off');

