clc();
clear();
clf();

xCanon = [0;0.1;0.1;0];
yCanon = [0;0;0.1;0.1];

% Left 
xL = xCanon(:,ones(1,6));
yL = [yCanon,yCanon*3+0.1,yCanon*3+0.3*1+0.1,yCanon*3+0.3*2+0.1,yCanon*3+0.3*3+0.1,yCanon+0.3*4+0.1];

%   Right;
xR = xL + 0.3;
yR = yL;

%   Bridges
xBrB = xCanon*2 + 0.1;
yBrB = yCanon;
xBrT = xBrB;
yBrT = yCanon + 0.3*4 + 0.1;

%   Global
xG = [xL,xBrB,xR,xBrT];
yG = [yL(:,end:-1:1),yBrB,yR,yBrT];

Qdots = [1,2,4,6];
spI   = {1,3,4,6};

minT = Inf;
maxT = 0;

var = 'state';
n = 1;

for k = 4:numel(Qdots)
    q = load(['qSat_',num2str(Qdots(k),'%02G'),'kW.mat']);
    minT = min([minT;q.(var){end}.x(:,end)]);
    maxT = max([maxT;q.(var){end}.x(:,end)]);
end
q = load(['qSat_',num2str(Qdots(end),'%02G'),'kW.mat']);
minT = min([minT;q.(var){end}.x(:,end)]);
maxT = max([maxT;q.(var){end}.x(:,end)]);
    
for k = 1:numel(Qdots(1:end-1))
    power = num2str(Qdots(k),'%02G');
    q = load(['qSat_',power,'kW.mat']);
    power = ['$',num2str(Qdots(k),'%2G'),'$ kW'];
    subplot(2,3,spI{k});
    
    x = q.(var){end}.x(:,end);
    p = patch(xG,yG,(x<0)*0+ (x>=-eps()).*x);
    p.Parent.TickLabelInterpreter = 'latex';
    title(power,'interpreter','latex','FontWeight','normal');
    caxis([0,maxT]);
    axis('tight');
    axis('equal');
    box('on');
end
power = num2str(Qdots(end),'%02G');
q = load(['qSat_',power,'kW.mat']);
    power = ['$',num2str(Qdots(end),'%2G'),'$ kW'];
subplot(2,3,spI{end});
x = q.(var){end}.x(:,end);
    p = patch(xG,yG,(x<0)*0+ (x>=-eps()).*x);
p.Parent.TickLabelInterpreter = 'latex';
title(power,'interpreter','latex','FontWeight','normal');
axis('tight');
axis('equal');
box('on');


subplot(2,3,6);
c = colorbar('AxisLocation','out','Location','north');
c.Limits = [0,maxT];
c.Label.Interpreter = 'latex';
switch(var)
    case('rho')
        c.Label.String = '$\rho$ [kg/m${}^{3}$]';
    case('T')
        c.Label.String = '$T$ [K]';
    case('P')
        c.Label.String = '$P$ [bar]';
    case('state')
        c.Label.String = '$x$ [-]';
end


