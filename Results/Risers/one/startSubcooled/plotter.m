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

Qdots = 2.^[0,1,2,3,4];
spI   = {1,5,6,8,10};

minT = Inf;
maxT = 0;

var = 'rho';
n = 1;

for k = 1:numel(Qdots)
    q = load(['q_',num2str(Qdots(k),'%02G'),'kW.mat']);
    minT = min([minT;q.(var){end,n}(:,end)]);
    maxT = max([maxT;q.(var){end,n}(:,end)]);
end
q = load(['q_',num2str(Qdots(end),'%02G'),'kW.mat']);
minT = min([minT;q.(var){end,n}(:,end)]);
maxT = max([maxT;q.(var){end,n}(:,end)]);
    
for k = 1:numel(Qdots(1:end-1))
    power = num2str(Qdots(k),'%02G');
    q = load(['q_',power,'kW.mat']);
    power = ['$',num2str(Qdots(k),'%2G'),'$ kW'];
    subplot(2,5,spI{k});
    p = patch(xG,yG,q.(var){end,n}(:,end));
    p.Parent.TickLabelInterpreter = 'latex';
    title(power,'interpreter','latex','FontWeight','normal');
    caxis([minT,maxT]);
    axis('tight');
    axis('equal');
    box('on');
end
power = num2str(Qdots(end),'%02G');
q = load(['q_',power,'kW.mat']);
    power = ['$',num2str(Qdots(end),'%2G'),'$ kW'];
subplot(2,5,spI{end});
p = patch(xG,yG,q.(var){end,n}(:,end));
p.Parent.TickLabelInterpreter = 'latex';
title(power,'interpreter','latex','FontWeight','normal');
axis('tight');
axis('equal');
box('on');


subplot(2,5,8);
c = colorbar('AxisLocation','out','Location','north');
c.Limits = [minT,maxT];
c.Label.Interpreter = 'latex';
c.Label.String = '$\rho$ [kg/m${}^{3}$]';


% 
% for k = 1:numel(Qdots)
%     q = load(['q_',num2str(Qdots(k),'%02G'),'kW.mat']);
%     minT = min([minT;q.(var){end,n}(:,end)]);
%     maxT = max([maxT;q.(var){end,n}(:,end)]);
% end
% q = load(['q_',num2str(Qdots(end),'%02G'),'kW.mat']);
% minT = min([minT;q.(var){end,n}(:,end)]);
% maxT = max([maxT;q.(var){end,n}(:,end)]);
%     
% for k = 1:numel(Qdots(1:end-1))
%     power = num2str(Qdots(k),'%02G');
%     q = load(['q_',power,'kW.mat']);
%     power = ['$',num2str(Qdots(k),'%2G'),'$ kW'];
%     subplot(2,5,spI{k});
%     p = patch(xG,yG,q.(var){end,n}(:,end));
%     p.Parent.TickLabelInterpreter = 'latex';
%     title(power,'interpreter','latex','FontWeight','normal');
%     caxis([minT,maxT]);
%     axis('tight');
%     axis('equal');
%     box('on');
% end
% power = num2str(Qdots(end),'%02G');
% q = load(['q_',power,'kW.mat']);
%     power = ['$',num2str(Qdots(end),'%2G'),'$ kW'];
% subplot(2,5,spI{end});
% p = patch(xG,yG,q.(var){end,n}(:,end));
% p.Parent.TickLabelInterpreter = 'latex';
% title(power,'interpreter','latex','FontWeight','normal');
% axis('tight');
% axis('equal');
% box('on');
% 
% 
% subplot(2,5,8);
% c = colorbar('AxisLocation','out','Location','north');
% c.Limits = [minT,maxT];
% c.Label.Interpreter = 'latex';
% c.Label.String = '$\rho$ [kg/m${}^{3}$]';

