clc();
clear();

l = sqrt(1-1/4);
x = [-0.5,0.5,0.5,1,1,-1,-1,-0.5];
y = [0,0,3/4,3/4,1,1,3/4,3/4];
P = [x;y];

R = @(theta) [cos(theta),-sin(theta);sin(theta),cos(theta)];

RP = R(pi/2)*P;
fill(RP(1,:),RP(2,:),'b');
axis('square');
axis('equal');
axis([-1,1,-1,1]);