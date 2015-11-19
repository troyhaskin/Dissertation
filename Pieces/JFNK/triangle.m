clc();
clear();

n     = 1;
lr0   = [1,1.5]                     ;
theta = atan(lr0(2)/lr0(1))         ;
lAKm  = [0.5*norm(lr0,2)*cos(theta),zeros(n,1)];
lrm   = [lr0(1)-lAKm(:,1),lr0(2)-lAKm(:,2)]                  ;
gamma = acos((lrm(:,1).*lAKm(:,1) + lrm(:,2).*lAKm(:,2)) ./...
        (sqrt(lrm(:,1).^2+lrm(:,2).^2).* sqrt(lAKm(:,1).^2+lAKm(:,2).^2)));

Show(lr0')    
Show(lAKm')
Show(lrm')
Show(gamma*180/pi')