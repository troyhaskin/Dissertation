% clc();
clear();

hem = testHEM();

%   Adjust evolver time stuff
dt = 1e-5  ;
hem.set('evolver','time.span'         , [0,5e-2]  );
hem.set('evolver','time.step.maximum' , 1                               );
hem.set('evolver','time.step.minimum' , 1E-7                            );
hem.set('evolver','time.step.goal'    , dt                              );
hem.set('evolver','saveRate'          , dt                              );
tic;
hem.run();
[q,Dq,t] = hem.results();
tElasped = toc;


%   Pull
mass      = q( (1:14) + 0   ,:) ;
energy    = q( (1:14) + 14  ,:) ;
volume    = q( (1:14) + 14*2,:) ;
momentum  = q( (1:14) + 14*3,:) ;
Dmass     = Dq( (1:14) + 0   ,:);
Denergy   = Dq( (1:14) + 14  ,:);
Dvolume   = Dq( (1:14) + 14*2,:);
Dmomentum = Dq( (1:14) + 14*3,:);
rho       = mass   ./ volume    ;
i         = energy ./ mass      ;
T         = Temperature(rho,i)  ;
P         = Pressure(rho,T)     ;


%%
%   Plot

mask   = 1:5000;
tinter = linspace(t(mask(1)),t(mask(end)),1E3);
mask2  = [2800:200:3900,4000];
Pinter = IntrepidTwilight.ConvenientMeans.hermiteInterpolation(t(mask2),momentum(1,mask2),[nan(1,numel(mask2)-1),Dmomentum(1,mask2(end-[0]))],tinter);

plot(t(mask),momentum(1,mask),tinter,Pinter,t(mask2),momentum(1,mask2),'o');
ylim([0.05,0.3]);




