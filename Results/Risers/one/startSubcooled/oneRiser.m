% clc();
clear();


%   Numbers for normals
r2 = cos(pi/4);
s2 = sqrt(2);

%   Simple geometry stuff
hcorn   = 0.1       ;
hhorz   = 0.2       ;
hvert   = 0.3       ;
Af      = 0.1^2     ;
volCorn = hcorn * Af;
volHorz = hhorz * Af;
volVert = hvert * Af;

%   Vertical / horizontal maps
vertMap = [-1;-1;-1;-1;-1; 0; 0;+1;+1;+1;+1;+1; 0; 0];
horzMap = [ 0; 0; 0; 0; 0;+1;+1; 0; 0; 0; 0; 0;-1;-1];

%   Path lengths / hydraulic diameters
L  = [   hcorn/2 + hvert/2   ;
         hvert*ones(3,1)     ;
         hcorn/2 + hvert/2   ];
L  = [L ; (hcorn + hhorz)*[1;1]/2 ; L ; (hcorn + hhorz)*[1;1]/2 ];
Dh = 0.1 * ones(14,1);

%   Define initial state
T0   = 300                      ; % [K]
P0   = 101325                   ; % [Pa]
x0   = 0.01                     ; % [-]
rho0 = Density(P0,T0)           ; % [kg/m^3]
i0   = InternalEnergy(rho0,T0)  ; % [J/kg]
v0   = 0.1                      ; % [m/s]

%   Get densities and energies assuming incompressible, hydrostatic pressures 
P     = P0 + 9.81*rho0*[0; hcorn/2 + hvert*(1:2:7).'/2 ; hvert*4 + hcorn ];
P     = [P ; P(end) ; P(end:-1:1) ; P(1)];
T     = P*0 + T0;
rho   = Density(P,T);
i     = InternalEnergy(rho,T);


%   Define extensive properties
volume       = [              volCorn           ;
                         volVert*ones(4,1)      ;
                    volCorn ; volHorz ; volCorn ;
                         volVert*ones(4,1)      ;
                         volCorn ; volHorz      ];
volumeMax    = volume                                                   ;
mass         = rho .* volume                                            ;
energy       = i   .* mass                                              ;
momentum     = mass .* v0                                               ;


%   Instantiate the simulation
hem = IntrepidTwilight.AdamantWave.HEM();

%  12 control volumes
hem.set('model','momentumCell.from', (1:14).');
hem.set('model','momentumCell.to'  , [2:14,1].');
%
%   momentum cell geometry
hem.set('model','momentumCell.directionX'         , horzMap         );
hem.set('model','momentumCell.directionY'         , vertMap         );
hem.set('model','momentumCell.volumeFractionFrom' , ones(14,1)/2    );
hem.set('model','momentumCell.volumeFractionTo'   , ones(14,1)/2    );
hem.set('model','momentumCell.flowArea'           , Af*ones(14,1)   );
hem.set('model','momentumCell.LoD'                , L./Dh           );
hem.set('model','momentumCell.source.friction'    , 0.1             );
%
%   Interfaces
hem.set('model','interface.up'  , (1:14).'  );
hem.set('model','interface.down', [2:14,1].');
hem.set('model','interface.normalX',       [  0 ;  0 ;  0 ;  0 ; +r2 ; +1 ; +r2 ;  0 ;  0 ;  0 ;  0 ; -r2 ; -1 ; -r2]);
hem.set('model','interface.normalY',       [ -1 ; -1 ; -1 ; -1 ; -r2 ;  0 ; +r2 ; +1 ; +1 ; +1 ; +1 ; +r2 ;  0 ; -r2]);
hem.set('model','interface.flowArea', Af * [ +1 ; +1 ; +1 ; +1 ; +s2 ; +1 ; +s2 ; +1 ; +1 ; +1 ; +1 ; +s2 ; +1 ; +s2]);

%   Sources
hem.set('model','controlVolume.source.mass'   , @(varargin) 0);
source = @(t) 1E3 * (1 + tanh(0.1*(t - 50)));
hem.set('model','controlVolume.source.energy' , @(~,~,~,~,t) ...
    [0; -source(t) ; zeros(6,1) ; source(t) ; zeros(5,1)]);
hem.set('model','momentumCell.source.momentum', ...
    @(~,~,momentum,TD,~)...
    [0;sign(momentum(2))*Af*(TD.P(2) - P0) ; zeros(12,1)] ...
);

%   State
hem.set('model','controlVolume.mass'          , mass      );
hem.set('model','controlVolume.energy'        , energy    );
hem.set('model','controlVolume.volume'        , volume    );
hem.set('model','controlVolume.maximumVolume' , volumeMax );
hem.set('model','momentumCell.momentum'       , momentum  );

%   Non-dimensionalizers
hem.set('model','dimensionalizer.mass'     , max(mass)      );
hem.set('model','dimensionalizer.energy'   , max(energy)    );
hem.set('model','dimensionalizer.volume'   , max(volume)    );
hem.set('model','dimensionalizer.momentum' , max(momentum)  );

%
hem.set('preconditioner','kind','block-jacobi');

%   Set evolution parameters
dt = 1e-4  ;
hem.set('evolver','time.span'         , [0,5e-2]  );
hem.set('evolver','time.step.maximum' , 1                               );
hem.set('evolver','time.step.minimum' , 1E-7                            );
hem.set('evolver','time.step.goal'    , dt                              );
hem.set('evolver','initialCondition'  , [mass;energy;volume;momentum]   );
hem.set('evolver','saveRate'          , dt                              );
tic;
hem.run();
[q{1},t{1}] = hem.results();
tElasped(1) = toc;


%% 
clear mass energy momentum volume P
n = 1;
mass{1,n}     = [];
energy{1,n}   = [];
momentum{1,n} = [];
volume{1,n}   = [];
P{1,n}        = [];
for k = 1:n
    mass{k}     = q{k}( (1:14) + 0   ,:)                    ;
    energy{k}   = q{k}( (1:14) + 14  ,:)                    ;
    volume{k}   = q{k}( (1:14) + 14*2,:)                    ;
    momentum{k} = q{k}( (1:14) + 14*3,:)                    ;
    rho         = bsxfun(@rdivide,mass{k},volume{k}(:,1))   ;
    i           = energy{k}./mass{k}                        ;
    T           = Temperature(rho,i)                        ;
    P{k}        = Pressure(rho,T)                           ;
end


%% 

%   Pressure plot
figure(1);
plot(...
    t{1}+1E-5,P{1}(1,:),...
    t{1}+1E-5,P{1}(2,:),'-o',...
    t{1}+1E-5,P{1}(3,:),'--+',...
    t{1}+1E-5,P{1}(4,:),'--s','LineWidth',2);
title('Location: Volume 1 (Top-Right Corner)');
xlabel('Simulation Time [s]');
ylabel('Pressure [Pa]');
legend('1','2','3','4');
grid('on');


%   Mass balance plot
figure(2);
% plot(...
%     t{1},sum(mass{1})/sum(mass{1}(:,1)) - 1,...
%     t{2},sum(mass{2})/sum(mass{2}(:,1)) - 1,'-o',...
%     t{3},sum(mass{3})/sum(mass{3}(:,1)) - 1,'-+',...
%     t{4},sum(mass{4})/sum(mass{4}(:,1)) - 1,'-s','LineWidth',2);
plot(...
    t{1}+1E-5,mass{1}(1,:),...
    t{1}+1E-5,mass{1}(2,:),'-o',...
    t{1}+1E-5,mass{1}(3,:),'--+',...
    t{1}+1E-5,mass{1}(4,:),'--s','LineWidth',2);
title('Location: Whole System');
xlabel('Simulation Time [s]');
ylabel('Relative Mass Defect [\Delta{kg}/kg]');
legend(strcat('\Delta{t} = ',num2str(dt.','%4.2E'),{'s;   '},'t_{run} = ',num2str(tElasped.','%5.2f'),'s'),'Location','SouthEast');
grid('on');

%   Energy balance plot
figure(3);
% plot(...
%     t{1},sum(energy{1})/sum(energy{1}(:,1))-1,...
%     t{2},sum(energy{2})/sum(energy{2}(:,1))-1,'-o',...
%     t{3},sum(energy{3})/sum(energy{3}(:,1))-1,'-+',...
%     t{4},sum(energy{4})/sum(energy{4}(:,1))-1,'-s','LineWidth',2);
plot(...
    t{1}+1E-5,energy{1}(1,:),...
    t{1}+1E-5,energy{1}(2,:),'-o',...
    t{1}+1E-5,energy{1}(3,:),'--+',...
    t{1}+1E-5,energy{1}(4,:),'--s','LineWidth',2);
title('Location: Whole System');
xlabel('Simulation Time [s]');
ylabel('Relative Energy Defect [\Delta{J}/J]');
legend(strcat('\Delta{t} = ',num2str(dt.','%4.2E'),{'s;   '},'t_{run} = ',num2str(tElasped.','%5.2f'),'s'),'Location','SouthEast');
grid('on');

% 
% % mass     = q( 1:12    ,2)   ;
% energy   = q((1:12)+12,2)   ;
% volume   = q((1:12)+24,2)   ;
% momentum = q((1:12)+36,2)   ;
% 
% hem.set('model','controlVolume.mass'          , mass      );
% hem.set('model','controlVolume.energy'        , energy    );
% hem.set('model','controlVolume.volume'        , volume    );
% hem.set('model','momentumCell.momentum'       , momentum  );
% 
% hem.set('model','dimensionalizer.mass'     , max(mass)      );
% hem.set('model','dimensionalizer.energy'   , max(energy)    );
% hem.set('model','dimensionalizer.volume'   , max(volume)    );
% hem.set('model','dimensionalizer.momentum' , max(momentum)  );
% 
% hem.set('evolver','time.span'         , [1E-3,5]                        );
% hem.set('evolver','time.step.maximum' , 1                               );
% hem.set('evolver','time.step.minimum' , 1E-7                            );
% hem.set('evolver','time.step.goal'    , 3E-5                            );
% hem.set('evolver','initialCondition'  , [mass;energy;volume;momentum]   );
% hem.set('evolver','saveRate'          , 0.1                            );
% 
% hem.run();
% [q,t] = hem.results();
% 


