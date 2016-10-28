clc();
clear();



% problem.miscellaneous.nCV      = max([problem.geometry.from;problem.geometry.to]);
% problem.miscellaneous.nMC      = length(problem.geometry.from);
% problem.miscellaneous.nInter   = length(problem.geometry.nx);
% problem.miscellaneous.nEq      = 2*problem.miscellaneous.nCV + problem.miscellaneous.nMC;

% problem.miscellaneous.sRho  = @(rho,rhoe,rhov,TD,t) 0;
% problem.miscellaneous.sRhov = @(rho,rhoe,rhov,TD,t) 0;

% energySource = (1:12)'*0;
% energySource(8) = +5E7;
% problem.miscellaneous.sRhoe = @(rho,rhoe,rhov,TD,t) 0;%((t > 0.01)&&(t<=0.02))*(energySource*(t-0.01)/(0.02-0.01)) + (t>0.02)*energySource;

%   Geometry

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
L = [   hcorn/2 + hvert/2   ;
        hvert*ones(3,1)     ;
        hcorn/2 + hvert/2   ;
        hcorn/2 + hhorz/2   ;
        hcorn/2 + hhorz/2   ;
        hcorn/2 + hvert/2   ;
        hvert*ones(3,1)     ;
        hcorn/2 + hvert/2   ;
        hcorn/2 + hhorz/2   ;
        hcorn/2 + hhorz/2   ];
Dh = 0.1 * ones(14,1);

%   Define initial state
P0   = 101325                   ; % [Pa]
T0   = 300                      ; % [K]
rho0 = Density(P0,T0)           ; % [kg/m^3]
i0   = InternalEnergy(rho0,T0)  ; % [J/kg]
v0   = 0                        ; % [m/s]

%   Get densities and energies assuming incompressible, hydrostatic pressures 
P   = P0 + 9.81 * [0; hcorn/2 + hvert*(1:2:7).'/2 ; hvert*4 + hcorn ];
P   = [P ; P(end) ; P(end:-1:1) ; P(1)];
T   = P*0 + T0;
rho = Density(P,T);
i   = InternalEnergy(rho,T);


%   Define extensive properties
volume       = [              volCorn ; 
                         volVert*ones(4,1) ; 
                    volCorn ; volHorz ; volCorn ; 
                         volVert*ones(4,1);
                         volHorz ; volCorn ];
volumeMax    = volume                                                   ;
mass         = rho .* volume                                            ;
energy       = i   .* mass                                              ;
momentum     = mass .* v0                                               ;
momentum(1)  = 1.1*momentum(1)                                          ;


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
hem.set('model','controlVolume.source.energy' , @(varargin) 0);
hem.set('model','momentumCell.source.momentum', @(varargin) 0);

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
hem.set('model','dimensionalizer.momentum' , 1  );

%   Set evolution parameters
hem.set('evolver','time.span'         , [0,0.1]                         );
hem.set('evolver','time.step.maximum' , 1                               );
hem.set('evolver','time.step.minimum' , 1E-7                            );
hem.set('evolver','time.step.goal'    , 1E-5                            );
hem.set('evolver','initialCondition'  , [mass;energy;volume;momentum]   );
hem.set('evolver','saveRate'          , 5E-4                            );

hem.run();
[q,~] = hem.results();





% mass     = q( 1:12    ,2)   ;
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


