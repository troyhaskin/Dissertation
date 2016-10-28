function hem = experimentLikeHEM(P0,T0,q0,mass,energy,~,momentum,T,P)
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
    volume       = [...
        volCorn                     ;
        volVert*ones(4,1)           ;
        volCorn ; volHorz ; volCorn ;
        volVert*ones(4,1)           ;
        volCorn ; volHorz           ];
    
    %   Define initial state
   
    if (nargin < 4)
        %   Get densities and energies assuming incompressible, hydrostatic pressures
        v0    = 0.1                                                     ; % [m/s]
        rho0  = Density(P0,T0,0)                                        ;
        P     = flipud(P0 + 9.81 * rho0 * cumsum(vertMap.*L,'reverse')) ;
        rho   = Density(P,P*0 + T0,P*0)                                 ;
        i     = InternalEnergy(rho,T0)                                  ;

        %   Define extensive properties
        volumeMax    = volume           ;
        mass         = rho .* volume    ;
        energy       = i   .* mass      ;
        momentum     = mass .* v0       ;
        
        tstart = 0;
        thold  = 5;
    else
        
        if (T0 < 373)
            deltaT = T0 - T(1)  ;
            T      = T + deltaT ;
            deltaP = P0 - P(1)  ;
            P      = P + deltaP ;
        else
            T      = T*0+T0     ;
            deltaP = P0 - P(1)  ;
            P      = P + deltaP ;
        end 
        rho    = Density(P,T)           ;
        i      = InternalEnergy(rho,T)  ;
        
        
        %   Define extensive properties
        volumeMax    = volume           ;
        mass         = rho .* volume    ;
        energy       = i   .* mass      ;
        
        tstart = -1 ;
        thold  = 0  ;
    end
    
    
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
    shape  = @(t) min([(t>=tstart)* (1-0)/(thold-tstart)*(t-tstart),1]);
    source = @(t) q0 * shape(t);
    hem.set('model','controlVolume.source.energy' , @(~,~,~,~,t) ...
        [0; 0 ; zeros(6,1) ; source(t) ; zeros(5,1)]);
    hem.set('model','momentumCell.source.momentum', @(varargin)0);
    
    %   State
    hem.set('model','controlVolume.mass'          , mass     );
    hem.set('model','controlVolume.energy'        , energy   );
    hem.set('model','controlVolume.volume'        , volume    );
    hem.set('model','controlVolume.maximumVolume' , volumeMax );
    hem.set('model','momentumCell.momentum'       , momentum );
    
    %   Non-dimensionalizers
    hem.set('model','dimensionalizer.mass'     , mass(:)          );
    hem.set('model','dimensionalizer.energy'   , energy(:)        );
    hem.set('model','dimensionalizer.volume'   , volume(:)        );
    hem.set('model','dimensionalizer.momentum' , 1          );
    
    %
    hem.set('semidiscretization','isDynamicVolume',[false;true(13,1)]);
    
    %   Set guard
    hem.set('residual','guard.step', @(q,dq) guardStep(q,dq));
    
    %   Preconditioning
    hem.set('preconditioner','kind','block-jacobi');
    
    %   Solver
%         [...
%             x(1:14*3) .* (norm(r(1:14*3),2)>1E-6) ;...
%             abs(r((14*3+1):end)) > abs(1E-4*x((14*3+1):end))]
    hem.set('solver','isNormNotDone',@(x,r) ...
        (norm(r(1:14*3),2) > 1E-6) || ... % Control volume relative measure
        any(abs(r((14*3+1):end)) > abs(1E-4*x((14*3+1):end))));
    hem.set('solver','maximumIterations' , 10);
   
    
    %   Set evolution parameters
    hem.set('evolver','initialCondition'  , [mass;energy;volume;momentum]   );
    hem.set('evolver','tolerance.relative',0.10);
    hem.set('evolver','tolerance.absolute',100 );
end

function dq = guardStep(q,dq)
    while any( abs(dq(1:14*3)) > abs(0.001*q(1:14*3)) )
        dq = 0.5 * dq;
    end
    while any( (q(1:14*3) - dq(1:14*3)) < 0 )
        dq = 0.5 * dq;
    end
end




