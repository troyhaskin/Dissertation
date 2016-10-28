function hem = testHEMPulse(P0,T0,q0,mass,energy,~,momentum,qPulse)
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
        v0    = 0.01                                                    ; % [m/s]
        rho0  = Density(P0,T0,0)                                        ;
        P     = flipud(P0 + 9.81 * rho0 * cumsum(vertMap.*L,'reverse')) ;
        rho   = Density(P,P*0 + T0,P*0)                                 ;
        i     = InternalEnergy(rho,T0)                                  ;

        %   Define extensive properties
        volumeMax    = volume           ;
        mass         = rho .* volume    ;
        energy       = i   .* mass      ;
        momentum     = mass .* v0       ;
        
        tstart = -1;
        thold  = 0;
        
    else
        
%         deltaP     = P0 - P(1)              ;
%         P          = P + deltaP             ;
%         deltaT     = T0 - T(1)              ;
%         T          = (T + deltaT)           ;
%         [~,state]  = Temperature(mass./volume,energy./mass);
%         rho        = Density(P,T,state.x)   ;
%         i          = InternalEnergy(rho,T)  ;
%         
        
        %   Define extensive properties
        volumeMax    = volume           ;
%         mass         = rho .* volume    ;
%         energy       = i   .* mass      ;
        
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
    source = @(t) q0 * shape(t) * ((t>tstart));
    hem.set('model','controlVolume.source.energy' , @(~,~,~,~,t) ...
        [-source(t); 0 ; zeros(6,1) ; source(t) + pulse() ; zeros(5,1)]);
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
        (norm(r(1:14*3),2)      > 1E-8 ) || ... % Control volume relative measure
        norm(r((14*3+1):end),1) > 1E-3 );
    hem.set('solver','maximumIterations' , 10);
   
    
    %   Set evolution parameters
    hem.set('evolver','initialCondition'  , [mass;energy;volume;momentum]   );
    hem.set('evolver','tolerance.relative',0.10);
    hem.set('evolver','tolerance.absolute',100 );
    
    
    function dq = guardStep(q,dq)
        dq(1:14*3) = dq(1:14*3).*[[false;true(13,1)];[false;true(13,1)];[false;true(13,1)]];
        while any( (q - dq) < 0 )
            dq = 0.5 * dq;
        end
    end
    
    
    hem.set('state','hook.preupdate',@(~,t,dt) activatePulse(t+dt));
    isPulsing    = false;
    hasNotPulsed = true;
    function [] = activatePulse(t)
        if (t>101) && hasNotPulsed
            isPulsing = true;
        end
    end
    hem.set('state','hook.postupdate',@(~,t,dt,stats) deactivatePulse(stats));
    function [] = deactivatePulse(stats)
        if isPulsing
            isPulsing = false;
            if strcmp(stats.returnStatus,'NormConverged')
                hasNotPulsed = false;
            end
        end
    end
    function q = pulse()
        if isPulsing
            q = qPulse;
        else
            q = 0;
        end
    end
    
end




