function hem = testHEM(P0,T0,q0,mass,energy,~,momentum,T,P)
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
    vertMap = [...
        -1 ; -1 ; -1 ; -1 ; -1  ;   ... % Downcomer
        0 ;  0 ;                   ... % 1st lower run
        +1 ; +1 ; +1 ; +1 ;        ... % 1st riser
        0 ;  0 ;                   ... % 2nd lower run
        +1 ; +1 ; +1 ; +1 ; +1  ;   ... % 2st riser
        0 ;  0 ;  0 ;  +1 ;  0  ];  % Upper run
    horzMap = [...
        0 ;  0 ;  0 ;  0 ;  0  ;   ... % Downcomer
        +1 ; +1 ;                   ... % 1st lower run
        0 ;  0 ;  0 ;  0 ;   ... % 1st riser
        +1 ; +1 ;                   ... % 2nd lower run
        0 ;  0 ;  0 ;  0 ;  0  ;   ... % 2st riser
        -1 ; -1 ; -1 ; 0 ; -1           ];  % Upper run
    
    %   Path lengths / hydraulic diameters
    L  = [   hcorn/2 + hvert/2   ;
        hvert*ones(3,1)     ;
        hcorn/2 + hvert/2   ];
    L  = [L ; (hcorn + hhorz)*[1;1]/2 ; L ; (hcorn + hhorz)*[1;1]/2;  L ; (hcorn + hhorz)*[1;1;1;1]/2 ];
    
    nMC = numel(vertMap);
    
    Dh = 0.1 * ones(nMC,1);
    volume       = [...
        volCorn                                 ;
        volVert*ones(4,1)                       ;
        volCorn ; volHorz ; volCorn             ;
        volVert*ones(4,1)                       ;
        volHorz ; volCorn                       ;
        volVert*ones(4,1)                       ;
        volCorn ; volHorz ; volCorn ; volHorz   ];
    
    
    
    %   Instantiate the simulation
    hem = IntrepidTwilight.AdamantWave.HEM();
    
    %  From/To Volumes
    nCV  = 22;
    from = [1:11 , 8  , 13:21 , 12 , 22].';
    to   = [2:12 , 13 , 14:22 , 21 ,  1].';
    volumeFractionFrom = [ones(7,1)/2 ; 2/8 ; ones(3,1)/2 ; 3/8 ; ones(8,1)/2 ; 3/8 ; 1/2 ; 1/2];
    volumeFractionTo   = [ones(6,1)/2 ; 3/8 ; ones(4,1)/2 ; ones(8,1)/2 ; 3/8 ; 0.5 ; 2/8 ; 1/2 ];
    
    %   Upwind/downwind cell
    up   = [1:5 , 6:10 , 7  , 12 , 12:20 , 11 , 20 , 22 , 21 , 23].';
    down = [2:6 , 7:11 , 12 , 8  , 13:21 , 22 , 22 , 21 , 23 ,  1].';
    nx   =    [ zeros(4,1) ; +r2 ; +1 ; +r2   ; zeros(3,1) ; +1  ; -r2  ; +1 ; +r2 ; zeros(4,1) ; -r2 ; -1 ; -1  ; 0 ; -r2  ; -r2  ; -1 ; -r2];
    ny   =    [ -ones(4,1) ; -r2 ;  0 ; +r2   ;  ones(3,1) ;  0  ; +r2  ;  0 ; +r2 ;  ones(4,1) ; +r2 ;  0 ;  0  ; 1 ; -r2  ; +r2  ;  0 ; -r2];
    aI   = Af*[  ones(4,1) ;  s2 ;  1 ;  s2/2 ;  ones(3,1) ; 0.5 ; s2/2 ;  1 ;  s2 ;  ones(4,1) ;  s2 ;  1 ; 1/2 ; 1 ; s2/2 ; s2/2 ;  1 ;  s2];
    
    
    hem.set('model','momentumCell.from', from);
    hem.set('model','momentumCell.to'  , to);
    %
    %   momentum cell geometry
    hem.set('model','momentumCell.directionX'         , horzMap             );
    hem.set('model','momentumCell.directionY'         , vertMap             );
    hem.set('model','momentumCell.volumeFractionFrom' , volumeFractionFrom  );
    hem.set('model','momentumCell.volumeFractionTo'   , volumeFractionTo    );
    hem.set('model','momentumCell.flowArea'           , Af*ones(nMC,1)      );
    hem.set('model','momentumCell.LoD'                , L./Dh               );
    hem.set('model','momentumCell.source.friction'    , 0.1                 );
    %
    %   Interfaces
    hem.set('model','interface.up'  , up );
    hem.set('model','interface.down', down);
    hem.set('model','interface.normalX',  nx);
    hem.set('model','interface.normalY',  ny);
    hem.set('model','interface.flowArea', aI);
    
    
    
    
    
    %   Define initial state
    if (nargin < 4)
        %   Get densities and energies assuming incompressible, hydrostatic pressures
        v0    = 0.3                                                    ; % [m/s]
        rho0  = Density(P0,T0,0)                                        ;
        Pmax  = P0 + 9.81 * rho0 * 1.3                                  ;
        P     = [...
            linspace(P0,Pmax,6).';Pmax;linspace(Pmax,P0+0.2*9.81*rho0,5).';...
            Pmax;linspace(Pmax,P0,6).';P0;P0;P0];
        T     = P*0 + T0;
        rho   = Density(P,P*0 + T0,P*0)                                 ;
        i     = InternalEnergy(rho,T0)                                  ;
        
        %   Define extensive properties
        volumeMax    = volume           ;
        mass         = rho .* volume    ;
        energy       = i   .* mass      ;
        momentum     = (volumeFractionFrom .* mass(from) + volumeFractionTo .* mass(to)) .* v0;
        
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
    hem.set('semidiscretization','isDynamicVolume',[false;true(nCV-1,1)]);
    
    %   Set guard
    hem.set('residual','guard.step', @(q,dq) guardStep(q,dq));
    
    %   Preconditioning
    hem.set('preconditioner','kind','block-jacobi');
    
    
    %   Sources
    hem.set('model','controlVolume.source.mass'   , @(varargin) 0);
    shape  = @(t) min([(t>=tstart)* (1-0)/(thold-tstart)*(t-tstart),1]);
    source = @(t) q0 * shape(t);
    hem.set('model','controlVolume.source.energy' , @(~,~,~,~,t) ...
        [zeros(8,1); source(t)/2 ; zeros(5,1) ; source(t)/2 ; zeros(7,1)]);
    hem.set('model','momentumCell.source.momentum', @(varargin)0);
    
    
    %   Solver
    hem.set('solver','isNormNotDone',@(x,r) ...
        (norm(r(1:nCV*3),2)      > 1E-5 ) || ... % Control volume relative measure
        norm(r((nCV*3+1):end),1) > 1E-3 );
    hem.set('solver','maximumIterations' , 10);
    
    
    %   Set evolution parameters
    hem.set('evolver','initialCondition'  , [mass;energy;volume;momentum]   );
    hem.set('evolver','tolerance.relative',0.10);
    hem.set('evolver','tolerance.absolute',100 );
    
    
    
    function dq = guardStep(q,dq)
        while any( (q(1:nCV*3) - dq(1:nCV*3)) < 0 )
            dq = 0.5 * dq;
        end
    end
end






