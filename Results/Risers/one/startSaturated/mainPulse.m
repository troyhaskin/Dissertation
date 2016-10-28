clc();
clear();
T0 = 372;
P0 = SaturationStateGivenTemperature(T0);

%   372K steady-state/ 1kW
data = load('qSat_372K-100seconds.mat');

mass0     = data.mass{1}(:,end);
energy0   = data.energy{1}(:,end);
volume0   = data.volume{1}(:,end);
momentum0 = data.momentum{1}(:,end);
Tg        = Temperature(mass0./volume0,energy0./mass0);
Pg        = Pressure(mass0./volume0,Tg);

%   Adjust evolver time stuff
dt          = 0.02              ;
qPulse      = [1,2,5,10,15]*1E3;
m           = numel(qPulse)     ;
tElapsed(m) = 0                 ;
q{m}        = []          ;
Dq{m}       = []          ;
t{m}        = []          ;
mass{m}     = []          ;
energy{m}   = []          ;
volume{m}   = []          ;
momentum{m} = []          ;
rho{m}      = []          ;
i{m}        = []          ;
T{m}        = []          ;
state{m}    = []          ;
P{m}        = []          ;
newstate    = {mass0,energy0,volume0,momentum0};
q0          = 16000        ;

tall = tic;
for k = 1:m
    %
    %   Setup
    hem = testHEMPulse(P0,T0,q0,newstate{:},qPulse(k));
    hem.set('evolver','time.span'         , [0,150] );
    hem.set('evolver','time.step.maximum' ,  [0,0.50;95,0.02;105,0.25] );
    hem.set('evolver','time.step.minimum' ,   1E-8   );
    hem.set('evolver','time.step.goal'    ,   dt     );
    hem.set('evolver','saveRate'          ,  0.01    );
    %
    %   Run
    trun = tic;
    hem.run();
    tElapsed(k) = toc(trun);
    [q{k},Dq{k},t{k}] = hem.results()       ;
    %
    %   Post-process
    mass{k}         = q{k}( (1:14) + 0   ,:)    ;
    energy{k}       = q{k}( (1:14) + 14  ,:)    ;
    volume{k}       = q{k}( (1:14) + 14*2,:)    ;
    momentum{k}     = q{k}( (1:14) + 14*3,:)    ;
    rho{k}          = mass{k}   ./ volume{k}    ;
    i{k}            = energy{k} ./ mass{k}      ;
    [T{k},state{k}] = Temperature(rho{k},i{k})  ;
    P{k}            = Pressure(rho{k},T{k})     ;
    
%     newstate = {mass{k}(:,end),energy{k}(:,end),volume{k}(:,end),momentum{k}(:,end),T{k}(:,end),P{k}(:,end)};
end
toc(tall);
