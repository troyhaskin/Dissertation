clc();
clear();
T0 = 372:0.2:390;
P0 = SaturationStateGivenTemperature(T0);

%   372K steady-state/ 1kW
a = [...
    +9.5917134654076430E-01	+2.8775184509942111E+00	+2.8775234831546337E+00	+2.8775285152819916E+00	+2.8775335473755033E+00	+9.5917896682425186E-01	+1.9183579274539335E+00	+9.5917896062858365E-01	+2.8770385737103275E+00	+2.8770335055450564E+00	+2.8770284378509818E+00	+2.8770233700419894E+00	+7.7716934633516255E-01	+1.2901988230822263E+00
    +3.9730252983823430E+05	+1.1919020160664497E+06	+1.1918956582313345E+06	+1.1918893003679758E+06	+1.1918829424768141E+06	+3.9729290196805558E+05	+7.9458581176287425E+05	+3.9729290979594988E+05	+1.1944773668821142E+06	+1.1944837689584994E+06	+1.1944901705916678E+06	+1.1944965725233327E+06	+3.2265189635973045E+05	+5.3561002290698711E+05
    +1.0000000000000002E-03	+3.0000000000000005E-03	+3.0000000000000005E-03	+3.0000000000000005E-03	+3.0000000000000005E-03	+1.0000000000000002E-03	+2.0000000000000005E-03	+1.0000000000000002E-03	+3.0000000000000005E-03	+3.0000000000000005E-03	+3.0000000000000005E-03	+3.0000000000000005E-03	+1.0000000000000002E-03	+2.0000000000000005E-03
    +2.4665494487187565E+00	+3.6998231541676998E+00	+3.6998231541415714E+00	+3.6998231541145343E+00	+2.4665473312634787E+00	+1.8499099555551224E+00	+1.8499099575456479E+00	+2.4662262605216250E+00	+3.6998166574140039E+00	+3.6998166547648172E+00	+3.6998166516669011E+00	+2.3496255013716696E+00	+1.6403298957363441E+00	+2.1501200189764029E+00];

mass0     = a(1,:).';
energy0   = a(2,:).';
volume0   = a(3,:).';
momentum0 = a(4,:).';
Tg        = Temperature(mass0./volume0,energy0./mass0);
Pg        = Pressure(mass0./volume0,Tg);

%   Adjust evolver time stuff
dt          = 50          ;
m           = numel(P0)   ;
tElapsed(m) = 0           ;
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
newstate    = {mass0,energy0,volume0,momentum0,Tg,Pg};
q0          = 16000        ;

tall = tic;
for k = 1:m
    %
    %   Setup
    hem = testHEMSat(P0(k),T0(k),q0,newstate{:});
    hem.set('evolver','time.span'         , [0,5000] );
    hem.set('evolver','time.terminator'   , ...
        constructTerminator(hem.get('evolver','initialCondition')));
    hem.set('evolver','time.step.maximum' ,    dt    );
    hem.set('evolver','time.step.minimum' ,   1E-12  );
    hem.set('evolver','time.step.goal'    ,   dt     );
    hem.set('evolver','saveRate'          ,  1E5 );
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
    
    newstate = {mass{k}(:,end),energy{k}(:,end),volume{k}(:,end),momentum{k}(:,end),T{k}(:,end),P{k}(:,end)};
end
toc(tall);
