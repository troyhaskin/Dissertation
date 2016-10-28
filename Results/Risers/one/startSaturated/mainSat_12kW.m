clc();
clear();
T0 = 372:0.2:390;
P0 = SaturationStateGivenTemperature(T0);

%   372K steady-state/ 1kW
a = [...
    +9.5917134654076430E-01	+2.8775177434587369E+00	+2.8775227807360935E+00	+2.8775278179805768E+00	+2.8775328551911383E+00	+9.5917873741172532E-01	+1.9183574706870063E+00	+9.5917873327480740E-01	+2.8773294903615367E+00	+2.8773244293572491E+00	+2.8773193684227301E+00	+2.8773143073949266E+00	+8.3989568010061444E-01	+1.5442872040460802E+00
    +3.9730252983823430E+05	+1.1919029099919682E+06	+1.1918965457039520E+06	+1.1918901813870338E+06	+1.1918838170421873E+06	+3.9729319183034258E+05	+7.9458638888792647E+05	+3.9729319705805468E+05	+1.1929503186420512E+06	+1.1929567170820334E+06	+1.1929631155505897E+06	+1.1929695142576748E+06	+3.4822046907912620E+05	+6.4024545964588365E+05
    +1.0000000000000002E-03	+3.0000000000000005E-03	+3.0000000000000005E-03	+3.0000000000000005E-03	+3.0000000000000005E-03	+1.0000000000000002E-03	+2.0000000000000005E-03	+1.0000000000000002E-03	+3.0000000000000005E-03	+3.0000000000000005E-03	+3.0000000000000005E-03	+3.0000000000000005E-03	+1.0000000000000002E-03	+2.0000000000000005E-03
    +2.0004044805334433E+00	+3.0006064504678576E+00	+3.0006064504449075E+00	+3.0006064504212033E+00	+2.0004031328305105E+00	+1.5003019098879700E+00	+1.5003019109654343E+00	+2.0002947743434678E+00	+3.0006011830320611E+00	+3.0006011808772586E+00	+3.0006011785879250E+00	+1.9382440148894284E+00	+1.4196171429009632E+00	+1.6214312153002022E+00];

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
q0          = 12000        ;

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
