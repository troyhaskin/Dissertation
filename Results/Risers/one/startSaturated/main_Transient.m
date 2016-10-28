clc();
% clear();
T0 = 372*ones(5,1);
P0 = SaturationStateGivenTemperature(T0);

%   372K steady-state/ 1kW
a = [...
    +9.5917134654076430E-01	+2.8775174744965946E+00	+2.8775225165173244E+00	+2.8775275585048377E+00	+2.8775326004586050E+00	+9.5917865372098077E-01	+1.9183573052158913E+00	+9.5917865149473458E-01	+2.8774621946172987E+00	+2.8774571380234497E+00	+2.8774520813597984E+00	+2.8774470245942414E+00	+8.9578460637339130E-01	+1.7475040557360018E+00
    +3.9730252983823430E+05	+1.1919032498110230E+06	+1.1918968795360909E+06	+1.1918905092326200E+06	+1.1918841389013105E+06	+3.9729329757621570E+05	+7.9458659796730848E+05	+3.9729330039076495E+05	+1.1922534399356653E+06	+1.1922598400572340E+06	+1.1922662404466625E+06	+1.1922726411497304E+06	+3.7116442544131016E+05	+7.2406592468353780E+05
    +1.0000000000000002E-03	+3.0000000000000005E-03	+3.0000000000000005E-03	+3.0000000000000005E-03	+3.0000000000000005E-03	+1.0000000000000002E-03	+2.0000000000000005E-03	+1.0000000000000002E-03	+3.0000000000000005E-03	+3.0000000000000005E-03	+3.0000000000000005E-03	+3.0000000000000005E-03	+1.0000000000000002E-03	+2.0000000000000005E-03
    +1.4444086829507281E+00	+2.1666129828997227E+00	+2.1666129828816034E+00	+2.1666129828630103E+00	+1.4444078115424694E+00	+1.0833055415169692E+00	+1.0833055419348971E+00	+1.4443796214382398E+00	+2.1666091786310040E+00	+2.1666091761827007E+00	+2.1666091736746220E+00	+1.4205518590045620E+00	+1.0655406006791779E+00	+1.1186020709106661E+00];




mass0     = a(1,:).';
energy0   = a(2,:).';
volume0   = a(3,:).';
momentum0 = a(4,:).';
Tg        = Temperature(mass0./volume0,energy0./mass0);
Pg        = Pressure(mass0./volume0,Tg);

%   Adjust evolver time stuff
dt          = 0.1         ;
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
q0          = [2:4:14,16]*1E3;

tall = tic;
for k = m:m
    %
    %   Setup
    hem = testHEMSat(P0(k),T0(k),q0(k),newstate{:});
    hem.set('evolver','time.span'         , [0,500] );
    hem.set('evolver','time.terminator'   , ...
        constructTerminator(hem.get('evolver','initialCondition')));
    hem.set('evolver','time.step.maximum' ,    dt    );
    hem.set('evolver','time.step.minimum' ,   1E-12  );
    hem.set('evolver','time.step.goal'    ,   dt     );
    hem.set('evolver','saveRate'          ,   0.05   );
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
