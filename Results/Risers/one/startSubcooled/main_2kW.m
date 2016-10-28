clc();
% clear();

T0 = linspace(300,372,25);%373,373.1];
P0 = linspace(1,2,25)*101325;
[P0,T0] = meshgrid(P0,T0);

%   Adjust evolver time stuff
dt            = 20          ;
[m,n]         = size(P0)    ;
tElapsed(m,n) = 0           ;
q{m,n}        = []          ;
Dq{m,n}       = []          ;
t{m,n}        = []          ;
mass{m,n}     = []          ;
energy{m,n}   = []          ;
volume{m,n}   = []          ;
momentum{m,n} = []          ;
rho{m,n}      = []          ;
i{m,n}        = []          ;
T{m,n}        = []          ;
state{m,n}    = []          ;
P{m,n}        = []          ;
newstate      = {}          ;
q0            = 2E3         ;

tall = tic;
for k = 1:1
    %
    %   Setup
    hem = testHEM(P0(k),T0(k),q0,newstate{:});
    hem.set('evolver','time.span'         , [0,2000] );
    hem.set('evolver','time.terminator'   , ...
        constructTerminator(hem.get('evolver','initialCondition')));
    hem.set('evolver','time.step.maximum' ,    dt    );
    hem.set('evolver','time.step.minimum' ,   1E-12  );
    hem.set('evolver','time.step.goal'    ,   dt     );
    hem.set('evolver','saveRate'          ,   10   );
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


tall = tic;
for l = n:n
    for k = m:m
    newstate = {mass{k,l-1}(:,end),energy{k,l-1}(:,end),volume{k,l-1}(:,end),momentum{k,l-1}(:,end),T{k,l-1}(:,end),P{k,l-1}(:,end)};
    %
    %   Setup
    hem = testHEM(P0(k,l),T0(k,l),q0,newstate{:});
    hem.set('evolver','time.span'         , [0,2000] );
    hem.set('evolver','time.terminator'   , ...
        constructTerminator(hem.get('evolver','initialCondition')));
    hem.set('evolver','time.step.maximum' ,    dt    );
    hem.set('evolver','time.step.minimum' ,   1E-12  );
    hem.set('evolver','time.step.goal'    ,   dt     );
    hem.set('evolver','saveRate'          ,   10     );
    %
    %   Run
    trun = tic;
    hem.run();
    tElapsed(k,l) = toc(trun);
    [q{k,l},Dq{k,l},t{k,l}] = hem.results() ;
    %
    %   Post-process
    mass{k,l}         = q{k,l}( (1:14) + 0   ,:)        ;
    energy{k,l}       = q{k,l}( (1:14) + 14  ,:)        ;
    volume{k,l}       = q{k,l}( (1:14) + 14*2,:)        ;
    momentum{k,l}     = q{k,l}( (1:14) + 14*3,:)        ;
    rho{k,l}          = mass{k,l}   ./ volume{k,l}      ;
    i{k,l}            = energy{k,l} ./ mass{k,l}        ;
    [T{k,l},state{k,l}] = Temperature(rho{k,l},i{k,l})  ;
    P{k,l}            = Pressure(rho{k,l},T{k,l})       ;
    end
end


%     mass     = mass(1:m,1:n);
%     energy   = energy(1:m,1:n);
%     volume   = volume(1:m,1:n);
%     momentum = momentum(1:m,1:n);
%     rho      = rho(1:m,1:n);
%     i        = i(1:m,1:n);
%     T        = T(1:m,1:n);
%     state    = state(1:m,1:n);
%     P        = P(1:m,1:n);



%
% P0 = linspace(101325,2*101325,100);
% newstate = {mass{1}(:,end),energy{1}(:,end),volume{1}(:,end),momentum{1}(:,end)};
% for k = 1:n
%     for m = 2:numel(P0)
%         %
%         %   Setup
%         hem = testHEM(P0(m),T0(k),1000,newstate{:});
%         hem.set('evolver','time.span'         , [0,1000] );
%         hem.set('evolver','time.terminator'   ,@(t,q,Dq) (max(abs(Dq)) < 10) && (t > 10) );
%         hem.set('evolver','time.step.maximum' ,    dt    );
%         hem.set('evolver','time.step.minimum' ,   1E-12  );
%         hem.set('evolver','time.step.goal'    ,   dt     );
%         hem.set('evolver','saveRate'          ,   5      );
%         %
%         %   Run
%         trun = tic;
%         hem.run();
%         tElapsed(k,m) = toc(trun);
%         [q{k,m},Dq{k,m},t{k,m}] = hem.results()       ;
%         %
%         %   Post-process
%         mass{k,m}           = q{k,m}( (1:14) + 0   ,:)      ;
%         energy{k,m}         = q{k,m}( (1:14) + 14  ,:)      ;
%         volume{k,m}         = q{k,m}( (1:14) + 14*2,:)      ;
%         momentum{k,m}       = q{k,m}( (1:14) + 14*3,:)      ;
%         rho{k,m}            = mass{k,m}   ./ volume{k,m}    ;
%         i{k,m}              = energy{k,m} ./ mass{k,m}      ;
%         [T{k,m},state{k,m}] = Temperature(rho{k,m},i{k,m})  ;
%         P{k,m}              = Pressure(rho{k,m},T{k,m})     ;
%
%         newstate = {mass{k,m}(:,end),energy{k,m}(:,end),volume{k,m}(:,end),momentum{k,m}(:,end)};
%     end
% end
