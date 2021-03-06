clear();
clc();

%   Mesh
nes = round(logspace(1,4,10))';
% nes = 6;
%
BC = 'DirichletDirichlet';
% BC = 'DirichletNeumann'  ;

nne = numel(nes);
errorfem(nne,1) = 0;
errorfdm(nne,1) = 0;
errorfvm(nne,1) = 0;

for k = 1:numel(nes)
    [solfem,solfdm,solfvm,Tana,qana,Tmax,qmax] = main(nes(k),BC);
    errorfem(k,1) = solfem.error;
    errorfdm(k,1) = solfdm.error;
    errorfvm(k,1) = solfvm.error;
    disp(k);
end


%   Reconstruct solutions
xre          = linspace(-1,1,100);
[Tfem,qfem]  = solfem.reconstruct(xre)  ;
Tfdm  = solfdm.reconstruct(xre)  ;
% qfdm         = -DTfdm                   ;

xfv  = linspace(-1,1,nes(end)+1);
xbar = (xfv(1:end-1) + xfv(2:end))/2;
Tfvm = solfvm.values;


xana = linspace(-1,1,100);

figure(1);
switch(BC)
    case('DirichletDirichlet')
        
        clf();
        h = plot(xana,Tana(xana),'k--',xre,Tfem,'-',xre,Tfdm,'-',-10,-10,'-');
        h(4).MarkerFaceColor = h(4).Color;
        h(4).LineWidth = 1.5;
        hold('on');
        for k = 1:nes(end)
            dx = diff(xfv(k:k+1));
            g                    = plot(xfv(k:k+1) + dx*[0.02,-0.02],Tfvm([k,k]),'-');%,mean(xfv(k:k+1)),Tfvm(k),'o','LineWidth',1.5);
            g(1).Color     = h(4).Color;
            g(1).LineWidth = 1.75;
%             g(2).Color           = h(4).Color;
%             g(2).MarkerFaceColor = h(4).Color;
        end
        hold('off');
        h(1).LineWidth = 2.5;
        h(2).LineWidth = 1.5;
        h(3).LineWidth = 1.5;
        h(2).MarkerFaceColor = h(2).Color;
        h(3).MarkerFaceColor = h(3).Color;
        grid('on');
        xlabel('$x$ [-]');
        ylabel('$T$ [$\Delta$K]');
        axis([-1,1,0,1.30]);
        legend(h,'Analytical','Finite Element','Finite Difference','Finite Volume','Location','SouthEast');
%         matlab2tikz('..\Graphics\Discrete_Example_Conduction.tikz','width','4in','showInfo',false,'parseStrings',false);
        
    case('DirichletNeumann')
        
        clf();
        subplot(2,1,1);
        h = plot(xana,Tana(xana),xre,Tfem,'o--',xre,Tfdm,'o--',NaN,NaN,'o--');
        hold('on');
        g = plot(xbar,Tfvm,'o',xgam,Tgam,'--');
        hold('off');
        g(1).Color              = h(4).Color;
        g(2).Color              = h(4).Color;
        g(1).MarkerFaceColor    = h(4).Color;
        
        subplot(2,1,2);
        g = plot(xana,qana(xana),xre,qfem,'o--',xre,qfdm,'o--');
        
        h(1).LineWidth = 2;
        g(1).LineWidth = 2;
        h(2).MarkerFaceColor = h(2).Color;
        h(3).MarkerFaceColor = h(3).Color;
        g(2).MarkerFaceColor = g(2).Color;
        g(3).MarkerFaceColor = g(3).Color;
end
%%

figure(2);
h = loglog(1E6,1E6,'k-o',nes,errorfem,'-o',nes,errorfdm,'-o',nes,errorfvm,'-o','LineWidth',1.75);
h(2).MarkerFaceColor = h(2).Color;
h(3).MarkerFaceColor = h(3).Color;
h(4).MarkerFaceColor = h(4).Color;
h(2).MarkerSize = 6.1;
h(2).LineWidth = 2;
h(1).Color = h(2).Color;
h(1).MarkerFaceColor = h(1).Color;
grid('on');
axis('square');
legend(h([1,3,4]),'Finite Element','Finite Difference','Finite Volume','Location','NorthEast');
xlabel('Number of Elements');
ylabel('L\subs[\;]{2}');
%
axis([1E1,1E4,1E-9,1]);
%
 matlab2tikz('..\Graphics\Discrete_Example_ConductionError.tikz','width','4.5in','showInfo',false,'parseStrings',false);


polyfit(log(nes),log(errorfdm),1)
polyfit(log(nes),log(errorfem),1)
polyfit(log(nes),log(errorfvm),1)