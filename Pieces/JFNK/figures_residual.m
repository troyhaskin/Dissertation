clc();
clear();

n     = 75  ;
aq     = 4   ;
x0    = 1   ;
rgoal = 0.75;
x     = linspace(-1,1.5,n);
%
r     = @(x)  exp(-(x+1/aq).^2) - rgoal      ;
Dr    = @(x)  -2*(x+1/aq).*exp(-(x+1/aq).^2)  ;
%
%
%   First step
m1    = Dr(x0)              ;
dx1   = r(x0)/m1            ;
x1    = x0 - dx1            ;
xnew1 = linspace(x0,x1)     ;
rlin1 = m1*(xnew1-x0)+r(x0) ;
%
%   Second step
m2    = Dr(x1)              ;
dx2   = r(x1)/m2            ;
x2    = x1 - dx2            ;
xnew2 = linspace(x1,x2)     ;
rlin2 = m2*(xnew2-x1)+r(x1) ;
%
%   Second step
m3    = Dr(x2)              ;
dx3   = r(x2)/m3            ;
x3    = x2 - dx3            ;
xnew3 = linspace(x2,x3)     ;
rlin3 = m3*(xnew3-x2)+r(x2) ;
%
limits = [-0.075,1.05,-0.6,0.3];
%
figure(1);
    plot([-10,10],[0,0],[0,0],[-10,10],'LineWidth',0.1,'Color',0.15*[1,1,1]);
    hold('on');
    lines = plot(x,r(x),...
        [x3,x3],[0,r(x3)],'--',x3,r(x3),'s',...
        xnew3,rlin3,[x2,x2],[0,r(x2)],'--',x2,r(x2),'s',...
        xnew2,rlin2,[x1,x1],[0,r(x1)],'--',x1,r(x1),'s',...
        xnew1,rlin1,[x0,x0],[0,r(x0)],'--',x0,r(x0),'s',...
        'LineWidth',2);
    colors      = get(gcf,'DefaultAxesColorOrder');
    lines(1).Color  = colors(1,:);
    lines(2).Color  = colors(2,:);
    lines(4).Color  = colors(3,:);
    lines(7).Color  = colors(4,:);
    lines(10).Color = colors(5,:);
    lines(3).Color  = lines(2).Color;
    lines(5).Color  = lines(4).Color;
    lines(6).Color  = lines(4).Color;
    lines(8).Color  = lines(7).Color;
    lines(9).Color  = lines(7).Color;
    lines(11).Color = lines(10).Color;
    lines(12).Color = lines(10).Color;
    %
    lines(3).MarkerFaceColor = lines(2).Color;
    lines(3).MarkerSize      = 9;
    lines(6).MarkerFaceColor = lines(4).Color;
    lines(6).MarkerSize      = 9;
    lines(9).MarkerFaceColor = lines(7).Color;
    lines(9).MarkerSize      = 9;
    lines(12).MarkerFaceColor = lines(10).Color;
    lines(12).MarkerSize      = 9;
    axis(limits);
    grid('on');
    xlabel('$x$','FontSize',16);
    ylabel('$r(x)$');
    legend(lines([1,12,9,6,3]),'Residual','Initial guess','First update','Second update','Third update','Location','SouthWest');
    f = gca;
    f.Children([10,7,4,1]) = f.Children([1,4,7,10]);
    matlab2tikz('.\JFNK_Example_GaussianNewtonSolution.tikz','width','5.5in','showInfo',false,'parseStrings',false);
    hold('off');


    
%%

clc();
clear();

%
%   Hand-coded stuff.  Abandoned because Hessian is too annoying to enter by hand, and
%   Mathematica's copy-paste sucks.
%
% % f     = @(x,y,mu) [exp(-mu.*(x+y).^2)./(2*mu);cosh(y)./(x.^2+1)];
% % r     = @(x,y,mu) f(x,y,mu) - f(zeros(size(x)),zeros(size(x)),mu);
% % Dr    = @(x,y,mu) [-(x+y).*exp(-mu.*(x+y).^2),-(x+y).*exp(-mu.*(x+y).^2)    ;...
% %                   -2*x.*cosh(y)./(x.^2+1).^2  , sinh(y)./(x.^2+1)            ];

% Calculate the determinant of the Hessian of r (and other stuff ... might as well).
syms x y mu f(x,y,mu) r(x,y,mu) rnorm(x,y,mu) 
f(x,y,mu)     = [exp(-mu*(x+y)^2)/(2*mu);cosh(y)/(x.^2+1)]  ;
r(x,y,mu)     = f(x,y,mu) - f(0,0,mu)                       ;
rnorm(x,y,mu) = r(x,y,mu)'*r(x,y,mu)                        ;
G             = gradient(rnorm(x,y,mu),[x,y])               ;
H             = hessian(rnorm(x,y,mu),[x,y])                ;
%
r     = matlabFunction(r(x,y,mu),'vars',[x,y,mu])                   ;
Dr    = matlabFunction(jacobian(r(x,y,mu),[x,y]),'vars',[x,y,mu])   ;
detHr = matlabFunction((H(1,1)*H(2,2)-H(1,2)*H(2,1))/(1+G(1)^2+G(2)^2)^2,'vars',[x,y,mu]);



%%

%   Setup
n     = 1000 ;  
x0    = 1;
y0    = 1;
mus   = linspace(0.19,1,4);
alpha = [linspace(0,1,n),linspace(1.001,1.5,n)];
%
%   Newton steps
dxs = arrayfun(@(mu) Dr(x0,y0,mu)\-r(x0,y0,mu),mus,'UniformOutput',false);
%
%   Residuals and Concavities
norm2 = @(v) sqrt(v(1,:).^2+v(2,:).^2);
rs    = cellfun(@(c) norm2(c),...
            arrayfun(@(k) ...
                r(x0+alpha*dxs{k}(1),y0+alpha*dxs{k}(2),mus(k)),...
                    1:numel(mus),'UniformOutput',false),'UniformOutput',false);
Hs = arrayfun( @(k) ...
            arrayfun( @(a) detHr(x0+a*dxs{k}(1),y0+a*dxs{k}(2),mus(k)),alpha),...
            1:numel(mus),'UniformOutput',false);
Hs = cellfun( @(c) c.*(abs(c)<10)+10*(abs(c)>=10) ,Hs,'UniformOutput',false);
dphi   = arrayfun( @(k) ...
            r(x0,y0,mus(k))'*Dr(x0,y0,mus(k))*dxs{k}/norm(r(x0,y0,mus(k))),1:numel(mus),'UniformOutput',false);
armijo = arrayfun( @(k) @(beta) ...
            norm(r(x0,y0,mus(k))) + ...
            beta * alpha*dphi{k},1:numel(mus),'UniformOutput',false);
% %
% %
% args  = [rs,rs];
% args(1:2:end) = {alpha};
% args(2:2:end) = rs;
% figure(1);
%     plot(args{:},'LineWidth',2);
%     axis([0,1,0,1.501]);
%     grid('on');
%     box('on');
%     legend(strcat('$\mu = ',num2str(mu'),'$'),'Location','NorthEast');
%     xlabel('$\alpha$');
%     ylabel('$\|r(\alpha)\|_2$','interpreter','none');
%     matlab2tikz('.\JFNK_Example_ResidualVsRelaxor.tikz','width','5.5in','showInfo',false,'parseStrings',false);
% %
% %
% args  = [Hs,Hs];
% args(1:2:end) = {alpha};
% args(2:2:end) = Hs;
% figure(2);
%     h = plot(args{:},'LineWidth',2);
%     axis([0,1,-2.16,1.65]);
%     grid('on');
%     box('on');
%     legend(strcat('$\mu = ',num2str(mu'),'$'),'Location','SouthEast');
%     xlabel('$\alpha$');
%     ylabel('$\det\left[H_r(\alpha)\right]$','interpreter','none');
%     matlab2tikz('.\JFNK_Example_ConvexityVsRelaxor.tikz','width','5.5in','showInfo',false,'parseStrings',false);
%
%
% 
% figure(3);
%     mask = @(x) x([1:10:n,n]);
%     h = plot(...
%         alpha,rs{1},                                alpha,rs{4},...
%         mask(alpha),mask(armijo{1}(1)),'--',mask(alpha),mask(armijo{4}(1)),'--',...
%         mask(alpha),mask(armijo{1}(0.5)),'-.',mask(alpha),mask(armijo{4}(0.5)),'-.',...
%         -1,0,'-.',-1,0,'--o','LineWidth',1.75);
%     h(3).Color = h(1).Color;
%     h(5).Color = h(1).Color;
%     h(2).Color = h(4).Color;
%     h(6).Color = h(4).Color;
%     %
%     h(3).Marker = 'o';
%     h(3).MarkerFaceColor = h(1).Color;
%     h(3).MarkerSize = 6;
%     %
%     h(4).Marker          = 'o';
%     h(4).MarkerFaceColor = h(4).Color;
%     h(4).MarkerSize      = 6;
%     %
%     scale= 0.40;
%     set(h([7,8]),'Color',scale*[1,1,1]);
%     h(8).MarkerFaceColor = [1,1,1]*scale;
%     %
%     h = legend(h([1,2,7,8]),{'$\mu = 0.19$','$\mu = 1$',...
%         'Armijo: $\beta=0.5$','Armijo: $\beta=1$'},'Interpreter','none');
%     axis([0,1,0,1.501]);
%     grid('on');
%     box('on');
%     xlabel('$\alpha$');
%     ylabel('$\|r(\alpha)\|_2$','interpreter','none');
%     matlab2tikz('.\JFNK_Example_ArmijoComparison.tikz','width','5.5in','showInfo',false,'parseStrings',false);
clc();

t      = 1      ;
mu     = mus(t) ;
dx     = dxs{t} ;
rv     = rs{t}  ;
%
alpha1 = 1;
r0     = rv(1);
r1     = rv(n);
%
aq = (r1 - rv(1) .* (1-alpha1))./alpha1.^2;
bq = rv(1);
Show(-bq/(2*aq));
%
dc     =  rv(1);
cc     = -rv(1);
alphab = 0.50  ;
rnormb = norm(r(x0+alphab*dx(1),y0+alphab*dx(2),mu),2);
A      = [alpha1^3,alpha1^2;alphab^3,alphab^2];
b      = [r1;rnormb] - cc*[alpha1;alphab] - dc;
ab     = A\b;
ac     = ab(1);
bc     = ab(2);
%
disp(' ');
Show([(-bc+sqrt(bc^2-3*ac*cc))/(3*ac),(-bc-sqrt(bc^2-3*ac*cc))/(3*ac)]);
%
epsilon = 1E-8;
alphah = 1;
r1v = r(x0+alphah*dx(1),y0+alphah*dx(2),mu);
r1  = norm(r1v,2);
m1  = (r(x0+(alphah+epsilon)*dx(1),y0+(alphah+epsilon)*dx(2),mu) - r1v)/epsilon;
m1  = r1v'*m1/r1;
hermiteC = @(t) (t-1).^2.*(t+1)*r0 + (-2*t.^3+3*t.^2)*r1 + (t.^3-t.^2)*alphah*m1;
%
disp(' ');
Show([(r0-3*r1+m1*alphah-sqrt((4*r0^2+(-3*r1+m1*alphah)^2+r0*(-12*r1+5*m1*alphah))))/(3*(r0-2*r1+m1*alphah)),...
      (r0-3*r1+m1*alphah+sqrt((4*r0^2+(-3*r1+m1*alphah)^2+r0*(-12*r1+5*m1*alphah))))/(3*(r0-2*r1+m1*alphah))]*alphah);
%
plot(alpha,rv,alpha,aq*alpha.^2 + bq*(1-alpha),alpha,ac*alpha.^3 + bc*alpha.^2 + cc*alpha + dc,alpha,hermiteC(alpha/alphah));
legend('True','Quadratic','Cubic','Hermite Cubic');
grid('on');
axis([0,1.5,0,2])










