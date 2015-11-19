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
syms x y mu f(x,y,mu) r(x,y,mu) rnorm(x,y,mu) Hr
f(x,y,mu)     = [exp(-mu*(x+y)^2)/(2*mu);cosh(y)/(x.^2+1)]  ;
r(x,y,mu)     = f(x,y,mu) - f(0,0,mu)                       ;
rnorm(x,y,mu) = r(x,y,mu)'*r(x,y,mu)                        ;
Hr            = hessian(rnorm(x,y,mu),[x,y])                ;
%
r     = matlabFunction(r(x,y,mu),'vars',[x,y,mu])                       ;
Dr    = matlabFunction(jacobian(r(x,y,mu),[x,y]),'vars',[x,y,mu])       ;
detHr = matlabFunction(Hr(1,1)*Hr(2,2)-Hr(1,2)*Hr(2,1),'vars',[x,y,mu]) ;


%%
%   Setup
n     = 200 ;  
x0    = 1;
y0    = 1;
mu    = linspace(0.19,1,4);
alpha = [linspace(0,0.25,n/2),linspace(0.251,1.1,n/2)];
%
%   Newton steps
dx    = arrayfun(@(mu) Dr(x0,y0,mu)\r(x0,y0,mu),mu,'UniformOutput',false);
%
%   Residuals and Concavities
norm2 = @(v) v(1,:).^2+v(2,:).^2;
rs    = cellfun(@(c) norm2(c),...
            arrayfun(@(k) ...
                r(x0-alpha*dx{k}(1),y0-alpha*dx{k}(2),mu(k)),...
                    1:numel(mu),'UniformOutput',false),'UniformOutput',false);
Hs = arrayfun( @(k) ...
            detHr(x0-alpha*dx{k}(1),y0-alpha*dx{k}(2),mu(k)),...
            1:numel(mu),'UniformOutput',false);
Hs = cellfun( @(c) c.*(abs(c)<10)+10*(abs(c)>=10) ,Hs,'UniformOutput',false);
%
%
%   Form quadratic fit
dxHat  = cellfun(@(c) c/norm(c),dx,'UniformOutput',false);
phi0s  = arrayfun(@(m) norm(r(x0,y0,m),2).^2,mu,'UniformOutput',false);
phi0ps = arrayfun(@(k) -2*r(x0,y0,mu(k)).'*(Dr(x0,y0,mu(k))*dxHat{k}),1:numel(mu),'UniformOutput',false);
phibs  = arrayfun(@(k) @(beta) norm2(r(x0-beta*dx{k}(1),y0-beta*dx{k}(2),mu(k))),1:numel(mu),'UniformOutput',false);
a      = cellfun( @(phi0,phi0p,phib) @(beta) (phib(beta) - (phi0 + phi0p*beta))/beta^2 ,phi0s,phi0ps,phibs,'UniformOutput',false);
p      = cellfun( @(a,phi0,phi0p) @(beta) a(beta)*alpha.^2 + phi0p*alpha + phi0 ,a,phi0s,phi0ps,'UniformOutput',false);
alphao = cellfun( @(phi0,phi0p,phib) @(beta) phi0p*beta^2/(2*(phi0+phi0p*beta-phib(beta))) ,phi0s,phi0ps,phibs,'UniformOutput',false);
%
%
%
args  = [rs,rs];
args(1:2:end) = {alpha};
args(2:2:end) = rs;
figure(1);
    beta = 1;
    n = 4;
    plot(args{:},alpha,p{n}(beta),'--','LineWidth',2);
    Show(alphao{n}(beta));
    axis([0,1.1,-2,3]);
    grid('on');
    box('on');
    legend(strcat('$\mu = ',num2str(mu'),'$'),'Location','NorthEast');
    xlabel('$\alpha$');
    ylabel('$\|r(\alpha)\|_2$','interpreter','none');


