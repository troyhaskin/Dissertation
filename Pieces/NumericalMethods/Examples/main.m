function [solfem,solfdm,solfvm,Tana,qana,Tmax,qmax] = main(ne,BC)
    
    nn = ne + 1             ;
    x  = linspace(-1,1,nn)' ;
    
    %   Analytical
    switch(BC)
        case('DirichletDirichlet')
            Tana = @(x) (1+1*x+4*cos(pi*x/2))/4;
            qana = @(x) -(1 - 2*pi*cos(pi*x/2))/8;
            
            %   Maximums
            xTmax = 2*asin(1/2/pi)/pi;
            Tmax  = Tana(xTmax);
            xqmax = 1;
            qmax  = qana(xqmax);
            
        case('DirichletNeumann')
            Tana = @(x) (4*pi - pi^2 + (4-pi)*pi*x + 8*cos(pi*x/2))/(2*pi^2);
            qana = @(x) -(4*pi - pi^2 - 4*pi*sin(pi*x/2))/(2*pi^2);
            
            %   Maximums
            xTmax = 2*asin(1-pi/4)/pi;
            Tmax  = Tana(xTmax);
            xqmax = 1;
            qmax  = qana(xqmax);
    end
    
    
    
    %   FEM
    solfem = femSolve(x,BC,Tana);
    
    %   FDM
    switch(BC)
        case('DirichletDirichlet')
            ul = [1,0,0]    ;
            ur = [1,0,2/4]  ;
        case('DirichletNeumann')
            ul = [1,0,0]    ;
            ur = [0,1,-0.5] ;
    end
    b      = -cos(pi*x/2)*pi^2/4        ;
    solfdm = fdmSolve(x,b,ul,ur,Tana)   ;
    
    
    switch(BC)
        case('DirichletDirichlet')
            ul = [1,0,0]    ;
            ur = [1,0,2/4]  ;
        case('DirichletNeumann')
            ul = [1,0,0]    ;
            ur = [0,1,-0.5] ;
    end
    Q      = @(x) -cos(pi*x/2)*pi^2/4   ;
    solfvm = fvmSolve(x,Q,ul,ur,Tana)   ;

end