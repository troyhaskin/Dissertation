function sol = fdmSolve(x,b,ul,ur,Tana)
    
    n  = numel(x)   ;
    dx = x(2) - x(1);
    
    %   FD Matrix
    A = spdiags(ones(n,1)*[1,-2,1]/dx^2,[-1,0,1],n,n);
    
    %   Left boundary
    dirchlet = ul(1);
    neumann  = ul(2);
    value    = ul(3);
    A(1,:)   = 0;
    if (neumann == 0)
        A(1,1)   = dirchlet ;
        b(1)     = value    ;
    else
        A(1,1)   = 2*(dirchlet/(dx*neumann) - 1./dx^2)  ;
        A(1,2)   = 2/dx^2                               ;
        b(1)     = b(1) + 2/(dx*neumann)*value          ;
    end
    
    %   Right boundary
    dirchlet = ur(1);
    neumann  = ur(2);
    value    = ur(3);
    A(n,:)   = 0;
    if (neumann == 0)
        A(n,n)   = dirchlet ;
        b(n)     = value    ;
    else
        A(n,n)   = -2*(dirchlet/(dx*neumann) + 1./dx^2) ;
        A(n,n-1) = 2/dx^2                               ;
        b(n)     = b(n) - 2/(dx*neumann)*value          ;
    end
    
    
    %   Initialization for closure
    u = A\b;


    %   Get Values
    sol.reconstruct = @(x) reconstruct(x);
    sol.error  = calculateError();

    
    function [ui,dui] = reconstruct(xi)
        
        if (nargin == 0) || isempty(xi)
            ui = u;
            
            if (nargout == 2)
                dui = [];
            end
            
            return;
        end


        ui    = xi;
        dui   = ui;
        
        for k = 2:n-1
            
            xkm1 = x(k-1)   ;
            xk   = x(k)     ;
            xkp1 = x(k+1)   ;
            
            mask = (xkm1 <= xi) & (xi <= xkp1);
            xv   = xi(mask);
            
            ui(mask) = (xv - xk  ).*(xv - xkp1)/(2*dx^2) * u(k-1) -...
                       (xv - xkm1).*(xv - xkp1)/(  dx^2) * u(k)   +...
                       (xv - xkm1).*(xv - xk  )/(2*dx^2) * u(k+1) ;
                   
            if (nargout == 2)
                dui(mask) = (2*xv - xk   - xkp1)/(2*dx^2) * u(k-1)  -...
                            (2*xv - xkm1 - xkp1)/(  dx^2) * u(k)    +...
                            (2*xv - xkm1 - xk  )/(2*dx^2) * u(k+1)  ;
            end
            
        end
    end
    
	function err = calculateError()
        err = 0;
        for k = 2:n-1
            
            xkm1 = x(k-1)   ;
            xk   = x(k)     ;
            xkp1 = x(k+1)   ;
            
            T = @(x)   (x - xk  ).*(x - xkp1)/(2*dx^2) * u(k-1) -...
                       (x - xkm1).*(x - xkp1)/(  dx^2) * u(k)   +...
                       (x - xkm1).*(x - xk  )/(2*dx^2) * u(k+1) ;

            integrand = @(x) abs(T(x) - Tana(x));
            err       = err + integral(integrand,x(k-1),x(k));
        end
        integrand = @(x) abs(T(x) - Tana(x));
        err       = err + integral(integrand,x(k),x(k+1));
        err       = (err);
    end
    
    
    
    
end