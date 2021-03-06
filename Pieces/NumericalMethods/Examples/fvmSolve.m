function sol = fvmSolve(x,Q,ul,ur,Tana)
    
    n   = numel(x)      ;
    ncv = n - 1         ;
    dx  = x(2) - x(1)   ;
    
    %   FD Matrix
    A = spdiags(ones(n,1)*[1,-2,1]/dx^2,[-1,0,1],ncv,ncv);
    
    %   Load
    b = Q((x(1:end-1) + x(2:end))/2);
    
    %   Left boundary
%     dirchlet = ul(1);
    neumann  = ul(2);
    value    = ul(3);
    A(1,:)   = 0;
    if (neumann == 0)
        A(1,1)   = -3/dx^2;
        A(1,2)   = 1/dx^2;
        b(1)     = b(1) - 2*value/dx^2;
%         A(1,1)   = -4/dx^2;
%         A(1,2)   = 4/(3*dx^2);
%         b(1)     = b(1) - 8/(3*dx^2)*value;
    else
        A(1,1)   = -1/dx^2              ;
        A(1,2)   = +1/dx^2              ;
        b(1)     = b(1) + value/dx     ;
    end
    
    %   Right boundary
%     dirchlet = ur(1);
    neumann  = ur(2);
    value    = ur(3);
    A(ncv,:)   = 0;
    if (neumann == 0)
        A(ncv,ncv)   = -3/dx^2;
        A(ncv,ncv-1) = 1/dx^2;
        b(ncv)       = b(ncv) - 2*value/dx^2;
%         A(ncv,ncv)   = -4/dx^2;
%         A(ncv,ncv-1) = 4/(3*dx^2);
%         b(ncv)       = b(ncv) - 8/(3*dx^2)*value;
    else
        A(ncv,ncv)   = +1/dx^2              ;
        A(ncv,ncv-1) = -1/dx^2              ;
        b(ncv)       = value/dx - b(ncv)  ;
    end
    
    
    %   Initialization for closure
    T = A\b;


    %   Get Values
    sol.values = T                  ;
    sol.error  = calculateError()   ;

    
	function err = calculateError()
        err = 0;
        for k = 1:ncv
            integrand = @(x) abs(T(k) - Tana(x));
            err       = err + integral(integrand,x(k),x(k+1));
        end
        err = (err);
    end
    
    
    
    
end