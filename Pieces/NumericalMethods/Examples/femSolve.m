function solve = femSolve(x,BCtype,Tana)
    
    dx = x(2) - x(1);
    n  = numel(x)   ;
    
    switch(lower(BCtype))
        case('dirichletneumann')
            [A,b]        = dirichletNeumann(x)  ;
            u            = A\b                  ;
            T            = u(2:2:end)           ;
            q            = u(1:2:end)           ;
            solve.error  = calculateError()     ;
            solve.reconstruct = @(x)  reconstructDN(x)  ;
            
        case('dirichletdirichlet')
            [A,b]        = dirichletdirichlet(x)        ;
            T            = A\b                          ;
            q            = [-diff(T)./diff(x);0]        ;
            solve.reconstruct = @(x)  reconstructDD(x)  ;
            solve.error  = calculateError()     ;
    end
    
    
    
    %
    
    
    
    function [Ti,qi] = reconstructDD(xi)
        
        if (nargin == 0) || isempty(xi)
            Ti = T;
            
            if (nargout == 2)
                qi = q;
            end
            
            return;
        end
        
        
        Ti = xi;
        qi = Ti;
        
        for k = 2:n
            
            xkm1 = x(k-1)   ;
            xk   = x(k)     ;
            
            mask = (xkm1 <= xi) & (xi <= xk);
            xv   = xi(mask);
            
            Ti(mask) = (xv - xkm1)/dx * T(k) - (xv - xk)/dx * T(k-1);
            
            if (nargout == 2)
                qi(mask) = -(T(k) - T(k-1))/dx;
            end
            
        end
    end
    
    
    function [Ti,qi] = reconstructDN(xi)
        
        if (nargin == 0) || isempty(xi)
            Ti = T;
            
            if (nargout == 2)
                qi = q;
            end
            
            return;
        end
        
        
        Ti = xi;
        qi = Ti;
        
        for k = 2:n
            
            xkm1 = x(k-1)   ;
            xk   = x(k)     ;
            
            mask = (xkm1 <= xi) & (xi <= xk);
            xv   = xi(mask);
            
            Ti(mask) = (xv - xkm1)/dx * T(k) - (xv - xk)/dx * T(k-1);
            
            if (nargout == 2)
                qi(mask) = (xv - xkm1)/dx * q(k) - (xv - xk)/dx * q(k-1);
            end
            
        end
    end
    
    
    
    function err = calculateError()
        err = 0;
        for k = 2:n
            integrand = @(xx) abs(((x(k) - xx)*T(k-1)/dx + (xx - x(k-1))*T(k)/dx) - Tana(xx));
            err       = err + integral(integrand,x(k-1),x(k));
        end
        err = (err);
    end
    
    
    
    function [A,b] = dirichletNeumann(x)
        
        ne   = numel(x) - 1 ;
        xim1 = x(1:end-1)   ;
        xi   = x(2:end)     ;
        dxi  = diff(x)      ;
        
        %   Stiffness
        Kij = @(k) [...
            dxi(k)/3,-1/2,dxi(k)/6,1/2;
            -1/2,0,1/2,0;
            dxi(k)/6,-1/2,dxi(k)/3,1/2;
            -1/2,0,1/2,0];
        
        
        %   Load
        b1 = -(4*cos(pi*xi/2) - 4*cos(pi*xim1/2) + 2*pi*dxi.*sin(pi*xim1/2))./(dxi*4);
        b3 =  (4*cos(pi*xi/2) - 4*cos(pi*xim1/2) + 2*pi*dxi.*sin(pi*xi/2))./(dxi*4);
        
        %   Assemble stiffness
        I = repmat(bsxfun(@plus,(1:4)',2*(0:ne-1)),4,1);
        I = I(:);
        J = repmat(1:4,4,1);
        J = bsxfun(@plus,J(:),2*(0:ne-1));
        J = J(:);
        S = zeros(4*4*ne,1);
        bottom = 1;
        top    = 16;
        for e = 1:ne
            Ke            = Kij(e) ;
            S(bottom:top) = Ke(:)  ;
            %
            bottom = top + 1;
            top    = bottom + 16 - 1;
        end
        A = sparse(I,J,S, (ne+1)*2,(ne+1)*2);
        
        %   Assemble load
        b((ne+1)*2,1) = 0;
        b(2:2:end-2) = b1;
        b(4:2:end) = b(4:2:end) + b3;
        
        %   B.C.s
        A(2,:) = 0;
        A(2,2) = 1;
        b(2)   = 0;
        A(end-1,:)     = 0;
        A(end-1,end-1) = 1;
        b(end-1)       = 1/2;
    end
    
    function [A,b] = dirichletdirichlet(x)
        
        ne   = numel(x) - 1 ;
        xim1 = x(1:end-1)   ;
        xi   = x(2:end)     ;
        dxi  = diff(x)      ;
        
        %   Stiffness
        Kij = @(k) [...
            +1/dxi(k),-1/dxi(k);...
            -1/dxi(k),+1/dxi(k)];
        nk = 2;
        
        %   Load
        b1 = -2*(2*cos(xi*pi/2)-2*cos(pi*xim1/2)+pi*dxi.*sin(pi*xim1/2))./(dxi*4);
        b2 =  2*(2*cos(xi*pi/2)-2*cos(pi*xim1/2)+pi*dxi.*sin(pi*xi/2))  ./(dxi*4);
        
        %   Assemble stiffness
        I = repmat(bsxfun(@plus,(1:nk)',(0:ne-1)),nk,1);
        I = I(:);
        J = repmat(1:nk,nk,1);
        J = bsxfun(@plus,J(:),(0:ne-1));
        J = J(:);
        S = zeros(nk^2*ne,1);
        bottom = 1;
        top    = nk^2;
        for e = 1:ne
            Ke            = Kij(e) ;
            S(bottom:top) = Ke(:)  ;
            %
            bottom = top + 1;
            top    = bottom + nk^2 - 1;
        end
        A = sparse(I,J,S, (ne+1),(ne+1));
        
        %   Assemble load
%         b = dx*cos(pi/2*x)*pi^2/4;
%         b([1,end]) = b([1,end])/2;
        b(ne+1,1) = 0;
        b(1:end-1) = b1;
        b(2:end) = b(2:end) + b2;
        
        %   B.C.s
        A(1,:)     = 0  ;
        A(1,1)     = 1  ;
        b(1)       = 0  ;
        A(end,:)   = 0  ;
        A(end,end) = 1  ;
        b(end)     = 2/4;
    end
    
end