function dq = guardStep(q,dq,qLo,qHi,problem)

    rho   = q(problem.miscellaneous.iRho)   .* problem.dimensionalizer.rho ;
    drho  = dq(problem.miscellaneous.iRho)  .* problem.dimensionalizer.rho ;
    
    rhoe  = q(problem.miscellaneous.iRhoe)  .* problem.dimensionalizer.rhoe ;
    drhoe = dq(problem.miscellaneous.iRhoe) .* problem.dimensionalizer.rhoe ;
    
    estar = 1 ./ rho /DimensioningInternalEnergy()  ;
    e     = rhoe  .* estar                          ;
    de    = drhoe .* estar                          ;
    
    rho   = rho  / CriticalDensity() ;
    drho  = drho / CriticalDensity() ;
    
    alpha = 1;
    while any(mayHaveIceDeltaIND(rho - alpha*drho,e - alpha*de))
        alpha = alpha * 0.5;
    end

    dq = alpha*dq;

end