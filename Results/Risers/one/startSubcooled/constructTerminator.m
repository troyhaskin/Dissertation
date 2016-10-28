function term = constructTerminator(q0)
    
    qs         = [repmat(q0*0,1,9),q0]          ;
    ts         = zeros(1,10)                    ;
    liveCounts = 1                              ;
    term       = @(t,q,Dq) terminator(t,q,Dq)   ;
    
    
    function stop = terminator(t,q,Dq)
        if (t > 2)
            qs         = [qs(:,2:10),q] ;
            ts         = [ts(2:10),t]   ;
            liveCounts = liveCounts + 1 ;
            valueWise  = ...
                (mean(max(abs(diff(qs,[],2)./qs(:,1:9)),[],1)) < 1E-3)  && ...
                (mean(diff(ts)) > 1E-1)    && (liveCounts >= 10)        ;
            slopeWise  = all(abs(Dq) < 10) && (liveCounts >= 1)         ;
            stop       = slopeWise                                      ;
        else
            stop = false;
        end
    end
    
end