function x = solve2by2(A,b)
    b = b./(A(1,1)*A(2,2) - A(1,2)*A(2,1));
    x = [A(2,2)*b(1)-A(1,2)*b(2) ; A(1,1)*b(2)-A(2,1)*b(1)];
end