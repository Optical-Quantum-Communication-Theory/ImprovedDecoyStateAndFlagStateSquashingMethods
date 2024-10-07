function unitaryTrans = pauliBasis(basis,revY)
unitaryTrans = eye(2);
switch basis
    case 1
        %Z basis, no transformation needed
    case 2
        %X basis, apply Hadamard
        unitaryTrans = [1,1;1,-1]/sqrt(2); %U|1> = |+>, U|2> = |->
    case 3
        % Y basis
        if revY
            unitaryTrans = [1,1;-1i,1i]/sqrt(2);% U|1> = |L>, U|2> = |R>
        else
            unitaryTrans = [1,1;1i,-1i]/sqrt(2);% U|1> = |R>, U|2> = |L>
        end
    otherwise
        error("Basis must be 1 (Z), 2 (X), or 3 (Y). Enterered: %d",basis);
end
end