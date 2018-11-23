function y = InterpolacionLagrange(X,pto)
    sz = size(X,1)

    for i = 1:sz
        nume = 1
        deno = 1
        for j = 1:sz
            if (i <> j)
                nume = (pto - X(j,1)) * nume 
            end
        end
        for j = 1:sz
            if (i <> j)
                deno = (X(i,1) - X(j,1)) * deno 
            end
        end
        l(i) = nume/deno;
    end
    y = 0;
    for i = 1:sz
        y = l(i) * X(i, 2) + y
    end
    
endfunction

//LAS CARIES NICO

function D = DDs(X) //asumimo que funca
    sz = size(X,1)
    for i = 1:sz
        D(i,1) = X(i,2)
    end
    
    for i = 2:sz
        for j = 1:sz-(i-1)
            D(j,i) = ( D(j+1,i-1) - D(j,i-1)) / ( X(j+i-1,1) - X(j,1) )
        end
    end
    disp(D)
endfunction


function y = InterpolacionNewton(X,pto)
    sz = size(X,1)
    y = X(1,2)
    D = DDs(X)
    for i = 2:sz
        pr = 1
        for j = 1:i-1
            pr = (pto - X(j,1)) * pr
        end 
        y = y + pr*D(1,i)
    end
endfunction
