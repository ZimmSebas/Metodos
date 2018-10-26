function radspec = jacobitest(A)
    sz = size(A,1)
    N = eye(sz,sz)
    I = N
    di = diag(A)
    for i = 1:sz
        N(i,i) = N(i,i) * di(i)
    end
    
    Nor = (I-inv(N)*A)
    
    radspec = max(abs(spec(Nor)))
    
    if(radspec < 1)
        mprintf("Converge \n")
    else
        mprintf("No converge\n")
    end
endfunction

function radspec = gausstest(A)
    sz = size(A,1)
    I = eye(sz,sz)
    N = zeros(sz,sz)
     
    for i = 1:sz
        for j = 1:sz
            if(i >= j)
                N(i,j) = A(i,j)
            end
        end
    end
    
    disp(I)
    disp(N)
    
    Nor = (I-inv(N)*A)
    
    radspec = max(abs(spec(Nor)))
    
    if(radspec < 1)
        mprintf("Converge \n")
    else
        mprintf("No converge\n")
    end
endfunction


function x = gausssolver(A,b,x0,eps)
    
    sz = size(A,1)
    x = x0
    xant = x
    suma = 0
    it = 1
    cont = 0
    
    while(abs(norm(x-xant)) > eps | cont == 0) 
    xant = x
        for i = 1:sz
            suma = 0
            for j = 1:i-1 
                suma = suma + A(i,j)*x(j)
            end
            
            for j = i+1:sz
                suma = suma + A(i,j)*x(j)
            end
            x(i) = 1/(A(i,i))*(b(i)-suma)
        end
     cont = cont + 1
    end
    
    mprintf("Cantidad de iteraciones: %d\n",cont);
    
endfunction

function x = jacobisolver(A,b,x0,eps) 
    
    sz = size(A,1)
    x = x0
    xant = x
    suma = 0
    it = 1
    cont = 0
    
    while(abs(norm(x-xant)) > eps | cont == 0) 
    xant = x
        for i = 1:sz
            suma = 0
            for j = 1:sz 
                if (i <> j)
                    suma = suma + A(i,j)*xant(j)
                end
            end
            x(i) = 1/(A(i,i))*(b(i)-suma)
        end
     cont = cont + 1
    end
    
    mprintf("Cantidad de iteraciones: %d\n",cont);
    
endfunction

