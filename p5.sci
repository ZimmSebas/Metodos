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

function x = gausssolver(A,b,x0,eps) //NO TERMINADO
    
    sz = size(A,1)
    x = 1:sz
    xant = 1:sz
    suma = 0
    
    while(x )
        for i = 1:sz    
            for j = 1:sz
                suma = suma + A(i,j)*x0(j)  
            end
            
            mprintf("%f\n",suma)
        
            x(i) = 1/(A(i,i))*(b(i)-suma) 
        end
    
    
    end
    
endfunction
