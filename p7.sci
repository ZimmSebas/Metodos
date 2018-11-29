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

function D = DDs(X)
    sz = size(X,1)
    for i = 1:sz
        D(i,1) = X(i,2)
    end
    
    for i = 2:sz
        for j = 1:sz-(i-1)
            D(j,i) = ( D(j+1,i-1) - D(j,i-1)) / ( X(j+i-1,1) - X(j,1) )
        end
    end
    //disp(D)
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


// Ejercicio 7

function [Q,R] = FactoRQ(A)   // A debe tener columnas LI
    sz = size(A,2)
    Q(:,1) = A(:,1)/norm(A(:,1))
    V(1) = norm(A(:,1))
    for i = 2:sz
        suma = 0
        for j = 1:i-1
            suma = suma + (A(:,i)'*Q(:,j))*Q(:,j)
        end
        Q(:,i) = A(:,i) - suma
        V(i) = norm(Q(:,i))
        Q(:,i) = Q(:,i)/V(i)  
    end
    
    R = diag(V)
    
    for i = 1:sz
        for j = i+1:sz
           R(i,j) = A(:,j)'*Q(:,i) 
        end
    end
endfunction

function [p, err] = minimoscuad(xi, y, gr)
    szy = size(y,1);
    n = gr+1;
    A = eye(szy, n);
    
    for j = 1:n
        for i = 1:szy
            A(i, j) = xi(i)**(j-1);    
        end
    end
    
    [Q,R] = FactoRQ(A);
    
    b = Q'*y
    
    sz = size(R,1)
    
    x(sz) = b(sz)/R(sz,sz)
    
    for i = 1:sz-1 //ResUELVE QR
        suma = 0
        for j = 1:i
           suma = suma + x(sz-j+1)*R(sz-i,sz-j+1)   
        end
        x(sz-i) = (b(sz-i)-suma)/R(sz-i,sz-i) 
    end
    
    p = poly(x, 'x', "coeff")
    
    E = A*x-y;
    err = E'*E;
    
endfunction


//ej9
function y = ej9(x)
    y = 1/(1+x.^2)
endfunction

function y = Interej9(pto)
    X2 = [1 0.500000;2 0.200000]
    X4 = [1 0.500000;2 0.200000;3 0.100000;4 0.058824]
    X6 = [1 0.500000;2 0.200000;3 0.100000;4 0.058824;5 0.038462;6 0.027027]
    X10 = [1 0.500000;2 0.200000;3 0.100000;4 0.058824;5 0.038462;6 0.027027;7 0.020000;8 0.015385;9 0.012195;10 0.009901]
    X14 = [1 0.500000;2 0.200000;3 0.100000;4 0.058824;5 0.038462;6 0.027027;7 0.020000;8 0.015385;9 0.012195;10 0.009901;11 0.008197;12 0.006897;13 0.005882;14 0.005076]
    y = InterpolacionNewton(X14,pto)
endfunction

function y = resta(x)
    y = abs(ej9(x) - Interej9(x))
endfunction


//ej10
function y = ej10(x)
    y = e.^x
endfunction


function ploty(fn,l,in,r) // Funcion, Limite Izq, Intervalo, Limite Der
//    xdel(winsid());
    x = [l:in:r];
//    n = size(x);
//   yy = zeros(1,n(2));
//   plot(x,yy)
    plot(x,fn)
    a = gca();
    a.auto_scale = "off"; 
endfunction
