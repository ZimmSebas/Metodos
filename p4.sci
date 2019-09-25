// Ejercicio 1
//Eliminación de Gauss con pivoteo parcial

function [P,L,U] = egpp(A)
    U = A;
    m = size(A, 1);
    L = eye(m,m);
    P = eye(m,m);
    //A(:,[1 3]) = A(:,[3 1]) <- permutar
    for k = 1:m-1
      ind = find(abs(U(k:m,k)) == max(abs(U(k:m,k))),1) //Inicio pivoteo parcial
      ind = ind + (k-1)
      U([k ind],k:m) = U([ind k],k:m)
      L([k ind],1:k-1) = L([ind k],1:k-1)
      P([k ind],:) = P([ind k],:)                       //Fin pivoteo parcial
      for j = k+1:m
          L(j,k) = U(j,k) / U(k,k)
          U(j,k:m) = U(j,k:m) - L(j,k)*U(k,k:m)
      end 
    end
endfunction

function x = Gesolver(A,b)
    [P, L, U] = egpp(A)
    sz = size(L,1)
    b = P*b
    c(1) = b(1)
    
    for i = 2:sz
        suma = 0
        for j = 1:i-1
           suma = suma + c(j)*L(i,j)   
        end
        c(i) = b(i)-suma 
    end
    
    x(sz) = c(sz)/U(sz,sz)

    for i = 1:sz-1
        suma = 0
        for j = 1:i
           suma = suma + x(sz-j+1)*U(sz-i,sz-j+1)   
        end
        x(sz-i) = (c(sz-i)-suma)/U(sz-i,sz-i) 
    end    
endfunction

function [L,U] = dolittle(A)
    sz = size(A,1)
    L = eye(sz,sz)
    U = zeros(sz,sz)
    
    U(1,:) = A(1,:)
    for k = 1:sz
           
        for j = k:sz
            suma = 0
            for m = 1:k-1
                suma = suma + L(k,m)*U(m,j)
            end
            U(k,j)=A(k,j) - suma
            
            for i = k+1:sz
                suma = 0
                for m = 1:k-1
                    suma = suma + L(i,m)*U(m,k)
                end    
                L(i,k)=(A(i,k) - suma)/U(k,k)
            end
        end
    end
endfunction

function x = dlsolver(A,b)
    sz = size(A,1)
    [L,U] = dolittle(A)
    c(1) = b(1)
    
    for i = 2:sz
        suma = 0
        for j = 1:i-1
           suma = suma + c(j)*L(i,j)   
        end
        c(i) = b(i)-suma 
    end
    
    x(sz) = c(sz)/U(sz,sz)

    for i = 1:sz-1
        suma = 0
        for j = 1:i
           suma = suma + x(sz-j+1)*U(sz-i,sz-j+1)   
        end
        x(sz-i) = (c(sz-i)-suma)/U(sz-i,sz-i) 
    end    
    
endfunction

function [U, ind] = Cholesky(A)
    eps = 1.0e-8
    n = size(A,1)
    U = zeros(n,n)
    for k = 1:n
        t = A(k,k) - U(1:k,k)'*(U(1:k,k))
        if (t <= eps)
            mprintf("Matriz no definida positiva.\n")
            ind = 0
            return
        end
        if (A <> A')
            mprintf("Matriz no simetrica. \n")
            ind = 0
            return 
        end
        disp(t);
        U(k,k)= sqrt(t)
        for j = k+1:n
            U(k,j) = ( A(k,j) - U(1:k,k)'*U(1:k,j) )/U(k,k)
        end
    end
    ind = 1
endfunction

function x = chsolver(A,b)
    [U,ind] = Cholesky(A)
    if (ind == 1)
        c(1) = b(1)
        L = U'
        
        for i = 2:sz
            suma = 0
            for j = 1:i-1
               suma = suma + c(j)*L(i,j)   
            end
            c(i) = b(i)-suma 
        end
        
        x(sz) = c(sz)/U(sz,sz)
    
        for i = 1:sz-1
            suma = 0
            for j = 1:i
               suma = suma + x(sz-j+1)*U(sz-i,sz-j+1)   
            end
            x(sz-i) = (c(sz-i)-suma)/U(sz-i,sz-i) 
        end  
    else
        mprintf("Error, no se puede \n");
    end
endfunction
