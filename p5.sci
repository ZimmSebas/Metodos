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
    //it = 1
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


// Ejercicio 3

function mat = ma3(n)
    mat = zeros(n,n)
    mat = mat + 2*eye(n,n)
    
    for i = 1:n
        for j = 1:n
            if(abs(i-j) == 1)
                mat(i,j) = -1
            end
        end
    end 
endfunction

function N = expl_gauss(A)
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
    
    Nor = (I-inv(N)*A)
    disp(Nor)
    
endfunction


// Ejercicio 4


function [P,L,U] = egpp(A)
    U = A;
    m = size(A, 1);
    L = eye(m,m);
    P = eye(m,m);
    //A(:,[1 3]) = A(:,[3 1]) <- permutar
    for k = 1:m-1
      ind = find(A(k:m,k) == max(A(k:m,k)),1)
      ind = ind + (k-1)
      U([k ind],k:m) = U([ind k],k:m)
      L([k ind],1:k-1) = L([ind k],1:k-1)
      P([k ind],:) = P([ind k],:)
      for j = k+1:m
          L(j,k) = U(j,k) / U(k,k)
          U(j,k:m) = U(j,k:m) - L(j,k)*U(k,k:m)
      end 
    end
endfunction

function x = lusolver(A,b)
    [P, L, U] = egpp(A);
    sz = size(L,1);
    b = P*b;
    c(1) = b(1);
    
    for i = 2:sz
        suma = 0;
        for j = 1:i-1
           suma = suma + c(j)*L(i,j)   
        end
        c(i) = b(i)-suma 
    end
    
    x(sz) = c(sz)/U(sz,sz);

    for i = 1:sz-1
        suma = 0
        for j = 1:i
           suma = suma + x(sz-j+1)*U(sz-i,sz-j+1)   
        end
        x(sz-i) = (c(sz-i)-suma)/U(sz-i,sz-i) 
    end    
endfunction


function [x,t] = res_cuatro(N)
    A = 8*eye(N,N) + 2*diag(ones(N-1, 1), 1) + 2*diag(ones(N-1, 1), -1) + diag(ones(N-3,1), 3) + diag(ones(N-3,1), -3)
    b = ones(N,1)
    tic();
    lusolver(A,b);
    t(1)=toc();
    
    eps = 10^-6;
    x0 = zeros(N,1);
    
    tic();
    gausssolver(A,b,x0,eps);
    t(2) = toc();
    
    eps = 10^-12;
    
    tic();
    gausssolver(A,b,x0,eps);
    t(3) = toc();
    
endfunction


// Ejercicio 5


function Nor = mat_jacobi(A)
    sz = size(A,1)
    N = eye(sz,sz)
    I = N
    di = diag(A)
    for i = 1:sz
        N(i,i) = N(i,i) * di(i)
    end
    
    Nor = (I-inv(N)*A)
    
endfunction

function vv = factorescala(Nor)
    radspec = max(abs(spec(Nor)))
    vv = 2/ (1 + sqrt(1 - radspec^2))
endfunction


function x = SOR(A, b, x0, vv, eps)
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
            x(i) = (1 - vv)*xant(i) + vv/(A(i,i))*(b(i)-suma)
        end
     cont = cont + 1
    end
    
    mprintf("Cantidad de iteraciones: %d\n",cont);
    
endfunction


// Bonus track: Factorización QR

function GS = Gram_Schmidt(A)   // A debe tener columnas LI
    sz = size(A,2)
    Qu(:,1) = A(:,1)/norm(A(:,1))
    for i = 2:sz
        suma = 0
        for j = 1:i-1
            suma = suma + (A(:,j)'*Qu(:,j))*Qu(:,j)
        end
        Qu(:,i) = A(:,i) - suma
        Qu(:,i) = Qu(:,i)/norm(Qu(:,i))  
    end
endfunction

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
