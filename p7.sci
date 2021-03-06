//InterpolacionLagrange(X,pto). Interpolacion de Lagrange, toma un conjunto de pares de elementos (x,f(x)) y un pto a aproximar
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

//DDs(X). Calcula diferencias divididas de Newton dado un conjunto de pares de elemntos (x,f(x))
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

//InterpolacionNewton(X,pto). Interpolacion de Newton, toma un conjunto de pares de elementos X (x,f(x)) y un pto a aproximar
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

// FactoRQ(A) Factoriza una matriz A en Q,R (usado en mínimos cuadrados)
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

//minimoscuad(xi,y,gr): Toma un conjunto de puntos xi, y un conjunto de f(xi)=y, y un grado de polinomio. Aproxima un polinomio por minimos cuadrados y devuelve el error.
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


// Ejercicio 8

X = [4 102.56; 4.2 113.18; 4.5 130.11; 4.7 142.05; 5.1 167.53; 5.5 195.14; 5.9 224.87; 6.3 256.73; 6.8 299.5; 7.1 326.72];
xi = X(:,1);
y = X(:,2);
scatter(xi,y, "fill")
[p1, e1] = minimoscuad(xi, y, 1);
[p2, e2] = minimoscuad(xi, y, 2);
[p3, e3] = minimoscuad(xi, y, 3);
plot(xi, [horner(p1,xi) horner(p2,xi) horner(p3,xi)])

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
    y = %e.^x
endfunction

//Chevyshev(n). Calcula las raices de un polinomio de Chevyshev para n nodos interpolantes 
function rot = Chevyshev(n)
    for i = 1:n
        rot(i) = cos( ((2*i)-1) *%pi/ (2*n) )  
    end
endfunction

//polinomioInter(fn,n). Dada una funcion fn y un nro de nodos n, calcula la Interpolacion de Newton de las raices de un polinomio de Chevyshev
function pf = polinomioInter(fn,n)
    rot = Chevyshev(n)
    sz = size(rot,1)
    
    for i = 1:sz
        X(i,1) = rot(i)
        X(i,2) = fn(rot(i))
    end
    
    D = DDs(X)
    
    sz = size(D,1)
    
    p(1) = 1
    pf = p(1) * D(1,1)
    for j = 2:sz
        p(j) = p(j-1) * poly([-X(j-1,1),1],'x','coeff')
        pf = pf + p(j) * D(1,j)
    end
endfunction


//ej11

function y = ej11(x)
    y = cos(x)
endfunction

//Chevyshev(n,a,b). Calcula las raices de un polinomio de Chevyshev para n nodos interpolantes en un intervalo general a,b. 
function rot = ChevyshevGen(n,a,b)
    rot = Chevyshev(n)
    
    for i = 1:n
        rot(i) = (a+b+rot(i)*(b-a) )/ 2
    end
endfunction

//polinomioInterGen(fn,n,a,b). Dada una funcion fn y un nro de nodos n, y par de puntos a,b de rango; calcula la Interpolacion de Newton de las raices de un polinomio de Chevyshev para un intervalo general[a,b]
function pf = polinomioInterGen(fn,n,a,b)
    rot = ChevyshevGen(n,a,b)
    sz = size(rot,1)
    
    for i = 1:sz
        X(i,1) = rot(i)
        X(i,2) = fn(rot(i))
    end
    
    D = DDs(X)
    
    sz = size(D,1)
    
    p(1) = 1
    pf = p(1) * D(1,1)
    for j = 2:sz
        p(j) = p(j-1) * poly([-X(j-1,1),1],'x','coeff')
        pf = pf + p(j) * D(1,j)
    end
endfunction

//polyChevyshev(n) Crea el polinomio de Chevyshev de n nodos y calcula sus raices
function rot = polyChevyshev(n)
    T1(1) = 1
    T2(1) = 0
    T2(2) = 1
    for i = 1:n-1
        sz2 = size(T2,1)
        for j = 1:sz2+1
            if (j == 1)
                T3(j) = 0
            else
                T3(j) = T2(j-1)*2
            end
        end
        sz1 = size(T1,1)
        for j = 1:sz1
            T3(j) = T3(j) - T1(j)
        end
        
        T1 = T2
        T2 = T3
    end

    pc = poly(T3','x','coeff')
    rot = roots(pc)
endfunction

//ploty(fn,l,in,r). Plotea con funcion, limite izq, intervalo, limite der
function ploty(fn,l,in,r)
//    xdel(winsid());
    x = [l:in:r];
//    n = size(x);
//   yy = zeros(1,n(2));
//   plot(x,yy)
    plot(x,fn)
    a = gca();
    a.auto_scale = "off"; 
endfunction
