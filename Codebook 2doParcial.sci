//Codebook de Scilab pt 2!

//Indice

//1ra parte: Aproximacion de Autovalores y Autovectores

//gershgorin(A). Muestra las cotas de los autovalores de A dada dicha matriz.
//potencia(A, z0, max_iter). Dada una matriz A y un vector z0 (estimación de posible autovector) y una cantidad de iteraciones max_iter aproxima al autovalor cuyo módulo es el radio espectral (mayor valor absoluto).

//2da parte: Interpolación y mínimos cuadrados

//InterpolacionLagrange(X,pto). Interpolacion de Lagrange, toma un conjunto de pares de elementos (x,f(x)) y un pto a aproximar
//DDs(X). Calcula diferencias divididas de Newton dado un conjunto de pares de elemntos (x,f(x))
//InterpolacionNewton(X,pto). Interpolacion de Newton, toma un conjunto de pares de elementos X (x,f(x)) y un pto a aproximar
// FactoRQ(A) Factoriza una matriz A en Q,R (usado en mínimos cuadrados)
//minimoscuad(xi,y,gr): Toma un conjunto de puntos xi, y un conjunto de f(xi)=y, y un grado de polinomio. Aproxima un polinomio por minimos cuadrados y devuelve el error.
//Chevyshev(n). Calcula las raices de un polinomio de Chevyshev para n nodos interpolantes 
//polinomioInter(fn,n). Dada una funcion fn y un nro de nodos n, calcula la Interpolacion de Newton de las raices de un polinomio de Chevyshev

//3ra parte: Metodos de Integración

//reglaTrapecio(fn,a,b). Aplica método del trapecio a una función fn de a hasta b.
//reglaSimpson(fn,a,b). Aplica método del Simpson a una función fn de a hasta b.
//metodoCompTrapecio(fn,a,b,n). Aplica método del trapecio a una función fn de a hasta b en n intervalos.
//metodoCompSimpson(fn,a,b,n). Aplica método del Simpson a una función fn de a hasta b en n intervalos.
//reglaTrapecioExt(fn,x1,x2,y1,y2). Aplica método del trapecio a una función fn de x1 a x2 y de y1 a y2.
//IntDosS(f,a,b,cx,dx,n,m). Aplica método de Simpson a una función fn de dos variables desde a hasta b, y desde cx a dx en N iteraciones de x e M iteraciones de Y.
//IntDosT(f,a,b,cx,dx,n,m). Aplica método de trapecio a una función fn de dos variables desde a hasta b, y desde cx a dx en N iteraciones de x e M iteraciones de Y.

//4ta parte: Adicionales!

//errorcalc(a,b) Calcula el error absoluto y el error relativo
//reverse(arr). Algoritmo que da vuelta los indices de un arreglo
//ploty(fn,l,in,r). Plotea con funcion, limite izq, intervalo, limite der
//comp_metodos(N) Comparador de tiempos de métodos armado.


// -------------------------------------- ------------------------------------------
//              1ra parte: Aproximacion de Autovalores y Autovectores
// -------------------------------------- ------------------------------------------

//gershgorin(A). Muestra las cotas de los autovalores de A dada dicha matriz.

function gershgorin(A)
    sz = size(A, 1);
    
    for i = 1:sz
        suma = 0;
        for j = 1:sz
            if (i <> j)
                suma = suma + abs(A(i,j));
            end;
        end;
            mprintf("|lambda - %f| <= %f\n", A(i,i), suma);
    end;
endfunction


// poly([A], "x") -> polinomio caracteristico de la matriz A 
//(det(lambda*I - A) = p(lambda))

// Dado un polinomio mónico (o normalizado) p(lambda) una matriz
// A es compañera del polinomio si det(lambda*I - A) = p(lambda),
// es decir, si el polinomio característico de A es p.


function ej_3(A)
   sz = size(A, 1);
   for k = 0:10
       mprintf("k = %d\n", k)
       
       A(sz,sz) = 1 + 0.1*k;
       p = poly([A], "x");
       x = roots(p);
       disp(x)
       av = spec(A);
       disp(av)
       gershgorin(A);
   end
endfunction

//potencia(A, z0, max_iter). Dada una matriz A y un vector z0 (estimación de posible autovector) y una cantidad de iteraciones max_iter
// aproxima al autovalor cuyo módulo es el radio espectral (mayor valor absoluto).
function rho = potencia(A, z0, max_iter)
    sz = size(z0, 1)
    for i = 1:max_iter
        w = A*z0
        //disp(w)    
        if(i <> max_iter)
            z = w / norm(w, %inf)
            z0 = z
        end
    end
    //Elegimos la componente de mayor valor absoluto
    k = 1
    for i = 2:sz
        if( abs(w(i)) > abs(w(k)))
            k = i
        end
    end
    rho = w(k) / z0(k)
endfunction



// -------------------------------------- ------------------------------------------
//              2da parte: Interpolación y mínimos cuadrados
// -------------------------------------- ------------------------------------------

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

// -------------------------------------- -----------------------------------------------------
//                          3ra parte: Metodos de Integración
// -------------------------------------- -----------------------------------------------------

//reglaTrapecio(fn,a,b). Aplica método del trapecio a una función fn de a hasta b.
function y = reglaTrapecio(fn, a, b)
    h = b-a;
    y = h/2 * (fn(a) + fn(b));
endfunction


//reglaSimpson(fn,a,b). Aplica método del Simpson a una función fn de a hasta b.
function y = reglaSimpson(fn, a, b)
    h = (b-a)/2;
    med = (b+a)/2;
    y = h/3 * (fn(a) + 4*fn(med) + fn(b)); 
endfunction

//metodoCompTrapecio(fn,a,b,n). Aplica método del trapecio a una función fn de a hasta b en n intervalos.
function y = metodoCompTrapecio(fn, a, b, n)
    h = (b-a)/n;
    suma = 0;
    for i = 0:n
        asd = a + i*h;
        if(i == 0 | i == n)
            suma = suma + fn(asd);
        else
            suma = suma + 2*fn(asd);
        end
    end
    y = suma * h/2;
endfunction

//metodoCompSimpson(fn,a,b,n). Aplica método del Simpson a una función fn de a hasta b en n intervalos.
function y = metodoCompSimpson(fn, a, b, n)
    h = (b-a)/n;
    suma = 0;
    for i = 0:n
        if(i == 0 | i == n)
            suma = suma + fn(a + i*h);
        else
            if (pmodulo(i,2) == 1)
                suma = suma + 4*fn(a + i*h);
            else
                suma = suma + 2*fn(a + i*h);
            end    
        end
    end
    y = suma * h/3;
endfunction

//reglaTrapecioExt(fn,x1,x2,y1,y2). Aplica método del trapecio a una función fn de x1 a x2 y de y1 a y2.
function y = reglaTrapecioExt(fn,x1,x2,y1,y2)
    h = (y2-y1)*(x2-x1)/4;
    y = h * (fn(x1,y1)+fn(x2,y1)+fn(x1,y2)+fn(x2,y2));
endfunction

//IntDosS(f,a,b,cx,dx,n,m). Aplica método de Simpson a una función fn de dos variables desde a hasta b, y desde cx a dx en N iteraciones de x e M iteraciones de Y.
function y = IntDosS(f,a,b,cx,dx,n,m)
    deff('z=aux1(y)','z=f(a,y)')
    deff('z=aux2(y)','z=f(b,y)')
    temp = metodoCompSimpson(aux1,cx(a),dx(a),m) + metodoCompSimpson(aux2,cx(b),dx(b),m)
    
    h = (b-a)/n
    for i=1:n-1
        xi = a+i*h
        deff('z=aux(y)','z=f(xi,y)')
        if pmodulo(i,2) == 0 then
            temp = temp + 2*(metodoCompSimpson(aux,cx(xi),dx(xi),m))
        else
            temp = temp + 4*(metodoCompSimpson(aux,cx(xi),dx(xi),m))
        end
    end
    y = (h/3) * temp
endfunction

//IntDosT(f,a,b,cx,dx,n,m). Aplica método de trapecio a una función fn de dos variables desde a hasta b, y desde cx a dx en N iteraciones de x e M iteraciones de Y.
function y = IntDosT(f,a,b,cx,dx,n,m)
    deff('z=aux1(y)','z=f(a,y)')
    deff('z=aux2(y)','z=f(b,y)')
    temp= (metodoCompTrapecio(aux1,cx(a),dx(a),m)/2) + (metodoCompTrapecio(aux2,cx(b),dx(b),m)/2)
    h = (b-a)/n
    for i=1:n-1
        xi = a+i*h
        deff('z=aux(y)','z=f(xi,y)')
        temp = temp + (metodoCompTrapecio(aux,cx(xi),dx(xi),m))
    end
    y = h * temp
endfunction


// -------------------------------------- -----------------------------------------------------
//                          4ta parte: Adicionales!
// -------------------------------------- -----------------------------------------------------

// Comparador de tiempos de métodos armado.
function [x,t] = comp_metodos(N)
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


// Calcula el error absoluto y el error relativo
function  y = errorcalc(a,b) 
    y(1) = abs(a - b)
    y(2) = abs(a - b)/abs(a)
    mprintf("error absoluto %0.15f \n", y(1))
    mprintf("error relativo %0.15f \n", y(2))
endfunction

// Algoritmo que da vuelta los indices de un arreglo
function y = reverse(arr) 
    n = length(arr)
    y = (1:n)
    for i = 1:n
        y(n+1-i) = arr(i)
    end
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

