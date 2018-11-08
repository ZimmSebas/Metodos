//Codebook de Scilab!

//Indice

//1ra parte: Derivadas y Polinomio de Taylor.
//dkf(x) siendo k el numero de derivada de f(x), deriva k veces a f en x
//misraices(p) Función que calcula raices de polinomio p de grado 2
//Horner(arr,x) algoritmo de horner para expresar el resultado de un polinomio en x. Toma un arr de índices y un x.
//hornerder(arr,x). Algoritmo de horner que calcula p(x) y p'(x).
//taylor(f,n,v0,v). Aproxima funciones con un Polinomio de Taylor. Recibe Función, nro de derivadas, punto inicial, punto buscado. 

//2da parte: Métodos de cálculo de raices en funciones.
//Ejemplo de serie
//bsnewton(mini,maxi,eps,epsf) Método de la bisección, toma minimo, maximo, funcion, epsilon y epsilon funcion
//secante(fst,snd,func,eps,epsf) Método Secante, toma primer elemento, segundo, funcion, epsilon y epsilon funcion
//falsapp(a,b,fun,eps,epsf). Método de la falsa posición, toma minimo, maximo, funcion, epsilon y epsilon funcion
//newt_mult(fn, X, N) Método de Newton para una o varias variables. Toma una función fn, un vector X y N iteraciones.
//newt_mult_fin(fn, X, eps) Método de Newton para una o varias variables con cota de finalización. Recibe función, par de elementos y épsilon.
//puntofijo(fn,x,eps,cantmax). Método de aproximación por punto fijo, tomo una función fn, un x inicial, un épsilon de corte y una cantmax de iteraciones

//3ra parte: Resolución de sistemas de ecuaciones lineales con métodos directos.
// egpp(A). Realiza la factorización de Gauss con Pivoteo Parcial. Recibe matriz A
// Gesolver(A,b) Resuelve un sistema con Gauss con Pivoteo Parcial. Recibe matriz A y vector b
// dolittle(A) Realiza la factorización de Doolittle. Recibe matriz A
// dlsolver(A,b) Resuelve un sistema con Doolittle. Recibe matriz A y vector b
// Cholesky(A) Realiza la factorización de Cholesky. Recibe matriz A.
// chsolver(A,b) Resuelve un sistema con factorización Cholesky. Recibe matriz A y vector b

//4ta parte: Resolución de sistemas de ecuaciones lineales con métodos iterativos.
//jacobitest(A). Test de Jacobi para saber si converge a la solución. Recibe matriz A
//gausstest(A). Test de Gauss-Seidel para saber si converge a la solución. Recibe matriz A
//gausssolver(A,b,x0,eps) Método de resolución de Gauss-Seidel. Recibe matriz A, vector b, un x0 inicial y un épsilon.
//jacobisolver(A,b,x0,eps) Método de resolución de Jacobi. Recibe matriz A, vector b, un x0 inicial y un épsilon.
//mat_jacobi(A) Crea un Nor con la matriz de jacobi
//factorescala(Nor) Factor escala w de método SOR, recibe una matriz Nor
//SOR (A,b,x0,vv,eps) Método de resolución de SOR. Recibe matriz A, vector b, recibe un x0, un factor escala w y un épsilon.
//Gram_Schmidt(A) Factorización Gram Schmidt. Recibe matriz A. A debe tener columnas LI
//FactoRQ(A). //Factorización QR. Recibe una matriz A. A debe tener Columnas LI

//5ta parte: Adicionales.
//errorcalc. Calcula error absoluto y relativo  y los imprime. Recibe a y b.
//reverse: Algoritmo auxiliar para girar un arreglo.
//comp_metodos(N). Comparador de tiempos de métodos armado.
//ploty. Grafica una función. Recibe Función, Límite izquierdo, Intervalo, Límite Derecho



// -------------------------------------- ------------------------------------------
//              1ra parte: Errores, Derivadas y Polinomio de Taylor
// -------------------------------------- ------------------------------------------


//Genera nder derivadas de f(x) 

function y=f(x) //f(x) se usa en varias funciones!
    y=%e**x
endfunction

old = 'f';
nder = 20
for i=1:nder
    new = 'd'+string(i)+'f';
    deff('y='+new+'(x)','y=numderivative('+old+',x,0.1)');
    old=new;
end

//Función que calcula raices de polinomio p de grado 2
function r = misraices(p)
    c = coeff(p, 0)
    b = coeff(p, 1)
    a = coeff(p, 2)
    
    if(b < 0)
        r(1) = (2*c)/(-b + sqrt(b**2 - 4*a*c))
        r(2) = (-b + sqrt(b**2 - 4*a*c))/(2*a)
    else
        r(1) = (-b - sqrt(b**2 - 4*a*c))/(2*a)
        r(2) = (2*c)/(-b - sqrt(b**2 - 4*a*c))
    end
endfunction

// Algoritmo de que aproxima funciones con un Polinomio de Taylor
function y = taylor(f,n,v0,v) //Funcion, numero de derivadas, punto inicial, punto a ver.
    coeficientes = (1:n)
    coeficientes(1) = 0
    for j = 1:n
        coeficientes(j+1)  = der(v0,j)/factorial(j)
        mprintf("coef: deriv: %f fact: %f  coef: %f\n ",der(v0,j),factorial(j),coeficientes(j+1));
    end
    y = Horner(reverse(coeficientes),v-v0) + f(v0)
endfunction

function y = der(x,k) //Para la funcion derivadas
    deff('y=foo(x)','y=d'+string(k)+'f(x)');
    y = foo(x)
endfunction

//Horner de polinomios, recibe x y el polinomio en forma de arreglo
function y = Horner(arr,x) 
    n = length(arr);
    y = arr(1);
    for j = 2:n
        y = y*x + arr(j);
    end
endfunction

// Algoritmo de horner que retorna el resultado y el resultado de la derivada
function y = hornerder(arr,x) //arreglo = an + an-1 + an-2..., esta devuelve p(x) y p'(x)
    funcprot(0);
    n = length(arr);
    y(1) = arr(1);
    
    if (n>1) 
        y(2) = arr(2);
    end
    
    for j = 2:n
        y(1) = y(1)*x + arr(j)
        if (n>1 & j>2)
            y(2) = y(2)*x + arr(j); 
        end
    end
endfunction 



// -------------------------------------- ------------------------------------------
//              2da parte: Cálculo de raices en funciones
// -------------------------------------- ------------------------------------------

//Ejemplo de serie
function y = serie(x, c, n) //cotas x = 2 y c = 1, que onda lo de sqrt(z)?
    if(n == 0) y = x; return;end;
    y = serie((x+c*((x^2)-5)),c,n-1);
    return;
endfunction


//Método de Bisección, toma minimo, maximo, funcion, epsilon y epsilon funcion
function med = bisecc(mini,maxi,fun,eps,epsf) 
    
    if(fun(maxi).*fun(mini) > 0)
      error('Intervalos del mismo signo');
    end;
    m = ((mini+maxi)/2);
    while(maxi-m > eps | abs(fun(m)) > epsf )
      m = ((mini+maxi)/2);
      if((fun(maxi) < 0 & fun(m) > 0) | (fun(maxi) > 0 & fun(m) < 0))
        mini = m;
      else
        maxi = m;
      end;
    end;
    med = m;
endfunction;


//Método Secante, toma primer elemento, segundo, funcion, epsilon y epsilon funcion
function seca = secante(fst,snd,func,eps,epsf)
    if(func(fst).*func(snd) > 0)
       error('Intervalo incorrecto')
    end
    
    seca = snd - (func(snd) * (snd-fst )/(func(snd)-func(fst)))
    fst = snd
    snd = seca
    
    while(abs(func(seca)) > epsf | abs(snd-fst) > eps)
        seca = snd - (func(snd) * (snd-fst )/(func(snd)-func(fst)))
        fst = snd
        snd = seca
    end
endfunction

//Método de la falsa posición, toma minimo, maximo, funcion, epsilon y epsilon funcion
function c = falsapp(a,b,fun,eps,epsf) 
    
    if(fun(b).*fun(a) > 0)
      error('Intervalos del mismo signo');
    end;
    
    c = b - fun(b)*(b-a)/(fun(b)-fun(a)) 
    while(b-c > eps | abs(fun(c)) > epsf )
      if((fun(b) < 0 & fun(c) > 0) | (fun(b) > 0 & fun(c) < 0))
        a = c;
      else
        b = c;
      end;
      c = b - fun(b)*(b-a)/(fun(b)-fun(a))
    end;
endfunction;

//Método de Newton para una o varias variables. Toma una función fn, un vector X y N iteraciones.
function y = newt_mult(fn, X, N)
    Xn = X;

    mprintf("X0 = %f\n", Xn)
    for i = 1:N
      J = numderivative(fn, Xn);
      J = 1/J;
      y = Xn - J*fn(Xn);
      Xn = y
      //mprintf("X%d = %0.5f |-| %0.5f\n", i, Xn(1), Xn(2))
    end
endfunction

//Método de Newton para multivariables con cota de finalización. Recibe función, par de elementos y épsilon.
function y = newt_mult_fin(fn, X, eps)  //Función, punto inicial, error
    y = X;
    Xn = X;
    i = 0;
    mprintf("X0 = %f\n", Xn)
    while(norm(y-Xn) > eps | i == 0)      // Norma euclideana
      Xn = y;
      J = numderivative(fn, Xn);
      J = 1/J;
      y = Xn - J*fn(Xn);
      i = i+1;
      mprintf("X%d = %0.12f |-| %0.12f\n", i, y(1), y(2))
    end
    
    [Jac, Hes] = numderivative(fn, y, [] , 2, "blockmat");
    mprintf("Una matriz es definida positiva si sus autovalores son positivos.\n")
    mprintf("Autovalores del hessiano de fn: ")
    disp(spec(Hes))   
endfunction

//Método de punto fijo
function y = puntofijo(fn,x,eps,cantmax)
    y = fn(x)
    cont = 0
    while (abs(y-x) > eps & cont < cantmax)
        cont = cont + 1
        x = y
        y = fn(x)
    end
    mprintf("La cantidad de iteraciones fue %d\n",cont);
endfunction



// -------------------------------------- --------------------------------------------------------
//              3ra parte: Resolución de sistemas de ecuaciones lineales con métodos directos
// -------------------------------------- ---------------------------------------------------------




//Realiza la eliminación de Gauss con Pivoteo Parcial. Recibe A matriz
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

//Resuelve un sistema con Gauss con Pivoteo Parcial. Recibe matriz A y vector b
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

//Método factorización Doolittle. Recibe matriz A 
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

//Resuelve un sistema con factorización Doolittle. Recibe matriz A y vector b
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

//Método de factorización Cholesky. Recibe matriz A.
function [U, ind] = Cholesky(A)
    eps = 1.0e-8
    n = size(A,1)
    U = zeros(n,n)
    for k = 1:n
        t = A(k,k) - U(1:k-1,k)'*(U(1:k-1,k))
        if (t <= eps)
            mprintf("Matriz no definida positiva.\n")
            ind = 0
            return
        end
        disp(t);
        U(k,k)= sqrt(t)
        for j = k+1:n
            U(k,j) = ( A(k,j) - U(1:k-1,k)'*U(1:k-1,j) )/U(k,k)
        end
    end
    ind = 1
endfunction

//Resuelve un sistema con factorización Cholesky. Recibe matriz A y vector b
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





// -------------------------------------- -----------------------------------------------------
//              4ta parte: Resolución de sistemas de ecuaciones lineales en forma iterativa
// -------------------------------------- -----------------------------------------------------



//Test de Jacobi para saber si converge a la solución. Recibe matriz A
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

//Test de Gauss-Seidel para saber si converge a la solución. Recibe matriz A
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

//Método de resolución de Gauss-Seidel. Recibe matriz A, vector b, un x0 inicial y un épsilon.
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

//Método de resolución de Jacobi. Recibe matriz A, vector b, un x0 inicial y un épsilon.
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

//Crea un Nor con la matriz de jacobi
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

//Factor escala de método SOR, recibe un Nor
function vv = factorescala(Nor)
    radspec = max(abs(spec(Nor)))
    vv = 2/ (1 + sqrt(1 - radspec^2))
endfunction

//Método de resolución de SOR. Recibe matriz A, vector b, recibe un x0, un factor escala w y un épsilon.
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


//Factorización Gram Schmidt. Recibe matriz A
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

//Factorización QR. Recibe una matriz A.
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






// -------------------------------------- -----------------------------------------------------
//                                      5ta parte: Adicionales!
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

//Graficadora de funciones
function ploty(fn,l,in,r) // Funcion, Limite Izq, Intervalo, Limite Der
    //xdel(winsid());
    x = [l:in:r];
    y = fn(x);
    n = size(x);
    yy = zeros(1,n(2));
    plot(x,yy)
    plot(x,y)
    a = gca();
    a.auto_scale = "off"; 
endfunction

