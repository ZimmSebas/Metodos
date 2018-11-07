//Método de Newton, toma minimo, maximo, funcion, epsilon y epsilon funcion
function med = bsnewton(mini,maxi,fun,eps,epsf) 
    
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

function y = serie5(x, n)
    if(n == 0) y = x; return;end;
    y = serie5(2** (x-1), n-1);
    return;
endfunction


function y = serie6(x, c, n) //cotas x = 2 y c = 1, que onda lo de sqrt(z)?
    if(n == 0) y = x; return;end;
    y = serie6((x+c*((x^2)-5)),c,n-1);
    return;
endfunction

function y = ide(x)
    y = x
endfunction

// Ejercicio 8
// ge(f,x,n) es la funcion generica para iterar sobre una funcion dada
// g1 a g4 son las funciones del ejercicio

function y = g1(x)
    y = %e^x/3
endfunction

function y = g2(x)
    y = (%e^x - x)/2
endfunction

function y = g3(x)
    y = log(3*x)
endfunction

function y = g4(x)
    y = %e^x - 2*x
endfunction

function y = comp(f, x)
    y = f(x)
endfunction


function y = ge(x, n)    // Funcion, Punto, Cantidad iteraciones
    if(n == 0) y = x; return;end;
    y = ge(g3(x), n-1);
    return; 
endfunction


// Ejercicio 9

function n = fnueve(X)
    n = [1 + X(1)^2 - X(2)^2 + %e^X(1)*cos(X(2)); 2*X(1)*X(2) + %e^X(1)*sin(X(2))];
endfunction

function y = newt_mult(fn, X, N)
    Xn = X;

    mprintf("X0 = %f\n", Xn)
    for i = 1:N
      J = numderivative(f, Xn);
      J = 1/J;
      y = Xn - J*fn(Xn);
      Xn = y
      mprintf("X%d = %0.5f |-| %0.5f\n", i, Xn(1), Xn(2))
    end
endfunction

// Ejercicio 10

function n = fdiez(X)
    n = [X(1)^2 + X(1)*X(2)^3 - 9; 3*X(1)^2*X(2) - 4 - X(2)^3];
endfunction

Xa = [1.2; 2.5]
Xb = [-2; 2.5]
Xc = [-1.2; -2.5]
Xd = [2; -2.5]


// Ejercicio 11

function y = fonce(X)
    y = 2*X(1) + 3*X(2)^2 + %e^(2*X(1)^2 + X(2)^2)
endfunction


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
