funcprot(0);

old = 'f';
n = 20
for i=1:n
    new = 'd'+string(i)+'f';
    deff('y='+new+'(x)','y=numderivative('+old+',x,0.1)');
    old=new;
end

function y=cuad(x)
    y=x**2
endfunction

function y=f(x)
    y=%e**x
endfunction

function y = derivar(f,x,n,h)
    funcprot(0);
    if(n==0) 
        y = f(x);
    else
        y =( (f(x+h)-f(x-h)) / 2*h);
    end
endfunction

function  y = errorcalc(a,b)
    y(1) = abs(a - b)
    y(2) = abs(a - b)/abs(a)
    mprintf("error absoluto %0.15f \n", y(1))
    mprintf("error relativo %0.15f \n", y(2))
endfunction

function comp(f,x,n,h)
    funcprot(0);
    approxder=derivar(f,x,n,h)
    approxnumder=numderivative(f,x,h,n)
    mprintf("derivar: %f \n",approxder)
    mprintf("numderiv: %f \n",approxnumder)
    mprintf("error: %f \n",approxder-approxnumder)
endfunction

function y = reverse(arr)
    funcprot(0);
    n = length(arr)
    y = (1:n)
    for i = 1:n
        y(n+1-i) = arr(i)
    end
endfunction

function y = taylor(f,n,v0,v) //Funcion, numero de derivadas, punto inicial, punto a ver.
    funcprot(0);
    coeficientes = (1:n)
    coeficientes(1) = 0
    for j = 1:n
        coeficientes(j+1)  = der(v0,j)/factorial(j)
        mprintf("coef: deriv: %f fact: %f  coef: %f\n ",der(v0,j),factorial(j),coeficientes(j+1));
    end
    y = horner(reverse(coeficientes),v-v0) + f(v0)
endfunction

function y = der(x,k)
    funcprot(0);
    deff('y=foo(x)','y=d'+string(k)+'f(x)');
    y = foo(x)
endfunction

function y = horner(arr,x)
    funcprot(0);
    n = length(arr);
    y = arr(1);
    for j = 2:n
        y = y*x + arr(j);
    end
endfunction

function y = hornerder(arr,x) //arreglo = an + an-1 + an-2...
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
endfunction  // retorna el resultado y el resultado de la derivada
