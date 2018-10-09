old = 'f';
n = 20
for i=1:n
    new = 'd'+string(i)+'f';
    deff('y='+new+'(x)','y=numderivative('+old+',x,0.1)');
    old=new;
end

function y=derivar2(f,x,n,h)
    if n==0 then y=f(x)
    else y=(derivar2(f,x+(h/2),h)-derivar2(f,x-(h/2),h))/h
    end
endfunction

function y=f(x)
    y=x**2
endfunction

function y=powE(x)
    y=%e**x
endfunction

function comp(f,x,n,h)
    approxder=derivar(f,x,n,h)
    approxnumder=numderivative(f,x,h,n)
    mprintf("derivar: %f \n",approxder)
    mprintf("numderiv: %f \n",approxnumder)
    mprintf("error: %f \n",approxder-approxnumder)
endfunction
//mprintf("%f",derivar(powE,1,2,0.01))
//comp(powE,1,4,0.01)

function y = horner(arr,x)
    n = length(arr);
    y = arr(1);
    for j = 2:n
        y = y*x + arr(j);
    end
endfunction

function y = reverse(arr)
    n = length(arr)
    y = (1:n)
    for i = 1:n
        y(n+1-i) = arr(i)
    end
endfunction

function y = taylor(f,n,v0,v)
    coeficientes = (1:n)
    for j = 1:n
//        coeficientes(j)  = derivar(f,v0,j,0.001)/factorial(j)
        coeficientes(j)  = der(v0,j)/factorial(j)

        mprintf("coef: %f %f %f\n ",der(v0,j),factorial(j),coeficientes(j));
    end
    y = horner(reverse(coeficientes),v-v0) + f(v0)

endfunction

function y = der(x,k)
    deff('y=foo(x)','y=d'+string(k)+'f(x)');
    y = foo(x)
endfunction

//mprintf("deruiver + fact: %f\n ",derivar(trip,2,2,0.1)/factorial(2))
v0=5

res=taylor(trip,2,5,v0)
vres=trip(v0)
mprintf("\n resultado: %f",res)
mprintf("\n verdadero resultado: %f",vres)
mprintf("\n vError: %f",(abs(vres-res)/vres)*100)

//comp(trip,2,1,0.000001)
