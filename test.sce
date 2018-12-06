function y = intBidimensionalS(fn,a,b,cx,dx,n,m)
    deff('z=aux1(y)','z=fn(a,y)')
    deff('z=aux2(y)','z=fn(b,y)')
    temp= SimpsonCompuesto(cx(a),dx(a),aux1,m) + SimpsonCompuesto(cx(b),dx(b),aux2,m)
    
    h = (b-a)/n
    for i=1:n-1
        xi = a+i*h
        deff('z=aux(y)','z=fn(xi,y)')
        if modulo(i,2) == 0 then
            temp = temp + 2*(SimpsonCompuesto(cx(xi),dx(xi),aux,m))
        else
            temp = temp + 4*(SimpsonCompuesto(cx(xi),dx(xi),aux,m))
        end
    end
    y = (h/3) * temp
endfunction

function y = SimpsonCompuesto(a,b,f,n)
    h = (b-a)/n
    
    for j=1:n-1
        X(j)= a+j*h
    end
    disp(X);
    
    resultado = 0
    
    for j=1:n-1
        if modulo(j,2) == 0 then
            resultado = resultado + 2*f(X(j))
        else
            resultado = resultado + 4*f(X(j))
        end
    end
    
    y = (h/3) * (f(a) + resultado + f(b))
    
endfunction

function y=cx1(x)
    y=-sqrt(2*x-x**2)
endfunction

function y=dx1(x)
    y=sqrt(2*x-x**2)
endfunction

function y=uno(x,y)
    y=1
endfunction

disp(intBidimensionalS(uno,0,2,cx1,dx1,10,10))
