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

function y=cx1(x)
    y=-sqrt(2*x-x**2)
endfunction

function y=dx1(x)
    y=sqrt(2*x-x**2)
endfunction

function z=uno(x,y)
    z=1
endfunction

function y = ln(x)
    y = log(x);
endfunction

function y = fb(x)
    y = x**(1/3);
endfunction

function y = fc(x)
    y = sin(x)**2;
endfunction

function y = dosb(x)
    y = x**3;
endfunction

function y = dosf(x)
    y = x**2 * %e**x;
endfunction
