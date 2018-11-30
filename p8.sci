function y = reglaTrapecio(fn, a, b)
    h = b-a;
    y = h/2 * (fn(a) + fn(b));
endfunction

function y = reglaSimpson(fn, a, b)
    h = (b-a)/2;
    med = (b+a)/2;
    y = h/3 * (fn(a) + 4*fn(med) + fn(b)); 
endfunction

function y = metodoCompTrapecio(fn, a, b, n)
    h = (b-a)/n;
    suma = 0;
    for i = 0:n
        if(i == 0 | i == n)
            suma = suma + fn(a + i*h);
        else
            suma = suma + 2*fn(a + i*h);
        end
    end
    y = suma * h/2;
endfunction

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
