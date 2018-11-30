function y = reglaTrapecio(fn, a, b)
    h = b-a;
    y = h/2 * (fn(a) + fn(b));
endfunction

function y = reglaSimpson(fn, a, b)
    h = (b-a)/2;
    med = (b+a)/2;
    y = h/3 * (fn(a) + 4*fn(med) + fn(b)); 
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
