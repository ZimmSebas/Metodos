function y = f(x)
    y = 2^(x-1)
endfunction

function y = serie5(x, n)
    if(n == 0) y = x; return;end;
    y = iterar(2** (x-1), n-1);
    return;
endfunction


function y = serie6(x, n)
    if(n == 0) y = x; return;end;
    y = iterar(2** (x-1), n-1);
    return;
endfunction

function ploty(fn,l,in,r)
    xdel(winsid());
    x = [l:in:r];
    y = fn(x);
    n = size(x);
    yy = zeros(1,n(2));
    plot(x,yy)
    plot(x,y)
    a = gca();
    a.auto_scale = "off"; 
endfunction
