//clc;
//clear;
xdel(winsid());

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

function y = f(x)
    y = (x**2)/4 - sin(x)
endfunction

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
