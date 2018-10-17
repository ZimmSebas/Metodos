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

function ploty(fn,l,in,r) // Funcion, Limite Izq, Intervalo, Limite Der
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
