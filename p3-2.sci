clc;
clear;
xdel(winsid());


function y = h(x)
    y = log(x)-x+1;
endfunction;

function y = g(x)
    y = %e^-x-x**4;
endfunction;

function y = f(x)
    y = sin(x)-(x^2)/2;
endfunction;

function med = bsrec(mini,maxi,fun,eps,epsf)
    m = ((mini+maxi)/2);
    if(maxi-m < eps & abs(fun(m)) < epsf )
        med = m;
    else 
      if(fun(maxi).*fun(m) < 0)
        med = bsrec(m,maxi,fun,eps,epsf);
      elseif(fun(m) == 0)
        med = m;
      else
        med = bsrec(mini,m,fun,eps,epsf);
      end
    end
endfunction

function med = bs(mini,maxi,fun,eps,epsf)
    
    if(fun(maxi).*fun(mini) > 0)
      error('Intervalos del mismo signo');
    end;
    
    m = ((mini+maxi)/2);
    while(maxi-m > eps | abs(fun(m)) > epsf )
      m = ((mini+maxi)/2);
      if((fun(maxi) < 0 & fun(m) > 0) | (fun(maxi) > 0 & fun(m) < 0))
        mini = m;
    //      elseif(fun(m) == 0)
    //        m = m
      else
        maxi = m;
      end;
    end;
    med = m;
endfunction;


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

