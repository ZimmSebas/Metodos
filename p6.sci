function gershgorin(A)
    sz = size(A, 1);
    
    for i = 1:sz
        suma = 0;
        for j = 1:sz
            if (i <> j)
                suma = suma + abs(A(i,j));
            end;
        end;
            mprintf("|lambda - %f| <= %f\n", A(i,i), suma);
    end;
endfunction


// poly([A], "x") -> polinomio caracteristico de la matriz A 
//(det(lambda*I - A) = p(lambda))

// Dado un polinomio mónico (o normalizado) p(lambda) una matriz
// A es compañera del polinomio si det(lambda*I - A) = p(lambda),
// es decir, si el polinomio característico de A es p.


function ej_3(A)
   sz = size(A, 1);
   for k = 0:10
       mprintf("k = %d\n", k)
       
       A(sz,sz) = 1 + 0.1*k;
       p = poly([A], "x");
       x = roots(p);
       disp(x)
       av = spec(A);
       disp(av)
       gershgorin(A);
   end
endfunction

function rho = potencia(A, x0, max_iter)
    sz = size(x0, 1)
    for i = 1:max_iter
        x1 = A*x0
        disp(x1)    
        if(i <> max_iter)
            x1 = x1 / norm(x1, %inf)
            x0 = x1
        end
    end
    //Elegimos la componente de mayor valor absoluto
    k = 1
    for i = 2:sz
        if( abs(x1(i)) > abs(x1(k)))
            k = i
        end
    end
    rho = x1(k) / x0(k)
endfunction
