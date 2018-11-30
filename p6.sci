//gershgorin(A). Muestra las cotas de los autovalores de A dada dicha matriz.

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

//potencia(A, z0, max_iter). Dada una matriz A y un vector z0 (estimación de posible autovector) y una cantidad de iteraciones max_iter
// aproxima al autovalor cuyo módulo es el radio espectral (mayor valor absoluto).
function rho = potencia(A, z0, max_iter)
    sz = size(z0, 1)
    for i = 1:max_iter
        w = A*z0
        //disp(w)    
        if(i <> max_iter)
            z = w / norm(w, %inf)
            z0 = z
        end
    end
    //Elegimos la componente de mayor valor absoluto
    k = 1
    for i = 2:sz
        if( abs(w(i)) > abs(w(k)))
            k = i
        end
    end
    rho = w(k) / z0(k)
endfunction
