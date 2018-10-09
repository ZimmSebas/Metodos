
//funcprot(0)

function y = horner(arr,x) //arreglo = an + an-1 + an-2...
    n = length(arr);
    y(1) = arr(1);
    
    if (n>1) 
        y(2) = arr(2);
    end
    
    for j = 2:n
        y(1) = y(1)*x + arr(j)
        if ( n>1 & j>2)
            y(2) = y(2)*x + arr(j); 
        end
    end
endfunction  // retorna el resultado y el resultado de la derivada

