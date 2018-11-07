function [P,L,U] = egpp(A)
    U = A;
    m = size(A, 1);
    L = eye(m,m);
    P = eye(m,m);
    //A(:,[1 3]) = A(:,[3 1]) <- permutar
    for k = 1:m-1
      ind = find(U(k:m,k) == max(U(k:m,k)),1)
      ind = ind + (k-1)
      U([k ind],k:m) = U([ind k],k:m)
      L([k ind],1:k-1) = L([ind k],1:k-1)
      P([k ind],:) = P([ind k],:)
      for j = k+1:m
          L(j,k) = U(j,k) / U(k,k)
          U(j,k:m) = U(j,k:m) - L(j,k)*U(k,k:m)
      end 
    end
endfunction

function x = lusolver(A,b)
    [P, L, U] = egpp(A)
    sz = size(L,1)
    b = P*b
    c(1) = b(1)
    
    for i = 2:sz
        suma = 0
        for j = 1:i-1
           suma = suma + c(j)*L(i,j)   
        end
        c(i) = b(i)-suma 
    end
    
    x(sz) = c(sz)/U(sz,sz)

    for i = 1:sz-1
        suma = 0
        for j = 1:i
           suma = suma + x(sz-j+1)*U(sz-i,sz-j+1)   
        end
        x(sz-i) = (c(sz-i)-suma)/U(sz-i,sz-i) 
    end    
endfunction

//function [L,U] = dolittle(A,b)
//    sz = size(A,1)
//    L = eye(sz,sz)
//    U = eye(sz,sz)
//    
//    for k = 1:sz
//        for i = 1:sz
//            for j = 1:sz
//                L(i,k)
//    
//endfunction

