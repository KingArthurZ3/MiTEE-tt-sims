function w_til = til(w)
% this function is to simplify the caculation. til is the matrix tilde operator in 
% linear algebra. it transfers a cross production to a matrix multiplication. 
  w1 = w(1);
  w2 = w(2); 
  w3 = w(3);
  
  w_til = [0 -w3 w2;
           w3 0 -w1;
           -w2 w1 0];
end