make_adj_lower_tri <- function(B, q, lett) {
  
  # Initialization
  adj_lowert = matrix(0, q, q)
  my_list = list()
  # create a list such that each element my_list[[i]] contains the parents of node i
  for (i in 1:q){
    my_list[[i]] = which(B[,i]!=0)
  }
  
  npa = colSums(B)
  # Initialize the vector containing the "new" column order so that in the first position there is the column having
  # 0 parents
  
  order_def = unname(c(which(npa==0)))
  
  for (m in order_def){
    my_list[[m]] = 0
  }
  
  k = length(order_def)
  
  while (k<q) {
    for (i in 1:q){
      
      if (length(setdiff(my_list[[i]], order_def[1:k]))==0){
        
        order_def[k+1] = i
        k = k+1
        my_list[[i]] = 0 
        break
      }
    }
    
  }
  
  adj_lowert[1:q,1:q] = (B[rev(order_def), rev(order_def)])    # columns and row exchange
  dimnames(adj_lowert) = list( lett[rev(order_def)], lett[rev(order_def)]) # assign new labels
  
  adj_lowert
  
}

# the function computes the causal effect Y|do(X=x'), starting from a 
# matrix of edges coefficients: B_hat.
# NB: B_hat must be provided in LT format

ce_from_B <- function(B_lt, x.pos, y.pos) {
  
  ce_est = B_lt[x.pos, y.pos]
  j = 1
  
  B_hatj = B_lt
  
  while(j <= x.pos-1){
    
    B_hatj = B_hatj  %*% B_lt
    ce_est = ce_est + B_hatj[x.pos, y.pos]
    
    j = j+1
  }
  
  return(ce_est)
}

