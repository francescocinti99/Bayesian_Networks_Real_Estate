##CAPITOLO 3
#Set Working Directory
#setwd("~/Desktop/Tesi")

#Install and Load the packages
#install.packages("bnlearn")
library(bnlearn)
#install.packages("readxl")
library(readxl)
#install.packages("Rgraphviz")
library(Rgraphviz)

source("auxiliary_functions.R")


#Load the dataset for the analysis
data<- read_excel("DATASET TESI .xlsx")
#Rename variables from a to j
names(data) <- letters[1:ncol(data)]
str(data)
summary(data)

#Since bnlearn works both with continuous and discrete data, but the majority of the scoring methods require either continuous 
#variables or normally distributed variables, I employ the scale function to standardize the data to have mean 0 and variance 1.
data <- scale(data)

#Generate data as dataframe df
df <- data.frame(data)

# Use hill climbing algorithm with Gaussian assumption and visualize the network
blacklist <- data.frame(from = "h", to = "a")
hc_result <- hc(df, score = "bic-g", blacklist = blacklist) # inserire eventuali blacklist (frecce non consentite in una direzione)

## Sul DAG stimato
DAG = hc_result

#DAG è il DAG stimato su cui abbiamo applicato blacklist
plot(DAG)
graphviz.plot(DAG)

# Fit the network using the learned structure
fitted_network <- bn.fit(DAG, data = df)
q = length(fitted_network)

# Print the fitted network to see the parameters, to obtain coefficients given conditional dependencies
print(fitted_network)

lett = letters[1:ncol(data)]

#Creiamo una matrice B dove l'elemento B[u,v] è il coefficiente di Xu->Xv
# B[,v] ha tutte le frecce che arrivano in v
# gli 0 sono senza frecce, quelli diversi da 0 sono i parents di v
B = matrix(0,length(fitted_network),length(fitted_network))
colnames(B) = rownames(B) = lett

for(j in 1:q){
  pa_j = which(!is.na(match(lett, names(fitted_network[[j]]$coefficients[-1]))))
  B[pa_j,j] = fitted_network[[j]]$coefficients[-1]
}

B

# Funzione per rendere la matrice B lower triangular

#Rendiamo la matrice B lower triangular per poter applicare la formula
# adj_matrix: p x p adjacency matrix of a DAG
# p   = number of nodes
# lab = the labels of the nodes

adj.mat = amat(fitted_network)
B
B_lt = make_adj_lower_tri(B, 10, lett); 
B_lt


# Suppose now that we want to estimate the causal effect of an intervention of node i 
# on node cj
# we have to find the new positions of node i and j in B_lt

ce_ge = ce_from_B(B_lt, x.pos = which(colnames(B_lt) == "g"), y.pos = which(colnames(B_lt) == "e"))
ce_ge

ce_gc = ce_from_B(B_lt, x.pos = which(colnames(B_lt) == "g"), y.pos = which(colnames(B_lt) == "c"))
ce_gc

ce_jd = ce_from_B(B_lt, x.pos = which(colnames(B_lt) == "j"), y.pos = which(colnames(B_lt) == "d"))
ce_jd

ce_cg = ce_from_B(B_lt, x.pos = which(colnames(B_lt) == "c"), y.pos = which(colnames(B_lt) == "g"))
ce_cg

library(pcalg)
out_all_dags = pdag2allDags(as(B_lt, "matrix"))

# Funzione per calcolare l'effetto causale a partire da B
# data la posizione di x e y nella matrice B lower triangular
# do(X=x) su y
# output è effetto causale di X su Y

# return(gamma(DAG))


##

## Su ogni DAG nella sua classe di equivalenza

ce_eg = ce_from_B(B_lt, x.pos = which(colnames(B_lt) == "e"), y.pos = which(colnames(B_lt) == "g"))
ce_eg

library(pcalg)

hc_result_amat = amat(fitted_network)
CPDAG = dag2cpdag(as(hc_result_amat, "graphNEL"))

plot(CPDAG)

bn_cpdag = empty.graph(nodes = lett)
amat(bn_cpdag) = as(CPDAG, "matrix")


## 2

DAG2 = amat(DAG)

DAG2[1,5] = 0
DAG2[5,1] = 1
DAG2

bn_dag_2 = empty.graph(nodes = lett)
amat(bn_dag_2) = DAG2

## Ripeto gli steps precedenti, ossia:

fitted_network_2 <- bn.fit(bn_dag_2, data = df)

B2 = matrix(0,length(fitted_network_2),length(fitted_network_2))
colnames(B2) = rownames(B2) = lett

for(j in 1:q){
  pa_j = which(!is.na(match(lett, names(fitted_network_2[[j]]$coefficients[-1]))))
  B2[pa_j,j] = fitted_network_2[[j]]$coefficients[-1]
}

B2 # obs non c'è più il coefficiente di a -> e, ma di e -> a

adj.mat = amat(fitted_network_2)
B_lt_2 = make_adj_lower_tri(B2, 10, lett); 
B_lt_2

ce_eg_2 = ce_from_B(B_lt_2, x.pos = which(colnames(B_lt_2) == "e"), y.pos = which(colnames(B_lt_2) == "g"))
ce_eg_2

## 3

DAG3 = amat(DAG)

DAG3[5,8] = 0
DAG3[8,5] = 1
DAG3[1,5] = 0
DAG3[5,1] = 1

bn_dag_3 = empty.graph(nodes = lett)
amat(bn_dag_3) = DAG3

## Ripeto gli steps precedenti, ossia:

fitted_network_3 <- bn.fit(bn_dag_3, data = df)

B3 = matrix(0,length(fitted_network_3),length(fitted_network_3))
colnames(B3) = rownames(B3) = lett

for(j in 1:q){
  pa_j = which(!is.na(match(lett, names(fitted_network_3[[j]]$coefficients[-1]))))
  B3[pa_j,j] = fitted_network_3[[j]]$coefficients[-1]
}

B3

adj.mat = amat(fitted_network_3)
B_lt_3 = make_adj_lower_tri(B3, 10, lett); 
B_lt_3

ce_eg_3 = ce_from_B(B_lt_3, x.pos = which(colnames(B_lt_3) == "e"), y.pos = which(colnames(B_lt_3) == "g"))
ce_eg_3



## 4

DAG4 = amat(DAG)

DAG4[5,1] = 1
DAG4[1,8] = 0
DAG4[5,8] = 0
DAG4[5,7] = 0
DAG4[1,7] = 0

bn_dag_4 = empty.graph(nodes = lett)
amat(bn_dag_4) = DAG4

## Ripeto gli steps precedenti, ossia:

fitted_network_4 <- bn.fit(bn_dag_4, data = df)

B4 = matrix(0,length(fitted_network_4),length(fitted_network_4))
colnames(B4) = rownames(B4) = lett

for(j in 1:q){
  pa_j = which(!is.na(match(lett, names(fitted_network_4[[j]]$coefficients[-1]))))
  B4[pa_j,j] = fitted_network_4[[j]]$coefficients[-1]
}

B4

adj.mat = amat(fitted_network_4)
B_lt_4 = make_adj_lower_tri(B4, 10, lett); 
B_lt_4

ce_eg_4 = ce_from_B(B_lt_4, x.pos = which(colnames(B_lt_4) == "e"), y.pos = which(colnames(B_lt_4) == "g"))
ce_eg_4


##5
DAG5 = amat(DAG)

DAG5[5,1] = 0
DAG5[5,7] = 0
DAG5[1,8] = 0
DAG5[8,7] = 0
DAG5[5,8] = 1
DAG5[8,6] = 1
DAG5[6,3] = 1
DAG5[3,7] = 1
DAG5[6,9] = 0

bn_dag_5 = empty.graph(nodes = lett)
amat(bn_dag_5) = DAG5

## Ripeto gli steps precedenti, ossia:

fitted_network_5 <- bn.fit(bn_dag_5, data = df)

B5 = matrix(0,length(fitted_network_5),length(fitted_network_5))
colnames(B5) = rownames(B5) = lett

for(j in 1:q){
  pa_j = which(!is.na(match(lett, names(fitted_network_5[[j]]$coefficients[-1]))))
  B5[pa_j,j] = fitted_network_5[[j]]$coefficients[-1]
}

B5

adj.mat = amat(fitted_network_5)
B_lt_5 = make_adj_lower_tri(B5, 10, lett); 
B_lt_5

ce_eg_5 = ce_from_B(B_lt_5, x.pos = which(colnames(B_lt_5) == "e"), y.pos = which(colnames(B_lt_5) == "g"))
ce_eg_5

##6

DAG6 = amat(DAG)

DAG6[5,1] = 0
DAG6[5,7] = 0
DAG6[1,8] = 0
DAG6[8,7] = 0
DAG6[5,8] = 1
DAG6[8,6] = 1
DAG6[6,3] = 0
DAG6[3,7] = 1
DAG6[6,9] = 1

bn_dag_6 = empty.graph(nodes = lett)
amat(bn_dag_6) = DAG6

## Ripeto gli steps precedenti, ossia:

fitted_network_6 <- bn.fit(bn_dag_6, data = df)

B6 = matrix(0,length(fitted_network_6),length(fitted_network_6))
colnames(B6) = rownames(B6) = lett

for(j in 1:q){
  pa_j = which(!is.na(match(lett, names(fitted_network_6[[j]]$coefficients[-1]))))
  B6[pa_j,j] = fitted_network_6[[j]]$coefficients[-1]
}

B6

adj.mat = amat(fitted_network_6)
B_lt_6 = make_adj_lower_tri(B6, 10, lett); 
B_lt_6

ce_eg_6 = ce_from_B(B_lt_6, x.pos = which(colnames(B_lt_6) == "e"), y.pos = which(colnames(B_lt_6) == "g"))
ce_eg_6

##7

DAG7 = amat(DAG)

DAG7[5,1] = 0
DAG7[5,7] = 0
DAG7[1,8] = 0
DAG7[8,7] = 0
DAG7[5,8] = 1
DAG7[8,6] = 0
DAG7[6,3] = 0
DAG7[3,7] = 1
DAG7[6,9] = 1

bn_dag_7 = empty.graph(nodes = lett)
amat(bn_dag_7) = DAG7

## Ripeto gli steps precedenti, ossia:

fitted_network_7 <- bn.fit(bn_dag_7, data = df)

B7 = matrix(0,length(fitted_network_7),length(fitted_network_7))
colnames(B7) = rownames(B7) = lett

for(j in 1:q){
  pa_j = which(!is.na(match(lett, names(fitted_network_7[[j]]$coefficients[-1]))))
  B7[pa_j,j] = fitted_network_7[[j]]$coefficients[-1]
}

B7

adj.mat = amat(fitted_network_7)
B_lt_7 = make_adj_lower_tri(B7, 10, lett); 
B_lt_7

ce_eg_7 = ce_from_B(B_lt_7, x.pos = which(colnames(B_lt_7) == "e"), y.pos = which(colnames(B_lt_7) == "g"))
ce_eg_7


mean(-0.5319276, -0.5319276, -0.02729919, 0.235926)






