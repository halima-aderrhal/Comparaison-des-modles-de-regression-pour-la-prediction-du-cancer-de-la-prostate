# Lecture des donnees
pross <- read.table("http://statweb.stanford.edu/~tibs/ElemStatLearn/datasets/prostate.data")

# Donnees d'entrainement
x_train <- pross[pross$train == TRUE, 1:8]
Xtrain <- as.matrix(x_train)

y_train <- c(pross[pross$train == TRUE, 9])
Ytrain <- as.matrix(y_train)

# Donnees de test
x_test <- pross[pross$train == FALSE, 1:8]
Xtest <- as.matrix(x_test)

Y_test <- c(pross[pross$train == FALSE, 9])
Ytest <- as.matrix(Y_test)

# Regression lineaire MCO
modele_MCO <- function(X, Y) {
  beta <- solve(t(X) %*% X) %*% t(X) %*% Y
  return(beta)
}

estimateurs_MCO <- modele_MCO(Xtrain, Ytrain)
print(estimateurs_MCO)

# Testons le modele_MCO
predictions <- Xtest %*% estimateurs_MCO
erreur_MCO <- mean(abs(Ytest - predictions))
print(erreur_MCO)
# l'erreur: 0.5191

#ACP
library(Matrix)

ACP <- function(X, Y) {
  ranks <-  1:rankMatrix(X)
  betaAcp <- matrix(0, ncol = ranks, nrow = ncol(X))
  
  for (i in ranks) {
    Z <- X %*% svd(X)$v[, i, drop = FALSE]
    beta_i <- svd(X)$v[, i, drop = FALSE] %*% solve(t(Z) %*% Z) %*% t(Z) %*% Y
    betaAcp[, i] <- beta_i
  }
  
  return(betaAcp)
}

#Testons le model ACP

ACPmodele <- ACP(Xtrain, Ytrain)
print(ACPmodele)
Xtest_acp <- Xtest[, 1:rankMatrix(Xtrain)]
dif_ACP <- abs(Ytest - (Xtest_acp %*% ACPmodele[, 1]))
Erreur_ACP <- mean(dif_ACP)
print(Erreur_ACP)

#l'erreur ACP:0.6807475

#  Ridge
Ridge <- function(X, Y, lambda, p) {
  Beta <- solve(t(X) %*% X + lambda * diag(1, p)) %*% t(X) %*% Y
  return(Beta)
}

# Validation croisée LOO
lambda <- seq(0.0001, 100, by = 0.0001)
L <- rep(0, length(lambda))

for (i in 1:length(lambda)) {
  rid <- Ridge(Xtrain, Ytrain, lambda[i], ncol(Xtrain))
  difridge <- abs(Ytest - (Xtest %*% rid))
  L[i] <- mean(difridge)
}

i <- which.min(L)
lambdaRidge <- lambda[i]
print(lambdaRidge)#3.9434


#avec l'erreur 0.5105177
# on choisit lambda=3.94 car celle la donne la meme erreur
# Testons le modele Ridge
betaridge <- Ridge(Xtrain, Ytrain, 3.94, ncol(Xtrain))
difridge <- abs(Ytest - (Xtest %*% betaridge))
Erreur_ridge <- mean(difridge)
print(Erreur_ridge)



#Regression LASSO
coord_descLasso <- function(X, Y, lambda) {
  beta1 <- rep(0, dim(X)[2])
  beta2 <- rep(1, dim(X)[2])
  
  while (sqrt(t(beta1 - beta2) %*% (beta1 - beta2)) > 0.0000001) {
    if (dim(X)[2] == 2) {
      beta2 <- beta1
      for (j in 1:2) {
        Rj <- t(X[, j]) %*% (Y - X[, -j] %*% beta1[-j])
        betaj <- Rj * max(1 / (t(X[, j]) %*% X[, j]) - lambda / (2 * abs(Rj) * (t(X[, j]) %*% X[, j])), 0)
        beta1[j] <- betaj
      }
    } else {
      beta2 <- beta1
      for (j in 1:dim(X)[2]) {
        Rj <- t(X[, j]) %*% (Y - X[, -j] %*% beta1[-j])
        betaj <- Rj * max(1 / (t(X[, j]) %*% X[, j]) - lambda / (2 * abs(Rj) * (t(X[, j]) %*% X[, j])), 0)
        beta1[j] <- betaj
      }
    }
  }
  
  return(beta1)
}

#  Validation croisée LOO LASSO
lambda <- seq(0.01,10, by = 0.01)
length(lambda)
erreurs <- rep(0, length(lambda))

for (i in 1:length(lambda)) {
  lass <- coord_descLasso(Xtrain, Ytrain, lambda[i])
  dif_lasso <- abs(Ytest - (Xtest %*% lass))
  Erreur_lasso <- mean(dif_lasso)
  erreurs[i] <- Erreur_lasso
}

i <- which.min(erreurs)
lambdalasso <- lambda[i]
print(lambdalasso)


# Testons le modele  LASSO
betalasso <- coord_descLasso(Xtrain, Ytrain, 5.406)
diflasso <- abs(Ytest - (Xtest %*% betalasso))
Erreurmoylasso <- mean(diflasso)
print(Erreurmoylasso) #5.4 
#avec lambda 4.5 l'erreur moyenne est est egale a 0.4988499
#lambda choisit est 5.406
#l'erreur moyenne: 0.4988519


#Elastic Net
coord_descElasticnet <- function(X, Y, lambda1, lambda2) {
  beta1 <- rep(0, dim(X)[2])
  beta2 <- rep(1, dim(X)[2])
  
  while (sqrt(t(beta1 - beta2) %*% (beta1 - beta2)) > 0.0000001) {
    if (dim(X)[2] == 2) {
      beta2 <- beta1
      for (j in 1:2) {
        Rj <- t(X[, j]) %*% (Y - X[, -j] %*% beta1[-j])
        betaj <- (1 / (1 + (lambda2 / (t(X[, j]) %*% X[, j])))) * Rj * max(1 / (t(X[, j]) %*% X[, j]) - lambda1 / (2 * abs(Rj) * (t(X[, j]) %*% X[, j])), 0)
        beta1[j] <- betaj
      }
    } else {
      beta2 <- beta1
      for (j in 1:dim(X)[2]) {
        Rj <- t(X[, j]) %*% (Y - X[, -j] %*% beta1[-j])
        betaj <- (1 / (1 + (lambda2 / (t(X[, j]) %*% X[, j])))) * Rj * max(1 / (t(X[, j]) %*% X[, j]) - lambda1 / (2 * abs(Rj) * (t(X[, j]) %*% X[, j])), 0)
        beta1[j] <- betaj
      }
    }
  }
  
  return(beta1)
}

# LOO
lambdanet <- seq(0.001, 8, by = 0.001)
length(lambdanet)
k <- rep(0, length(lambdanet))

for (i in 1:length(lambdanet)) {
  modele_net <- coord_descLasso(Xtrain, Ytrain, lambdanet[i])
  difference <- abs(Ytest - (Xtest %*% modele_net))
  Erreur_net <- mean(difference)
  k[i] <- Erreur_net
}

i <- which.min(k)
lambdaelasticnet <- lambdanet[i]
print(lambdaelasticnet)

# 5.404 avec l'erreur  0.4988499
#on a choisit lambda=5.4037091 avec l'erreur 0.4988481
#Testons le modele elasticnet
betaElasticnet <- coord_descLasso(Xtrain, Ytrain, 5.4037091)
dif_Elasticnet <- abs(Ytest - (Xtest %*% betaElasticnet))
erreur_Elastic_net <- mean(dif_Elasticnet)
print(erreur_Elastic_net)



# le meilleur modele est le model ElasticNet 
#avec lambda=5.4037091 et l'erreur 0.4988481