library(glmnet)
data <-read.csv(file=file.choose(),h=TRUE,sep=";")
grille = seq(0, 1, 0.1)


x = model.matrix(y~., data)[,-1]
y = data$y


##Question 10
lasso.mod = glmnet(x, y, alpha = 1, lambda = grille, intercept = FALSE)
coeff = coef(lasso.mod)


nbind = nrow(x)
ind.train = sample.int(nbind,round(0.7*nbind))
x.train = x[ind.train,]
y.train = y[ind.train]
x.test = x[-ind.train,]
y.test = y[-ind.train]



#############################
####Code sans fonction R#####
#############################

#on standardise la fonction pur pouvoir faire le Lasso
standardize <- function(x){
  x=scale(x,T,F)
  v=sqrt(apply(x**2, 2, function(x) mean(x)))
  x=sweep(x, 2, v, "/")
  return(x)
}

#calcul du betachap pour le lasso
beta <- function(y, x, betai, lambda){
  err = 1
  temp= betai
  while (err>10^(-3)) {
    
    betai = temp
    for (j in (1:length(x[1,]))) {
      n = length(x[,1])
      vec = rep(0,n )
      
      for (i in (1:length(x[,1]))){
        a = x[i,j]
        if (j>1){
        b1 = sum(x[i,1:j-1]*temp[1:j-1])}
        else{b1=0}
        
        if (j<length(x[1,])){
        b2 = sum(x[i,(j+1):length(x[1,])]*temp[(j+1):length(x[1,])])
        }else{b2=0}
        
        vec[i] = a*(y[i]-b1-b2)
      }
      
      rj = (1/n)*sum(vec)
      temp[j] =  rj*max(1-lambda/abs(rj),0)
      }
    err = sqrt(sum((temp - betai)^2))
  }
  return(temp)
}

betaopt = beta(y,standardize(x),rep(1,length(x[1,])), 0.5)
#notre fonction nous sort le vecteur (2.78 0.31 0 00 1.35 0 0 0)
#ce qui est relativement proche de la fonction glmnet qui sort
# (2.67 0.54 0 0 1.5 0 0 0)


#fonction de validation croisée
val_cr<-function(y,x, lambda){
  test_lamb = rep(0, length(lambda))
  for (j in (1:length(lambda))){
    test = rep(0, length(y))
    beta_i = 0
    ychap_i = 0
    for (i in (1:length(x[,1]))){
      beta_i = beta(y[-i],standardize(x[-i,]), rep(1,length(x[1,])), lambda[j])
      ychap_i = sum(x[i,]*beta_i)
      test[i] = (y[i]-ychap_i)^2
    }
    res = sum(test)*(1/length(y))
    test_lamb[j] = res
    
  }
  pos = which.min(test_lamb)
  return(lambda[pos])
  
}


#resultat donnée pour ind.train = (6  4  1  8  3  7 10 19 18  2  9 13 11  5)
lambda_cal = val_cr(y.train, x.train, seq(0, 4, 0.1)) #lambda optimal = 0.7
lambda_cal # = 0.2


betachap = beta(y.train, standardize(x.train), rep(1,length(x.train[1,])) , lambda_cal)
betachap #= (3.2918158 0.9408017 0.0000000 0.0000000 1.6536833 0.0000000 0.0000000 0.0000000)

beta0 = c(3,1.5,0,0,2,0,0,0)
#Question 13
lambdaa = seq(0,4,0.01)
grille1 = lambdaa[order(lambdaa, decreasing = TRUE)]

fonc1 = rep(0, length(grille1))
a = c(standardize(x.test) %*% beta0)
for (i in (1:length(grille1))){
  bet = beta(y.train, standardize(x.train), rep(1,length(x.train[1,])) , grille1[i])
  b = c(standardize(x.test) %*% bet)
  fonc1[i] = sum((a-b)^2)
}
grille1[which.min(fonc1)]#la lambda qui minimise la Q13 est égale 0.23

plot(grille1, fonc1 ,type="b",cex=0.2, col = "mediumseagreen", xlab = "lambda" , ylab = "norme", main = "Norme de l'erreur fonction de Lambda")
points(grille1[which.min(fonc1)], fonc1[which.min(fonc1)], col = "firebrick" , pch = 10 , cex = 2)


a_opt1= c(standardize(x.test) %*% beta0)
b_opt1= c(standardize(x.test) %*% beta(y.train, standardize(x.train), rep(1,length(x.train[1,])) , grille1[which.min(fonc1)]))
sum((a_opt1-b_opt1)^2) #1.30

#les résultats trouvés dépendent beaucoup du ind.train pris.

#############################
####Code avec fonction R#####
#############################
#Question 10 : validation croisée
cv.out = cv.glmnet(x.train,y.train,alpha=1, intercept=FALSE,standardize = FALSE)
plot(cv.out, main="lasso")

lambda_opt_lasso = cv.out$lambda.min #donne la valeur de lambda optimal par validation croisée
lambda_opt_lasso #0.81


lambda_opt_cv= glmnet(x.train, y.train, alpha = 1, lambda=lambda_opt_lasso, intercept = FALSE, standardize = FALSE)
coef(lambda_opt_cv)

#Nous souhaitons choisir correctement le paramètre lambda. On sépare l'échantillon d'apprentissage en deux parties : 
#- Un premier échantillon comprenant environ 90% de l'échantillon d'apprentissage à partir duquel nous allons calculer l'estimateur pour différentes valeurs de lambda.
#- Un second échantillon contenant le reste de l'échantillon d'apprentissage qui nous servira à évaluer les performances des différents modèles obtenus.

#La fonction cv.glmnet() répète cette opération 10 fois en sélectionnant le premier souséchantillon
#de façon à ce que chaque individu soit sélectionné une (et une seule) fois. À la
#fin de l'opération, nous obtenons une prédiction de y pour tous les individus de l'échantillon
#et pour tous les lambda de la grille et nous sélectionnons la valeur de lambda pour
#laquelle l'erreur quadratique moyenne est minimale.


#question 12
lambdaa = seq(0,4,0.001)
grille1 = lambdaa[order(lambdaa, decreasing = TRUE)]
lasso.mod1 = glmnet(x.train, y.train, alpha = 1, lambda = grille1, intercept = FALSE, standardize =FALSE)
coeff1 = coef(lasso.mod1)
norm_0 = apply(as.matrix(coeff1!=0), 2, sum)
plot(grille1, norm_0 , xlab = "lambda" ,col = "mediumseagreen", ylab = "coeff non nul", main = "Nombre coefficients non nuls en fonction de Lambda")

plot(lasso.mod1, xvar="lambda") #graph des valeurs des coefficients en fonction du log(lambda)

#question 13

beta0 = c(3,1.5,0,0,2,0,0,0)
coefff = coeff1[2:9,]
fonc = rep(0, length(grille1))
a = c(x.test %*% beta0)
for (i in (1:length(grille1))){
  b = c(x.test %*% coefff[,i])
  fonc[i] = sum((a-b)^2)
}
grille1[which.min(fonc)]

plot(grille1, fonc ,type="b",cex=0.2, col = "mediumseagreen", xlab = "lambda" , ylab = "norme", main = "Norme de l'erreur fonction de Lambda")
points(grille1[which.min(fonc)], fonc[which.min(fonc)], col = "firebrick" , pch = 10 , cex = 2)

#pour cette question la lambda optimal est en fait 0.143 alors que 
#le lambda optimal obtenu par validation croisée était légèrement différent (0.81)

lambda_opt_13= glmnet(x.train, y.train, alpha = 1, lambda=grille1[which.min(fonc)], intercept = FALSE, standardize = FALSE)
coef(lambda_opt_13) #coefficient pour la lambda qui minimise l'erreur

a_opt= c(x.test %*% beta0)
b_opt= c(x.test %*% coef(lambda_opt_13)[2:9])
sum((a_opt-b_opt)^2) #0.96

#l'erreur trouvé avec la lambda de la validation croisée est : 
a_cv= c(x.test %*% beta0)
b_cv = c(x.test %*% coef(lambda_opt_cv)[2:9])
sum((a_cv-b_cv)^2) #20.37


#question bonus
lambda_en = seq(0,10,0.01)
grille2 = lambda_en[order(lambda_en, decreasing = TRUE)]
E_N.mod = glmnet(x.train, y.train, alpha = 0.75, lambda=grille2, intercept = FALSE)
coeff_en = coef(E_N.mod)
norm_0_en = apply(as.matrix(coeff_en!=0), 2, sum)
plot(grille2, norm_0_en , xlab = "lambda" ,col = "mediumseagreen", ylab = "coeff non nul", main = "Nombre coefficients non nuls en fonction de lambda")
#ici nous avons tracé à alpha fixé le nombre de variable selectionnée par le modèle en fonction de lambda.
plot(E_N.mod, xvar="lambda")

#validation croisée pour le choix du lambda à alpha fixé pour l'elastic net
cv_en = cv.glmnet(x.train,y.train,alpha=0.75, intercept=FALSE)
plot(cv_en, main="Elastice Net")

lambda_opt_en = cv_en$lambda.min #donne la valeur de lambda optimal par validation croisée pour l'elastic net
lambda_opt_en #0.76

E_N.opt= glmnet(x.train, y.train, alpha = 0.75, lambda=lambda_opt_en, intercept = FALSE)
coef(E_N.opt)

#l'erreur trouvé (q13) avec la lambda de la validation croisée est : 
a_en = c(x.test %*% beta0)
b_en = c(x.test %*% coef(E_N.opt)[2:9])
sum((a_en-b_en)^2) #16.63
