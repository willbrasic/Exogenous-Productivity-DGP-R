###############################################################################
### William Brasic ############################################################
###############################################################################
###############################################################################


###############################################################################
### Clear Plot Pane, Environment, and Console #################################
###############################################################################
graphics.off()
rm(list=ls())
cat("\014")


###############################################################################
### Set working Directory ## ##################################################
###############################################################################
setwd("C:/Users/wbras/OneDrive/Desktop/UNLV/GA_Project")


###############################################################################
### Packages ##################################################################
###############################################################################
require(Matrix)



################
### Functions ##
################

## function that rearranges a matrix into a vector ##
long <- function(X) {
  X2 <- c()
  for(i in 1:nrow(X)) {
    X2 <- c(X2, X[i,])
  }
  return(X2)
}




##########################
### Defining Parameters ##
##########################

## firms, time periods, and number of simulations ##
N.values <- c(50, 100, 200, 400)
T_ <- 11 # T_ = # of periods + 1 for the lag
R <- 5


## productivity parameters ## 
rho_0 <- 0.2; rho_1 <- 0.7


## price of materials parameters ## 
phi_0 <- 0.55


## price of energy parameters ## 
nu_0 <- 0.55


## production function parameters ##
beta_K <- 0.25; beta_M <- 0.3; beta_E = 0.4


## investment function parameters ##
theta_0 <- 0.5; theta_1 <- 0.4; theta_2 <- 0.1


## output AR1 coefficient ##
tau_0 <- 0.8


## true values of cost-funtion parameters ##
alpha_Y <- 1/(beta_M + beta_E)
alpha_K <- -beta_K/(beta_M + beta_E)
alpha_M <- beta_M/(beta_M + beta_E)


## number of parameters we are mainly interested ##
number_par <- 4


## null matrix for result storage ##
result_store <- array(NA, c(number_par, R, length(N.values)))


## Evaluation Metric Storage for W_hat ##
MB_W_Store <- matrix(NA, R, length(N.values))
RMSE_W_Store <- matrix(NA, R, length(N.values))
MAE_W_Store <- matrix(NA, R, length(N.values))


## IM Rank Store ##
IM_Rank_Store <- matrix(NA, R, length(N.values))


#######################
### Begin Main Loop  ##
#######################

## Start the loop ##
for(N in N.values){
  set.seed(N)
  
  n.value <- which(N.values==N) # the order number of the N.value currently used
  
  ## depreciation rates ##
  de <- sample(c(0.05, 0.075, 0.1, 0.125, 0.15), N, replace = TRUE)
  
  for(j in 1:R) {
    
    ## initial values ##
    W0 <- runif(N, 1, 3)
    ln.Y0 <- runif(N, 2, 5)
    K0 <- runif(N, 10, 200)
    ln.PM0 <- runif(N, 0, 1)
    ln.PE0 <- runif(N, 0, 0.5)
    I0 <- (exp(ln.Y0)^theta_0)*(K0^theta_1)*(exp(theta_2*W0))
    M0 <- ((exp(ln.Y0)/exp(W0))^(1/(beta_M + beta_E)))*((1/(K0^beta_K))^(1/(beta_M+beta_E)))*((exp(ln.PE0)*beta_M)/(exp(ln.PM0)*beta_E))^(beta_E/(beta_M + beta_E))
    E0 <- ((exp(ln.Y0)/exp(W0))^(1/(beta_M + beta_E)))*((1/(K0^beta_K))^(1/(beta_M+beta_E)))*((exp(ln.PM0)*beta_E)/(exp(ln.PE0)*beta_M))^(beta_M/(beta_M + beta_E))
    
    ## errors ##
    eps <- matrix(rnorm(N*T_, 0, sqrt(0.1)), N, T_)
    eta1 <- matrix(rnorm(N*T_, 0, sqrt(0.07)), N, T_)
    eta2 <- matrix(rnorm(N*T_, 0, sqrt(0.07)), N, T_)
    xi   <- matrix(rnorm(N*T_, 0, sqrt(0.04)), N, T_)   
    
    
    # Evolution processes #
    ln.Y <- K <- M <- E <- I <- W <- ln.PE <- ln.PM <- PM_tilde <- C <- C_tilde <- S <- varkappa <- matrix(NA, N, T_)
    
    for(i in 1:T_){
      if(i == 1) {
        W[,i] <- rho_0 + rho_1*W0 + xi[,i]
        ln.PM[,i] <- phi_0*ln.PM0 + eta1[,i]
        ln.PE[,i] <- nu_0*ln.PE0 + eta2[,i]
        ln.Y[,i] <- tau_0*ln.Y0 + eps[,i]
        K[,i] <- I0 + (1 - de)*K0
        I[,i] <- (exp(ln.Y[,i])^theta_0)*(K[,i]^theta_1)*(exp(theta_2*W[,i]))
      }
      else {
        W[,i] <- rho_0 + rho_1*W[,i-1] + xi[,i]
        ln.PM[,i] <- phi_0*ln.PM[,i-1] + eta1[,i]
        ln.PE[,i] <- nu_0*ln.PE[,i-1] + eta2[,i]
        ln.Y[,i] <- tau_0*ln.Y[,i-1] + eps[,i]
        K[,i] <- I[,i-1] + (1 - de)*K[,i-1]
        I[,i] <- (exp(ln.Y[,i])^theta_0)*(K[,i-1]^theta_1)*(exp(theta_2*W[,i-1]))
      }
      M[,i] <- ((exp(ln.Y[,i])/exp(W[,i]))^(1/(beta_M + beta_E)))*((1/(K[,i]^beta_K))^(1/(beta_M+beta_E)))*((exp(ln.PE[,i])*beta_M)/(exp(ln.PM[,i])*beta_E))^(beta_E/(beta_M + beta_E))
      E[,i] <- ((exp(ln.Y[,i])/exp(W[,i]))^(1/(beta_M + beta_E)))*((1/(K[,i]^beta_K))^(1/(beta_M+beta_E)))*((exp(ln.PM[,i])*beta_E)/(exp(ln.PE[,i])*beta_M))^(beta_M/(beta_M + beta_E))
    }
    
    
    for(i in 1:T_){
      C[,i] <- exp(ln.PM[,i])*M[,i] + exp(ln.PE[,i])*E[,i]
      C_tilde[,i] <- C[,i] / exp(ln.PE[,i])
      PM_tilde[,i] <- exp(ln.PM[,i])/exp(ln.PE[,i])
      S[,i] <- (exp(ln.PM[,i])*M[,i]) / C[,i]
      varkappa[,i] <- log(C_tilde[,i]) - S[,i]*log(PM_tilde[,i])
    }
    
    
    
    ## defining alpha_M_hat ##
    alpha_M_hat <- mean(S)
  
    
    
    ########################
    ### Step 2 Estimation ##
    ########################
    
    
    lagged_K <- K[,-T_]; lagged_ln.Y <- ln.Y[,-T_]; lagged_varkappa <- varkappa[,-T_]
    
    
    df <- data.frame(id = rep(1:N, each = (T_-1)), year=rep(1:(T_-1), N), ln.Y=long(ln.Y[,-1]), lagged_ln.Y = long(lagged_ln.Y), 
                     K=long(K[,-1]), lagged_K = long(lagged_K), M=long(M[,-1]), E = long(E[,-1]), 
                     I = long(I[,-1]), W = long(W[,-1]), ln.PM = long(ln.PM[,-1]), ln.PE = long(ln.PE[,-1]), 
                     varkappa = long(varkappa[,-1]), lagged_varkappa = long(lagged_varkappa), 
                     C = long(C[,-1]), C_tilde = long(C_tilde[,-1]), PM_tilde = long(PM_tilde[,-1]),
                     S = long(S[,-1]))
    
    
    ## Define minimization algorithm ##
    OBJ <- function(par) {
      sum((par[1]*log(df$K) + par[2]*df$ln.Y - df$varkappa - par[2]*par[3]
           - par[4]*(par[1]*log(df$lagged_K) + par[2]*df$lagged_ln.Y
                     - df$lagged_varkappa))^2)
      
      # par[1] = alpha_K; par[2] = alpha_Y; par[3] = lambda_0; par[4] = lambda_1
    }
    
    
    ## Find parameters that minimize OBJ ##
    result <- optim(par = c(alpha_K, alpha_Y, rho_0, rho_1) + rnorm(n = 4, mean = 0, sd = 0.1), 
                    fn = OBJ, method = 'BFGS')
    
    
    ## Store Step-2 estimates of alpha_K, alpha_Y, rho_0, rho_1 ##
    result_store[,j,n.value] <- c(result$par)
    
    
    #############################
    ### Productivity Estimates ##
    #############################
    
    ## productivity prediction ##
    W_hat <- (1/result$par[2])*(result$par[1]*log(df$K) + alpha_M_hat*log(df$PM_tilde)
                                                + result$par[2]*df$ln.Y
                                                - log(df$C_tilde)) 
    
    ## storing MB, RMSE, and MAE for W prediction ##
    MB_W_Store[j, n.value] <- (1/(T_-1)*N)*sum(W_hat - df$W)
    RMSE_W_Store[j, n.value] <- sqrt((1/(T_-1)*N)*sum((W_hat - df$W)^2))
    MAE_W_Store[j, n.value] <- (1/(T_-1)*N)*sum(abs(W_hat - df$W))
    
    
    ########################
    ## Information Matrix ##
    ########################
    
    ## defining z as defined in paper ##
    z <- alpha_K*log(df$lagged_K) + alpha_Y*df$lagged_ln.Y + df$lagged_varkappa
    
    ## IM row 1 ##
    IM_1.1 <- alpha_Y^2
    IM_1.2 <- z
    IM_1.3 <- -alpha_Y*log(df$K) + rho_1*log(df$lagged_K) 
    IM_1.4 <- -alpha_Y*(df$ln.Y + rho_0)
    
    
    ## IM row 2 ##
    IM_2.1 <- IM_1.2 
    IM_2.2 <- ((1/alpha_Y)*z)^2
    IM_2.3 <- z*log(df$lagged_K)*((-alpha_Y*log(df$K) + rho_1)/alpha_Y^2) 
    IM_2.4 <- z*(-(df$ln.Y - rho_0)/alpha_Y)
    
    
    ## IM row 3 ##
    IM_3.1 <- IM_1.3
    IM_3.2 <- IM_2.3
    IM_3.3 <- (log(df$K) - log(df$lagged_K)*(rho_1/alpha_Y))^2
    IM_3.4 <- (log(df$K) - log(df$lagged_K)*(rho_1/alpha_Y))*(df$ln.Y - rho_0)
    
    
    ## IM row 4 ##
    IM_4.1 <- IM_1.4
    IM_4.2 <- IM_2.4
    IM_4.3 <- IM_3.4
    IM_4.4 <- (df$ln.Y - rho_0)^2
    
    
    ## creating IM ##
    IM <- matrix(data = c(mean(IM_1.1), mean(IM_1.2), mean(IM_1.3), mean(IM_1.4),
                          mean(IM_2.1), mean(IM_2.2), mean(IM_2.3), mean(IM_2.4),
                          mean(IM_3.1), mean(IM_3.2), mean(IM_3.3), mean(IM_3.4),
                          mean(IM_4.1), mean(IM_4.2), mean(IM_4.3), mean(IM_4.4)),
                 nrow = 4, ncol = 4, byrow = TRUE)
    
    
    ## storing rank of IM ##
    IM_Rank_Store[j,n.value] <- rankMatrix(IM)[1]
    
  }
  
  print(n.value)
  
}
  



###############################################
## Summarize the results (excluding alpha_M) ##
###############################################

true_store <- rbind(rep(1, R))%x%cbind(c(alpha_K,alpha_Y,rho_0,rho_1))  # matrix of the replicated true values 

print('Evaluation Metrics for Parameters of Main Interest')
for (N in N.values) {
  
  n.value <- which(N.values==N) # the order number of the N.value currently used
  
  ## Mean Bias ##
  MB_par <- apply(result_store[,,n.value]-true_store, 1, function(g) { 1/R*sum(g) })
  Mean_MB_W <- apply(MB_W_Store, 2, FUN = mean)
  
  ## RMSE ##
  RMSE_par <- apply(result_store[,,n.value]-true_store,1, function(g) { sqrt(1/R*sum(g^2)) })
  Mean_RMSE_W <- apply(RMSE_W_Store, 2, FUN = mean)
  
  ## MAE ##
  MAE_par <- apply(result_store[,,n.value]-true_store,1, function(g) { 1/R*sum(abs(g)) })
  Mean_MAE_W <- apply(MAE_W_Store, 2, FUN = mean)
  
  
  par_metrics <- cbind(MB_par, RMSE_par, MAE_par)  
  colnames(par_metrics) <- c('MB','RMSE','MAE')
  rownames(par_metrics) <- c('alpha_K','alpha_Y','rho_0','rho_1')
  
  cat('N=',N,'\n')
  print(round(par_metrics, 5))
  cat('\n')
  
  
}

cat('\n')
W_Metrics <- cbind(Mean_MB_W, Mean_RMSE_W, Mean_MAE_W)
colnames(W_Metrics) <- colnames(par_metrics)
rownames(W_Metrics) <- c('N = 50', 'N = 100', 'N = 200', 'N = 400')
print('Evaluation Metrics for W')
print(W_Metrics)
cat('\n')


cat('\n')
print('Insufficient Rank Rate')
for(N in N.values){
  
  n.value <- which(N.values==N)
  
  IM_Rank_Check <- length(which(IM_Rank_Store[, n.value] != 4))
  
  IM_Rank_Error_Rate <- IM_Rank_Check / R
  
  cat('N=',N,'\n')
  print(round(IM_Rank_Error_Rate, 5))
  cat('\n')
  
}



