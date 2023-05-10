############################################################################################
###Note: The R codes are written based on MATLAB codes prepared by Garry Koop and vailable #
###at https://sites.google.com/site/garykoop/home/computer-code-2.                         #
############################################################################################
# Latest revision 17-04-2023

rm(list=ls()) 
# Load the packages required for the programme
library("largeList")
library("MASS")

# load the fuctions defined in TVP_funs.r
source("TVP_funs_R1.r")

# Load the data
Y <- read.csv("data.csv",header=T)
Y_lab <- as.Date(strptime(Y[,1], "%d-%b-%y"))
Y_lab <-format(Y_lab,"%b-%y")
Y <- as.matrix(Y[,-1]) 

#Number of observations and dimension of X and Y
t <- nrow(Y)
M <- ncol(Y)

#Number of factors & lags:
tau <- 40 # tau is the size of the training sample
p <- 2 # p is number of lags in the VAR part
numa <- M*(M-1)/2 # Number of lower triangular elements of A_t (other than 0's and 1's)

#Generate lagged Y matrix. This will be part of the X matrix
ylag <- cbind(makelag(Y,1), makelag(Y,2))[-c(1:(tau+p)),]
K <- M + p*(M^2) # K is the number of elements in the state vector
Z <- matrix(0, ncol=K, nrow=(t-tau-p)*M) #Create Z_t matrix

for (i in 1:(t-tau-p))
{
    ztemp <- diag(1,M)
    for (j in 1:p)
        {       
        xtemp = ylag[i,(((j-1)*M+1):(j*M))]
        xtemp = t(kronecker(diag(1,M),xtemp))
        ztemp = cbind(ztemp, xtemp)  
        }
    Z[(((i-1)*M+1):(i*M)),] <- ztemp
}

y <- t(Y[((tau+p+1):t),])
y_lab <-Y_lab[((tau+p+1):t)] 
t <- ncol(y)

####################### Preliminaries ######################
#Set some Gibbs - related preliminaries
nrep = 30000  # Number of replications
nburn = 10000 # Number of burn-in-draws
bn <- 3000

####################### PRIORS ######################
# the ts_prior() function is given seperately
priors <- ts_prior(Y,tau,M,p)
A_OLS <- priors$A_OLS
B_OLS <- priors$B_OLS
VA_OLS <- priors$VA_OLS
VB_OLS <- priors$VB_OLS
sigma_OLS <- priors$sigma_OLS
sigma_prmean <- sigma_OLS
sigma_prvar <- 4*diag(M)
S_prmean <- mm<-matrix(list(), M-1, 1)
S_prvar <- mm<-matrix(0, M-1, 1)

####################### INITIALIZE MATRICES  ######################
consQ <- 0.0001
consS <- 0.0001
consH <- 0.01
consW <- 0.0001

# Set some hyperparameters 
k_Q <- 0.01
k_S <- 0.1
k_W <- 0.01

sizeW <- M # Size of matrix W
sizeS <- 1:M #Size of matrix S

B_0_prmean <- B_OLS
B_0_prvar <-  4*VB_OLS

# A_0 ~ N(A_OLS, 4Var(A_OLS))
A_0_prmean <- A_OLS
A_0_prvar <- 4*VA_OLS

# Q ~ IW(k2_Q*size(subsample)*Var(B_OLS),size(subsample))
Q_prmean <- ((k_Q)^2)*tau*VB_OLS
Q_prvar <- tau

# W ~ IG(k2_W*(1+dimension(W))*I_n,(1+dimension(W)))
W_prmean <- ((k_W)^2)*matrix(1,nrow=M,ncol=1)
W_prvar = 2

ind <-1
for (ii in 2:M){
    #S is block diagonal 
S_prmean[[ii-1]] <-((k_S)^2)*(1 + sizeS[(ii-1)])*VA_OLS[((ii-1)+(ii-3)*(ii-2)/2):ind,((ii-1)+
                                                            (ii-3)*(ii-2)/2):ind]
S_prvar[ii-1] = 1 + sizeS[ii-1]
    ind=ind+ii}

Ht <- kronecker(matrix(1,nrow=t,ncol=1),consH*diag(M)) 
Htchol <- kronecker(matrix(1,nrow=t,ncol=1),sqrt(consH)*diag(M))# Cholesky of Htdraw defined above

Qdraw <- consQ*diag(K) # Initialize Qdraw, a draw from the covariance matrix Q
Sdraw <- consS*diag(numa)# Initialize Sdraw, a draw from the covariance matrix S

Sblockdraw <- matrix(list(), M-1, 1)
ijc = 1
for (jj in 2:M){
    Sblockdraw[[jj-1]] <-Sdraw[((jj-1)+(jj-3)*(jj-2)/2):ijc,((jj-1)+(jj-3)*(jj-2)/2):ijc]
    ijc <- ijc + jj
}

Wdraw <- consW*matrix(1,nrow=M,ncol=1) # Initialize Wdraw, a draw from the covariance matrix W
Btdraw <-matrix(0,nrow=K,ncol=t) # Initialize Btdraw, a draw of the mean VAR coefficients, B(t)
Atdraw <-matrix(0,nrow=numa,ncol=t) # Initialize Atdraw, a draw of the non 0 or 1 elements of A(t)
Sigtdraw <-matrix(0,nrow=t,ncol=M) # Initialize Sigtdraw, a draw of the log-diagonal of SIGMA(t)
# Matrix of the exponent of Sigtdraws (SIGMA(t))
sigt <- kronecker(matrix(0,nrow=t, ncol=1),0.01*diag(M)) 
statedraw <- 5*matrix(1,nrow=t,ncol=M) # initialize the draw of the indicator variable 
                               # (of 7-component mixture of Normals approximation)
Zs <-kronecker(matrix(1,nrow=t,ncol=1),diag(M))

#  Storage matrices for posteriors and stuff
Bt_postmean <-matrix(0,nrow=K,ncol=t)#  regression coefficients B(t)
At_postmean <-matrix(0,nrow=numa,ncol=t)#  lower triangular matrix A(t)
Sigt_postmean <-matrix(0,nrow=t,ncol=M) #  diagonal std matrix SIGMA(t)
Qmean <-matrix(0,nrow=K,ncol=K) #  covariance matrix Q of B(t)
Smean <-matrix(0,nrow=numa,ncol=numa)#  covariance matrix S of A(t)
Wmean <-matrix(0,nrow=M,ncol=1)#  covariance matrix W of SIGMA(t)

sigmean <-matrix(0,nrow=t,ncol=M)#  mean of the diagonal of the VAR covariance matrix
cormean <-matrix(0,nrow=t,ncol=numa) #  mean of the off-diagonal elements of the VAR cov matrix
sig2mo <-matrix(0,nrow=t,ncol=M)# squares of the diagonal of the VAR covariance matrix
cor2mo <-matrix(0,nrow=t,ncol=numa)#  squares of the off-diagonal elements of the VAR cov matrix

###### IMPULSE RESPONSES Storage matrices######
istore <-1
if (istore == 1) { 
    nhor <- 21; # Impulse response horizon
    bigj <- matrix(0, nrow=M, ncol=M*p)
    bigj[1:M,1:M] = diag(M)
    }
imprespt <- matrix(0,nrow=M*t,ncol=M*nhor)

################################# END OF PRELIMINARIES #################################
####################################START SAMPLING#####################################
# Depending on the system requirements, this could take a while.  Once executed these codes will save 1 the impulse response coefficients of each simulations after the burn-in period as' .llo' files

# Progress Bar
total <- nrep + nburn
pb <- txtProgressBar(min = 0, max = total, style = 3)
for ( irep  in 1:(nrep + nburn)){   
Sys.sleep(0.1)
    ################  STEP I: Sample B from p(B|y,A,Sigma,V)#################

# Accept draw
 Btdraw <- carter_kohn(y,Z,Ht,Qdraw,K,M,t,B_0_prmean,B_0_prvar)

#=====| Draw Q, the covariance of B(t) (from iWishart)
# Take the SSE in the state equation of B(t)

Btemp <- t(Btdraw[,2:t]) - t(Btdraw[,1:t-1])

sse_2 <- matrix(0,nrow=K,ncol=K)
for (i in 1:t-1){
    sse_2 <- sse_2+as.matrix(Btemp[1,])%*%t(Btemp[1,])
}

#  ...and subsequently draw Q, the covariance matrix of B(t)
Qinv <- solve(sse_2+Q_prmean,tol=1e-21)

Qinvdraw <- wish(Qinv,t-1+Q_prvar)
Qdraw <- solve(Qinvdraw) # this is a draw from Q
    
    ######################  STEP II: Draw A(t) from p(At|y,B,Sigma,V) ######################

# Drwas alpha

yhat <- matrix(0, nrow=M,ncol=t)
for (i in 1: t){
    yhat[,i] <- y[,i]-Z[(((i-1)*M+1):(i*M)),]%*%Btdraw[,i]
}

Zc <- -t(yhat)
sigma2temp <- exp(Sigtdraw)

 # Draw each block of A(t)
Atdraw <- carter_kohn(t(as.matrix(yhat[2,])),as.matrix(Zc[,1]),as.matrix(sigma2temp[,2]),
            Sblockdraw[[1]],sizeS[1],1,t,A_0_prmean[1,],A_0_prvar[1,1])
       

ind <- 3
for (ii in 3:M){
Atblockdraw <- carter_kohn(t(as.matrix(yhat[ii,])),Zc[,1:(ii-1)],as.matrix(sigma2temp[,ii]),
            Sblockdraw[[ii-1]],sizeS[ii-1],1,t,A_0_prmean[(((ii-1)+(ii-3)*(ii-2)/2):ind),],
                A_0_prvar[(((ii-1)+(ii-3)*(ii-2)/2):ind),(((ii-1)+(ii-3)*(ii-2)/2):ind)]);
        Atdraw <- rbind(Atdraw,Atblockdraw) # Atdraw is the final matrix of draws of A(t)
        ind = ind + ii;
    }
####### Draw S, the covariance of A(t) (from iWishart)
    #### Take the SSE in the state equation of A(t)

Attemp <- t(Atdraw[,2:t]) - t(Atdraw[,(1:(t-1))])
sse_2 <- matrix(0,nrow=numa,ncol=numa)

for (i in 1:(t-1)){
sse_2+(Attemp[i,])%*%t(Attemp[i,])
}

ijc <-1
for(jj in 2:M){
Sinv <- solve(sse_2[(((jj-1)+(jj-3)*(jj-2)/2):ijc),(((jj-1)+(jj-3)*(jj-2)/2):ijc)] + 
                   S_prmean[[jj-1]])

Sinvblockdraw <- wish(h=Sinv,n=t-1+S_prvar[jj-1])

Sblockdraw[[jj-1]] = solve(Sinvblockdraw)
ijc <- ijc+jj
}
    
    ########################### STEP III: Draw diagonal VAR covariance matrix log-SIGMA(t)#############
# First create capAt, the lower-triangular matrix A(t) with ones on the
# main diagonal. 

capAt <- diag(M)
aatemp <- as.matrix(Atdraw[,1])
    ic=1
    for (j in 2:M){
        capAt[j,(1:(j-1))] <- t(aatemp[(ic :(ic+j-2)),1])
        ic <- ic+j-1
    }

for (i in 2:t){
    capatemp <- diag(M)
    aatemp <- as.matrix(Atdraw[,i])
     ic=1
    for (j in 2:M){
        capatemp[j,(1:(j-1))] <- t(aatemp[(ic :(ic+j-2)),1])
        ic <- ic+j-1
    }
    capAt<-rbind(capAt,capatemp)
}

ytemps = capAt[(1:M),]%*%yhat[,1]
y2 <- (ytemps^2)


for (i in 2:t){
        ytemps <-  capAt[((i-1)*M+1):(i*M),]%*%yhat[,i]
        y2 <- cbind(y2, (ytemps^2))
  }

yss <-t(log(y2+ 10^(-6)))

for (j in 1:M){
Sigtdraw[,j] <- SVRW2(as.matrix(yss[,j]),Sigtdraw[,j],Wdraw[j,],sigma_prmean[j],
                      sigma_prvar[j,j],1)  
    }

sigt <- exp(0.5*Sigtdraw)
e2 <- Sigtdraw[-1,] - Sigtdraw[-(nrow(Sigtdraw)),]
W1 <- W_prvar + t - p - 1
W2 <- W_prmean + rowSums(t(e2^2))
Winvdraw <- matrix(rgamma(n=M,shape=W1/2,scale=2/W2))
Wdraw <- 1/Winvdraw

#  Create the VAR covariance matrix H(t). It holds that:
#  A(t) x H(t) x A(t)' = SIGMA(t) x SIGMA(t) '

Ht <- matrix(0,nrow=M*t,M)
Htsd <- matrix(0,nrow=M*t,M)

for (i in 1:t){
    inva <- solve(capAt[(((i-1)*M+1):(i*M)),])
    stem <- diag(sigt[i,])
    Hsd <- inva%*%stem
    Hdraw <- Hsd%*%t(Hsd)
    Ht[(((i-1)*M+1):(i*M)),] <- Hdraw #H(t)
    Htsd[(((i-1)*M+1):(i*M)),] <- Hsd  #Cholesky of H(t)
}

##########SAVE AFTER-BURN-IN DRAWS AND IMPULSE RESPONSES ##############
if (irep > nburn){               
Bt_postmean <- Bt_postmean + Btdraw # regression coefficients B(t)
At_postmean <- At_postmean + Atdraw # lower triangular matrix A(t)
Sigt_postmean <- Sigt_postmean + Sigtdraw  # diagonal std matrix SIGMA(t)
Qmean <- Qmean + Qdraw# covariance matrix Q of B(t)

ikc <-1
for(kk in 2:M){
    Sdraw[(((kk-1)+(kk-3)*(kk-2)/2):ikc),(((kk-1)+(kk-3)*(kk-2)/2):ikc)] <-Sblockdraw[[kk-1]]
    ikc <- ikc + kk
}
            
Smean <-  Smean + Sdraw #covariance matrix S of A(t)
Wmean <- Wmean + Wdraw #covariance matrix W of SIGMA(t)
#  Get time-varying correlations and variances

# for the first period
stemp6 <- matrix(0,nrow=M,ncol=1)
stemp8 <- corrvc(Ht[1:M,])
stemp7a <- matrix(0,nrow=M,ncol=0)
ic <-1
for(j in 1:M){
    if (j>1){
       stemp7a <- c(stemp7a,t(stemp8[j,(1:ic)])) 
       ic <- ic+1             
    } 
    stemp6[j,1] <- sqrt(Ht[j,j])
     
    }

stemp5 <- t(stemp6)
stemp7 <- t(stemp7a)                    
                   
for (i in 2:t) {
stemp6 <- matrix(0,nrow=M,ncol=1)
stemp8 <- corrvc(Ht[((i-1)*M+1):(i*M),])
stemp7a <- matrix(0,nrow=M,ncol=0)
ic <-1
for(j in 1:M){
    if (j>1){
       stemp7a <- c(stemp7a,t(stemp8[j,(1:ic)])) 
       ic <- ic+1             
    } 
    stemp6[j,1] <- sqrt(Ht[((i-1)*M+j),j])
     
    }
 stemp5 <- rbind(stemp5,t(stemp6))
 stemp7 <- rbind(stemp7,t(stemp7a))  
}                   

sigmean <- sigmean + stemp5 # diagonal of the VAR covariance matrix
cormean <- cormean + stemp7 # off-diagonal elements of the VAR cov matrix
sig2mo <- sig2mo + stemp5^2
cor2mo <- cor2mo + stemp7^2

# Impulse response analysis. 
if (istore==1){
biga <- matrix(0,nrow=M*p,ncol=M*p) # M= number of variables, P is nu
for(j in 1:(p-1)){
    biga[(j*M+1):(M*(j+1)), (M*(j-1)+1):(j*M)] <- diag(M)
}
for (i in 1:t){ #Get impulses recurssively for each time period
     bbtemp <- Btdraw[((M+1):K),i]
     splace <- 0
     for (ii in 1:p){
         for (iii in 1:M){
             biga[iii,(((ii-1)*M+1):(ii*M))] <- bbtemp[((splace+1):(splace+M))]
            splace <- splace + M 
         }
     }
#Identification code: 
# St dev matrix for structural VAR
Hsd <- Htsd[(((i-1)*M+1):(i*M)),(1:M)]  
Hsd
#  First shock is the Cholesky of the VAR covariance
diagonal = diag(diag(Hsd))
Hsd = solve(diagonal)%*%Hsd   # Unit initial shock

# Now get impulse responses for 1 through nhor future periods
impresp <- matrix(nrow=M,ncol=M*nhor)
impresp[1:M,1:M] <-Hsd

# First shock is the Cholesky of the VAR covariance
bigai <-biga
for (j in 1:(nhor-1)){
    impresp[,(j*M+1):((j+1)*M)] <- bigj%*%bigai%*%t(bigj)%*%Hsd
    bigai <- bigai%*%biga
}
if(i==1){ # please check
imprespt[1:M,] <-impresp
    } else{
    imprespt[((i-1)*M+1):(i*M),] <-impresp
}
}
#### the draws are saved as different files using large list package #####  
    
block <- floor((irep - nburn - 1)/bn) + 1
name1 <- paste0("S", irep - nburn)
ar12 <- list(A = imprespt)
names(ar12)[1] <- name1
compress <- (irep - nburn - bn*(block - 1) == 1)
file_name <- paste0("irfmcmc", block, ".llo")
saveList(object = ar12, file = file_name, append = !compress, compress = compress)

} #END the impulse response calculation section 
} #END saving after burn-in results 

# update progress bar
   setTxtProgressBar(pb, irep)
} #END main Gibbs loop (for irep = 1:nrep+nburn)
close(pb) 
###############################End of Sampling###########################

########################################################################
#### Codes to replicate the figures in the paper########################
########################################################################

# Please make sure that  the .llo files named (irfmcmc1.llo,...,irfmcmc10.llo) that we generated is saved in your home directory

# The following code use fuctions from 'largeList' packages and user deined fuctions from 'TVP_funs_R1.r'
# Please uncomment the lines below to load these files

#library(largeList)
#source (source("TVP_funs_R1.r")

# Load the packages required for replicating the figures
library("xts")
library("ggplot2")
library("scales")
library("plotly")
library("rgl")
library("scales")
library("patchwork")
library("reshape2")


# replication of 'figure 2', 'figure 3' and 'figure 2A' in appendix (This may take some time to execute)
results_figures(var_i_s = 1, var_i_r = 4, np = 4) 
    
# replication of 'figure 4', 'figure 5' and 'figure 3A' in appendix
results_figures(var_i_s = 2, var_i_r = 4, np = 4)   

# replication of 'figure 6', 'figure 7' and 'figure 4A' in appendix
results_figures(var_i_s = 3, var_i_r = 4, np = 4)   

# replication of 'figure 8', 'figure 9' and 'figure 5A' in appendix
results_figures(var_i_s = 4, var_i_r = 1, np = 4)    
    
### The following section reproduces figure 1 and fiure 1A
# load data on FSI 
FSIdata <- read.csv("FSI_data.csv", header=TRUE)

# Define the sequence of dates with month granularity
Period <- seq(as.Date("2001/1/1"), as.Date("2022/3/1"), "month")
   
# Create the data frame for plotting
fig1_df <- data.frame(Period, FSI.VEW = FSIdata$FSI.VEW, FSI.PCA = FSIdata$FSI.PCA)
fig1_df <- melt(fig1_df, id="Period")
options(repr.plot.width = 9, repr.plot.height =5)

# Create the plot with modified x-axis labeling
fig1 <- ggplot(fig1_df, aes(x=Period, y=value)) + 
  geom_line(aes(linetype=variable), linewidth=1) +
  scale_linetype_manual(values = c(3,1)) +
  scale_x_date(labels = date_format("%b-%y"), 
               breaks = date_breaks("29 months")) +  # adjust breaks to "1 year"
  ggtitle(paste("Figure 1: Financial stress index for India for the period 2001 to 2022")) +
  xlab("") +
  ylab("") +
  theme_classic() +
  theme(axis.text = element_text(size = 14, color = "black"), 
        axis.title = element_text(size = 14), 
        plot.title = element_text(size = 16),
        legend.text = element_text(size = 14, color = "black"),
        legend.position = "bottom") +
  guides(linetype = guide_legend(title = NULL))

# Reproduces the figre 1 in the paper
fig1
    
# For the plots in figure 1A (appendix)
create_graph2 (df=FSIdata, 1) 
create_graph2 (df=FSIdata, 2)
create_graph2 (df=FSIdata, 3) 
create_graph2 (df=FSIdata, 4) 
create_graph2 (df=FSIdata, 5)     
    