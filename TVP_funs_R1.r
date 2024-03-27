#########################################################################################################
##Function 1 lags
makelag<- function(y,nlag) 
{
  #arg y is a matrix
  # nlag number of lages
  nr=nrow(y)
  nc=ncol(y)
  y1 <- y[1:(nr-nlag),]
  ylag <- rbind(matrix(0,nrow=nlag,ncol=nc),y1)
  return(ylag)
} 
#########################################################################################################
# Wishart distribution
wish<-function(h,n){
    m <- nrow(h)
    A <- t(chol(h))%*% matrix( rnorm(m*n,mean=0,sd=1), m, n) 
    A <- A%*%t(A)
}

#########################################################################################################
# Function 2  Priors

ts_prior <- function(rawdat,tau,p,plag){
yt <- t(rawdat[((plag+1):(tau+plag)),])
m  <- p + (plag*(p^2)) # number of state vectors
Zt <-matrix(ncol=m) # Zt
for (i in (plag+1):(tau+plag))
{
    ztemp <- diag(1,p)
    for (j in 1:plag)
        {       
        xlag = rawdat[(i-j),(1:p)]
        xtemp = matrix(0,nrow=p, ncol=p*p)
        for(jj in 1:p){
         xtemp[jj,(((jj-1)*p)+1):(jj*p)] = xlag  
        }
    ztemp <- cbind(ztemp,xtemp)
    }
    Zt <-rbind(Zt,ztemp)
}

Zt <-Zt[-1,]

   
# see the main document for details
vbar <- matrix(0,nrow=m,ncol=m) 
xhy <- matrix(0,ncol=1,nrow=m) 

for (i in 1:tau){
    zhat1 <- Zt[(((i-1)*p+1):(i*p)),]
    vbar <- vbar + (t(zhat1)%*%zhat1)
    xhy <- xhy + t(zhat1)%*%yt[,i]
    }

vbar <- solve(vbar)
aols  <- vbar%*%xhy

sse2 <- matrix(0,nrow=p,ncol=p)

for (i in 1:tau){
    zhat1 <- Zt[((i-1)*p+1):(i*p),]
    sse2 = sse2 + (yt[,i] - (zhat1%*%aols))%*%t(yt[,i] - zhat1%*%aols)
}

hbar = sse2/tau # element wise right division

vbar <- matrix(0,ncol=m, nrow=m)

for (i in 1 : tau){
    zhat1 = Zt[((i-1)*p+1):(i*p),]
    vbar = vbar + t(zhat1)%*%solve(hbar)%*%zhat1
    }

library(Matrix) # for chol()

vbar <- solve(vbar)
achol <- t(chol(hbar))

ssig  <- matrix(0,nrow=p,ncol=p)


for( i in 1:p){
    ssig[i,i] = achol[i,i] 
    for (j in 1:p){
        achol[j,i] = achol[j,i]/ssig[i,i]
 }
}

achol  <- solve(achol)
numa <-  p*(p-1)/2
a0 <- matrix(0,nrow=numa,ncol=1)

ic = 1

for (i in 2:p){
    for (j in 1:(i-1)){
        a0[ic,1] <- achol[i,j]
        ic <- ic+1
    }
}

ssig1 <- matrix(0,nrow=p,ncol=1)

for (i in 1:p){
    ssig1[i,1] <-log(ssig[i,i]^2)

}

hbar1 <- solve(tau*hbar)

hdraw <- matrix(0,nrow=p,ncol=p)
a02mo <- matrix(0,nrow=numa,ncol=numa)
a0mean <-matrix(0,nrow=numa,ncol=1)

for(irep in 1:4000){
 hdraw <-  wish(hbar1,tau) 
  hdraw <- solve(hdraw)
achol <- t(chol(hdraw))
ssig <- matrix(0,nrow=p,ncol=p)
for (i in 1:p){
    ssig[i,i] <- achol[i,i]
    for (j in 1:p){
        achol[j,i] <- achol[j,i]/ssig[i,i]
    }
}
achol <- solve(achol)
a0draw <- matrix(0,nrow=numa,ncol=1)
ic=1
for (i in 2:p){
    for (j in 1:(i-1))
        {
        a0draw[ic,1] <- achol[i,j]
        ic <- ic+1
        
    }
}
a02mo <- a02mo + a0draw%*%t(a0draw)
    a0mean <- a0mean + a0draw
}

a02mo  <-  a02mo/4000
a0mean <- a0mean/4000;
a02mo <- a02mo - a0mean%*%t(a0mean)
results <- list(   B_OLS = aols,VB_OLS =vbar,A_OLS=a0,sigma_OLS =ssig1,VA_OLS =a02mo)
return(results)
}
#########################################################################################################
# Carter and Kohn (1994), On Gibbs sampling for state space models.
carter_kohn <- function(y,Z,Ht,Qt,m,p,t,B0,V0){

# Kalman Filter
bp <- B0

Vp <- V0
bt <- matrix(0,nrow=t, ncol=m)
Vt <- matrix(0, nrow=m^2,ncol=t)
log_lik <- 0
R <- matrix(0,nrow=p,ncol=p)
H <- matrix(0,nrow=p,ncol=m)

for (i in 1:t){

R[1:p,] <- Ht[((i-1)*p+1):(i*p),]
H[1:p,] <- Z[((i-1)*p+1):(i*p),]
cfe <- y[,i] - H%*%bp   #conditional forecast error
f<- H%*%Vp%*%t(H)+R #variance of the conditional forecast error
inv_f<- solve(f)
#log_lik<-log(det(f))+t(cfe)%*%inv_f%*%cfe
btt <- bp+Vp%*%t(H)%*%inv_f%*%cfe
 Vtt <- Vp -Vp%*%t(H)%*%inv_f%*%H%*%Vp
if(i<t){
    bp <- btt
    Vp <- Vtt+Qt
    
}
bt[i,] <- btt
Vt[,i] <- c(Vtt)
}

# draw Sdraw(T|T) ~ N(S(T|T),P(T|T))
bdraw <-matrix(0,nrow=t,ncol=m)

bdraw[t,] <- mvrnorm(mu=btt,Sigma=Vtt,n=1)

# Backward recurssions
for (i in 1:(t-1)){
    bf <- (bdraw[(t-i+1),])
    btt <- (bt[(t-i),])
    Vtt <- matrix(Vt[,t-i],nrow=m)
    f <- Vtt + Qt
    inv_f <- solve(f)
    cfe <- bf - btt
    bmean <- btt + Vtt%*%inv_f%*%cfe
    bvar <- Vtt - Vtt%*%inv_f%*%Vtt
    bdraw[t-i,] = mvrnorm(mu=bmean,Sigma=bvar,n=1); #bmean' + randn(1,m)*chol(bvar);
}
bdraw <- t(bdraw)
    return(bdraw)
    }
#########################################################################################################
# Carter and Kohn (1994), On Gibbs sampling for state space models.
carter_kohn2 <- function(y,Z,Ht,Qt,m,p,t,B0,V0,kdraw){

# Kalman Filter

bp <- B0

Vp <- V0
bt <- matrix(0,nrow=t, ncol=m)
Vt <- matrix(0, nrow=m^2,ncol=t)
log_lik <- 0

# R <- matrix(0,nrow=p,ncol=p)
 H <- matrix(0,nrow=p,ncol=m)

for (i in 1:t){

# R[1:p,] <- Ht[((i-1)*p+1):(i*p),]
H[1:p,] <- Z[((i-1)*p+1):(i*p),]
# F <- diag(m)    
cfe <- y[,i] - H%*%bp   #conditional forecast error
f<- H%*%Vp%*%t(H)+Ht[,,i] #variance of the conditional forecast error
# inv_f<- solve(f)
inv_f<- t(H)/f    
#log_lik<-log(det(f))+t(cfe)%*%inv_f%*%cfe
btt <- bp+Vp%*%inv_f%*%cfe
 Vtt <- Vp -Vp%*%inv_f%*%H%*%Vp
if(i<t){
    bp <- btt
    Vp <- Vtt+kdraw[i,]%*%Qt
    
}
bt[i,] <- t(btt)# check the dim
Vt[,i] <- c(Vtt)
}

# draw Sdraw(T|T) ~ N(S(T|T),P(T|T))
bdraw <-matrix(0,nrow=t,ncol=m)

bdraw[t,] <- mvrnorm(mu=btt,Sigma=Vtt,n=1)

# Backward recurssions
for (i in 1:(t-1)){
    bf <- (bdraw[(t-i+1),])
    btt <- (bt[(t-i),])
    Vtt <- matrix(Vt[,t-i],nrow=m)
    f <-  Vtt+kdraw[t-i,]%*%Qt
#     inv_f <- solve(f)
    inv_f <- Vtt/f
    cfe <- bf - btt
    bmean <- btt + inv_f%*%cfe
    bvar <- Vtt - inv_f%*%Vtt
    bdraw[t-i,] = mvrnorm(mu=bmean,Sigma=bvar,n=1); #bmean' + randn(1,m)*chol(bvar);
}
bdraw <- t(bdraw)
    return(bdraw)
   }

## function Repeats columns, and rows in a matrix
rep.row<-function(x,n){
   matrix(rep(x,each=n),nrow=n)
}

rep.col<-function(x,n){
   matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

SVRW2 <- function(Ystar,h,sig,sigma_prmean,sigma_prvar,TVP_Sigma){
    

T <- length(h)

# normal mixture
pi  <- c(0.0073, .10556, .00002, .04395, .34001, .24566, .2575)
mi <-c(-10.12999, -3.97281, -8.56686, 2.77786, .61942, 1.79518, -1.08819) - 1.2704 
# means already adjusted!! %%
sigi <- c(5.79596, 2.61369, 5.17950, .16735, .64009, .34023, 1.26261)
sqrtsigi <- sqrt(sigi)


# Sample S from a 7-point distrete distribution
temprand <- as.matrix(runif(T))




q <-rep.row(pi,T)*dnorm(rep.col(Ystar,7),(rep.col(h,7)+rep.row(mi,T)),rep.row(sqrtsigi,T)) 
# element wise multiplication
q <- q/rep.col(rowSums(q),7)

S <- 7-as.matrix(rowSums(1*(rep.col(temprand,7)<t(apply(q,1,cumsum)))))+1

vart <- array(0,c(1,1,T))
yss1 <- matrix(0,nrow=T,ncol=1)

for (i in 1:T){
    imix <-S[i]
    vart[1,1,i] <-sigi[imix]
    yss1[i,1] <- Ystar[i,1]-mi[imix]
}

h <- carter_kohn2(t(yss1),matrix(1,nrow=T,ncol=1),
                  vart,sig,m=1,p=1,t=T,sigma_prmean,sigma_prvar,kdraw =(TVP_Sigma*matrix(1,nrow=T,ncol=1)))


return(h)
}
########################################################################################################
#Computes a correlation matrix from a
#               variance-covariance matrix.
corrvc <- function(vc){
if (is.complex(vc) == TRUE){       
    print("ERROR: Not implemented for complex arguments.")
} else {
std <- sqrt(as.matrix(diag(vc)))    
}

     output = vc/(std%*%t(std))

    return(output)
}

#########################################################################################################
# Function to reproduce the figures in the paper
results_figures <- function(var_i_s,var_i_r,np ){# var_i_r is the reponse variables due to shocks in var_i_s, and np is number of variablbe in the VAR
 ##############################    
 # the function given below read simulated impluse respsonse functions saved in .llo files and rearranges 
read_data_matrix <- function(file_name, var_i_s, var_i_r, np) { 
    read_llo <- readList(file_name)
    col_nm <- length(read_llo[[1]][1,]) # length of the selected row
    row_nm <- length(read_llo[[1]][,1]) # length of the selected column
    N_sim <- length(read_llo) # number of simulations
    select_col <- seq(var_i_s, col_nm, by = np)
    select_row <- seq(var_i_r, row_nm, by = np)
    irfmcmc_by_time <- list()
    counter <- 1
    for (j in select_row){ # j is actualy var_i_r at different time periods
       
        irf_sim <- matrix(0, nrow = N_sim, ncol = length(select_col))
            for (i in 1:N_sim) {
                irf_sim[i,] <- read_llo[[i]][j, select_col]
                              }
        irfmcmc_by_time[[counter]] <- irf_sim 
        counter <- counter + 1
                        }
    return(irfmcmc_by_time)
}
#################################
# the function retruns quantile(qntl) impluse respsonses of var_i_r due to shocks in var_i_s
calculate_imp_res_quantile <- function( np, var_i_r, var_i_s,  qntl){ #qntl =0.5 is median, .16 is 16th and .84 is 84th quntiles
    
file_names <- paste0("irfmcmc", 1:10, ".llo")

# Define a function to read a data matrix from a file
read_data <- function(file_name) {
  read_data_matrix(file_name, var_i_s, var_i_r, np)
}

# Read the data matrices from the files
dd_list <- lapply(file_names, read_data)
t <-length(dd_list[[1]])
imp_res <- matrix(0, nrow = t, ncol = 21)
for (i in 1:length(dd_list[[1]])) {
  imp_res_arrange <- do.call(rbind, lapply(dd_list, "[[", i))
  imp_res[i,] <- apply(imp_res_arrange, 2, quantile, probs = qntl)
} 
    return(imp_res)
}
#####################################
# function to calculate row sums (num_col =1+ nth period ahead)
calculate_cols <- function(df, num_cols) {
  cols <- rowSums(df[, c(2:(1+num_cols))])
 # colnames(cols) <- c("L", "M", "U")
  return(cols)
}
#######################################
# Define a function to create figure 1 below
create_plot <- function(df) {
  
  ggplot(df, aes(x = period, y = M)) +
    geom_line() +
    geom_ribbon(aes(ymin = L, ymax = U), alpha = 0.3) +
    scale_x_date(labels = date_format("%b-%Y"), breaks = date_breaks("3 year")) +
     xlab("") + ylab("") +
    theme(panel.background = element_rect(fill = NA),
          legend.position = "bottom", 
          axis.line.y = element_line(color = "black"), 
          axis.text = element_text(size = 14, color="black")) +
    geom_hline(yintercept = 0)
}
########################################
#Prepare the data from .llo files    
    
dd16 <- data.frame(calculate_imp_res_quantile(var_i_s = var_i_s, var_i_r = var_i_r, np = np, qntl = 0.16))
dd50 <- data.frame(calculate_imp_res_quantile(var_i_s = var_i_s, var_i_r = var_i_r, np = np, qntl = 0.5))
dd84 <- data.frame(calculate_imp_res_quantile(var_i_s = var_i_s, var_i_r = var_i_r, np = np, qntl = 0.84))

Time <- seq(as.Date("2004/7/1"), as.Date("2022/2/1"), "month")
y_lab <-format(Time,"%b-%y")
dd16 <- cbind(y_lab,dd16)
dd50 <- cbind(y_lab,dd50)
dd84 <- cbind(y_lab,dd84)

colnames(dd16)<-c("period",0:20)
colnames(dd50)<-c("period",0:20)
colnames(dd84)<-c("period",0:20)

#### Defines R and S for plot title
lookup_s <- c("oil supply", "oil demand", "risk-related", "FSI")
lookup_r <- c("oil prices", "RXOI", "VIX", "FSI")

S <- lookup_s[var_i_s]
R <- lookup_r[var_i_r]

df4<- data.frame(Time,calculate_cols(dd16,5),calculate_cols(dd50,5),calculate_cols(dd84,5))
colnames(df4) <- c("period","L","M","U")

df8<- data.frame(Time,calculate_cols(dd16,9),calculate_cols(dd50,9),calculate_cols(dd84,9))
colnames(df8) <- c("period","L","M","U")

df12<- data.frame(Time,calculate_cols(dd16,13),calculate_cols(dd50,13),calculate_cols(dd84,13))
colnames(df12) <- c("period","L","M","U")
 # Figure 1
figure1_T <- paste("Accumulated impulse responses of", R, "to", S, "shocks")


figure1a <- create_plot(df4)+ggtitle("(a)") + theme(plot.title = element_text(size = 16, hjust = 0.5))
figure1b <- create_plot(df8)+ggtitle("(b)") + theme(plot.title = element_text(size = 16, hjust = 0.5))
figure1c <- create_plot(df12) +ggtitle("(c)") + theme(plot.title = element_text(size = 16, hjust = 0.5))  
figure1 <-(figure1a/ figure1b/ figure1c + plot_layout(ncol = 1,
                widths = c(9.75, 9.75,9.75), heights = c(3.2, 3.2, 3.2)))

#Figure2
dd_fig2 <-dd50[dd50$period  %in%  c("Aug-04","Feb-09","Jul-16","Dec-21"),] #Selecting specific values based on Month and year
dd_rname <-as.vector(dd_fig2[,1])
dd_fig2 <- dd_fig2[,-1]
new_dd_fig2 <- apply(dd_fig2, 1, cumsum)
colnames(new_dd_fig2) <- dd_rname
Periods <- 0:20
new_dd_fig2 <- data.frame(Periods,new_dd_fig2)

d <- melt(new_dd_fig2, id="Periods")
{options(repr.plot.width = 9.75, repr.plot.height =4.5)

#png(filename="figure3.png", width=800, height=450)
figure2 <- ggplot(d, aes(x=Periods, y=value,)) + 
geom_line(aes(linetype=variable), linewidth=1) +
scale_linetype_manual(values = c(3,1,4,2))+
ggtitle(paste("Accumulated impulse responses of", R, "to", S, " shocks at specific time periods")) +
  xlab("Periods") +
  ylab("Accumulated impulse coefficients")+ theme_classic()+ 
    theme(axis.text = element_text(size = 14, color = "black"), 
          axis.title = element_text(size = 16), 
          plot.title = element_text(size = 16),
         legend.text = element_text(size = 14, color = "black"),
         legend.position = "bottom")+
    guides(linetype = guide_legend(title = NULL))}
#dev.off()

# Figure3
figure3_T <- paste("Median impulse responses of", R, "due to", S, "shocks")
IRF <- as.matrix(dd50[,-1])

axy <- list(
  ticketmode = 'array',
  ticktext = c("Jul 2004", "Aug 2008", "Oct 2012", "Dec 2016", "Feb 2021"),
  tickvals = c(1,50,100,150,200),
  range = c(-5,212),title = ""
)

axx <- list(
  nticks = 5,
  range = c(0,21),title = ""
)

# axz <- list(  nticks = 5,  range = c(0,.9))

figure3 <- plot_ly(z = ~IRF)
figure3 <- figure3 %>%add_surface(colorbar = list(title = "IRF")) %>% layout(scene = list(xaxis=axx,yaxis=axy)) 

cat(figure1_T, "\n\n")
par(mfrow= c(3,1))
print(figure1a)
print(figure1b)
print(figure1c)
print(figure2) 
cat(figure3_T, "\n\n")
figure3
}
######################################################################
create_graph2 <- function(df,col_n){

Period <- seq(as.Date("2001/1/1"), as.Date("2022/3/1"), "month")
lookup_v <- c("FSI_MM", "FSI_DM", "FSI_EM", "FSI_B","FSI_FX")
V <- lookup_v[col_n]
fig2_df <- data.frame(Period,value=FSIdata[[V]])

time_series_plot <- ggplot(fig2_df, aes(x = Period, y = value))+
  geom_line() +
  geom_ribbon(aes(ymin = 0, ymax = value), fill = "grey", alpha = 0.2)+
theme_classic() +
  theme(axis.line.x = element_blank()) +
  labs(x = "", y = "", title = V)+
  theme(axis.line.y = element_line(color = "black"), 
        axis.text = element_text(size = 12, color="black"),
        plot.title = element_text(hjust = 0.5))+
  geom_rect(aes(xmin = min(Period) - 0.5, xmax = max(Period) + 
                0.5, ymin = -Inf, ymax = Inf),
            fill = NA, color = "black", linewidth = 1)+
  geom_segment(aes(x = min(Period) - 0.5, y = 0, xend = max(Period) + 0.5, yend = 0),
               color = "black", linewidth = .5) # Add vertical lines
date_list <- c("2002-01-01", "2004-01-01", "2006-01-01", "2008-01-01", "2010-01-01", 
               "2012-01-01", "2014-01-01", "2016-01-01", "2018-01-01", "2020-01-01",
               "2022-01-01")

vline_list <- lapply(date_list, function(x) {
  geom_vline(xintercept = as.numeric(as.Date(x)), color = "grey", linewidth = 0.75)
})

time_series_plot <- time_series_plot + vline_list

# Rotate x-axis tick labels
time_series_plot <- time_series_plot + 
  scale_x_date(labels = date_format("%b-%y"), breaks = date_breaks("29 months"))

#### For box plot

melted_data <- melt(fig2_df, id.vars = "Period")



# Create boxplot legend
boxplot_legend <- ggplot(melted_data, aes(x = variable, y = value)) +
  geom_boxplot(width = 0.5, fill = "white", colour = "black", outlier.size = 1) +
  theme_void() +
  labs(x = NULL, y = "Value", title = "")

# Combine plots using patchwork
(time_series_plot + boxplot_legend + plot_layout(ncol = 2,
                widths = c(11, 1), heights = c(.75, 1)))
}
