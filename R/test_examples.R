x_mat <- function(J,nper=1){
  t <- J+1
  X <- kronecker(rep(1,J),diag(t))
  #XJ <- kronecker(diag(J),rep(1,t))
  int <- c()
  for(i in 1:J){
    int <- c(int,rep(c(rep(0,t-(i)),rep(1,i)),nper))
  }
  if(nper>1){
    XJ <- kronecker(diag(nper),XJ)
    X <- kronecker(rep(1,nper),X)
  }
  #X <- cbind(X,XJ)
  
  X <- cbind(X, int)
  return(X)
}

nJ <- 6
nT <- nJ+1
M <- 20

df <- expand.grid(t=1:nT,J = 1:nJ)
df <- df[rep(1:nrow(df),each=M),]

X <- x_mat(nJ)
X <- X[rep(1:nrow(X),each=M),]
#some random parameters for the example, assume cluster means are zero
beta <- c(rep(0,nT),0.5)

data <- as.data.frame(X)
colnames(data) <- c(paste0("t",1:(nT)),"int")

# you can change the "indicator" variance below as needed
Sout <- gen_re_mat(df,
                   dims = list(1,2),
                   funs = list(
                     list("exp_power",c(1)),
                     list("indicator",c(0.01))
                   ), parallel = TRUE)

ss <- 100
d <- sample(c(rep(1,ss),rep(0,nrow(df)-ss)),nrow(df))
idx_in <- which(d==1)

dat1 <- create_data(formula(paste0("~Linear(",paste0("t",1:nT,collapse=", "),", int)")),
                    family=gaussian(),
                    data=data,
                    theta=list(
                      list(beta)
                    ),
                    Z=Sout[[1]],
                    D=Sout[[2]],
                    C=matrix(c(rep(0,nT),1)),
                    var_par = 1)



d1 <- grad_robust(idx_in,
                  C=list(dat1[[1]]),
                  X_list=list(dat1[[2]]),
                  sig_list = list(dat1[[3]]),
                  tol=0,
                  trace=T)

## plot results for d
Xq <- data
Xq <- cbind(Xq,incl=0)
Xq[sort(d1),ncol(Xq)] <- 1
Xq <- cbind(Xq,cl=rep(1:nJ,each=nT*M))
Xq <- cbind(Xq,t=rep(rep(1:nT,each=M),nJ))
Xqa <- aggregate(Xq[,c('incl','int')],list(Xq[,'cl'],Xq[,'t']),sum)
Xqa$int <- Xqa$int/M
Xqa$lab <- paste0("N=",Xqa$incl)

pcl <- ggplot(data=Xqa,aes(x=Group.2,y=Group.1,fill=incl,color=factor(int)))+
  geom_tile()+
  geom_label(aes(label=incl))+
  scale_fill_viridis_c(name="N per block",option="A")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="T",y="Cluster")+
  scale_x_continuous(breaks = 1:nT)+
  scale_y_continuous(breaks = 1:nJ)+
  scale_color_discrete(name="Intervention",labels=c("No","Yes"))

pcl100 <- optim_cl_des(gaussian(link="identity"),
                       Sout = Sout,
                       idx_in = idx_in,
                       data=data,
                       nT=nT,nJ=nJ,M=M)
pcl100 <- pcl100+ggtitle(TeX(r'($\eta =0.25$, $r=1$)'))
pcl100