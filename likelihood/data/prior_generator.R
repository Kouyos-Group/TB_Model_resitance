library(lhs)

nb = 30 #number of line parameter
lhs<-randomLHS(nb, 8) #optimization of the combinations of parameters that we want to infer with the ABC


#infection
inf<-10
sup<-10
lhs[,1]<-10^(sup*lhs[,1]-inf)

#s
a <-0.1 #inf
b <-1 #sup
lhs[,2]<-(b-a)*lhs[,2]+a

#beta
inf<-2
sup<-2
lhs[,3]<-10^(sup*lhs[,3]-inf)

#tau_rate
inf<-3
sup<-3
lhs[,4]<-10^(sup*lhs[,4]-inf)

#I
I<-lhs[,5]*1/(lhs[,5]+lhs[,6]+lhs[,7]+lhs[,8])*100000

#E
E<-lhs[,6]*1/(lhs[,5]+lhs[,6]+lhs[,7]+lhs[,8])*100000

#L
L<-lhs[,7]*1/(lhs[,5]+lhs[,6]+lhs[,7]+lhs[,8])*100000


theta<-data.frame(infection = lhs[,1],
         s = lhs[,2],
         beta = lhs[,3],
         t = lhs[,4],
         I = I,
         E = E,
         L = L
         )


write.table(theta,"prior.txt",sep="\t",row.names=F)

