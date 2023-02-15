###Load packages
library(TiPS) #Tuto to use TiPS: https://cran.r-project.org/web/packages/TiPS/vignettes/TiPS.html


###Write Reactions
#####
#INFECTED STRAIN (0: not resistant, 1:resistant)
###first X---: rifampicin resitant
###second -X--: isoniazid resistant
###third --X-: aminoglycosides resistant
##fourth ---X: belaquilin resistant
#paste0("s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"), each=8))

#TREATMENT
#t1=>1st line treatment (rifampicin and isoniazid) : effective and low toxicity. Developement of resistance 
#t2=>old 2nd line (aminoglycosides) : not very effective and toxic
#t3=>new 2nd line (bedaquiline) : effective and toxic. Developpement of resistance but less as compare to the 1st line
#t4=> 3rd line : I don't know really 
#paste0("t",c("1","2","3","4"))



#####
reactions <- c(
  #####birth
  paste0("0 [ pi * N ] -> S"),
         
  
  ####infection
  paste0("S [", 
         " S * ", 
         "alphas", c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8), " * ",
         c(paste0("I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8)), 
           paste0("T","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8))),
         " ] -> ", paste0("E","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8))
         ),
  
  ####Progression from early latent to late latent
  paste0(paste0("E","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8)), " [ ", 
         paste0("E","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8)), " * ",
         "epsilon ",
         "] -> ", paste0("L","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8))
         ),
  
  ####Progression from latent to active disease
  paste0(c(paste0("E","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8)),
           paste0("L","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8))), " [ ", 
         rep(c("delta1","delta2"), each=8), " * ",
         c(paste0("E","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8)),
           paste0("L","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8))),
         " ] -> ",paste0("I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8))
  ),
  
  ####disease progression from active disease to treatment
  paste0(paste0("I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8))," [ ",
         "psi * (1 - sigma) * ",
         paste0("I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8)),
         " ] -> ", paste0("T","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8)), "t1"
         ), #detected but not test to know with which strain
  
  paste0(paste0("I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8))," [ ",
         "psi * sigma * (1 - omega) * ",
         "iota", rep(c("1","2","3","4"),each=16),
         " * ",paste0("I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8)),
         " ] -> ",paste0("T","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8), 
                         "t", rep(c("1","2", "3", "4"), each=16))),
  #detected, tested but wrong result
  
  paste0(paste0("I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8))," [ ",
         "psi * sigma * omega *",
         " beta", "s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8), 
         "t", rep(c("1","2", "3", "4"), each=16),
         " * ", paste0("I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8)),
         " ] -> ",paste0("T","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8), 
                         "t", rep(c("1","2", "3", "4"), each=16))), 
  #detected, tested, good result
  
  
  ####re-infection from Late latent and recovered to early latent stage
  paste0(paste0("L","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8))," [", 
         paste0("L","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8)),"* ", 
         "LambdaL", rep(c("0","1"),each=16), rep(c("0","1"),each=32),rep(c("0","1"),each=64),rep(c("0","1"),each=128), " * ",
         c(paste0("I","s",rep(c("0","1"),each=16),rep(c("0","1"),each=32),rep(c("0","1"),each=64), rep(c("0","1"),each=128)),
           paste0("T","s",rep(c("0","1"), each=16),rep(c("0","1"),each=32),rep(c("0","1"),each=64), rep(c("0","1"),each=128))), 
         " ] -> ", "Es", rep(c("0","1"),each=16), rep(c("0","1"),each=32),rep(c("0","1"),each=64),rep(c("0","1"),each=128)
           ), #for late latent
  
  
  paste0(paste0("R","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8))," [", 
         paste0("R","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8)),"* ", 
         "LambdaR", rep(c("0","1"),each=16), rep(c("0","1"),each=32),rep(c("0","1"),each=64),rep(c("0","1"),each=128), " * ",
         c(paste0("I","s",rep(c("0","1"),each=16),rep(c("0","1"),each=32),rep(c("0","1"),each=64), rep(c("0","1"),each=128)),
           paste0("T","s",rep(c("0","1"), each=16),rep(c("0","1"),each=32),rep(c("0","1"),each=64), rep(c("0","1"),each=128))), 
         " ] -> ", "Is", rep(c("0","1"),each=16), rep(c("0","1"),each=32),rep(c("0","1"),each=64),rep(c("0","1"),each=128)
  ), #for recovered
  
  
  ####successful treatment
  paste0(paste0("T","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8), 
                "t", rep(c("1","2", "3", "4"), each=16))," [ ",
         "tau","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8), 
         "t", rep(c("1","2", "3", "4"), each=16),
         " * ",paste0("T","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8), 
                      "t", rep(c("1","2", "3", "4"), each=16)), 
         " ] -> R","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8)),
  
  ####cured by itself
  paste0(c(paste0("I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8)), 
           paste0("T","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8), 
                  "t", rep(c("1","2", "3", "4"), each=16))),
         " [ omicron * ",
         c(paste0("I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8)), 
           paste0("T","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8), 
                  "t", rep(c("1","2", "3", "4"), each=16))), 
         " ] -> R","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),rep(c("0","1"),each=8)),
  
  ####treatment fail, strains become resistant
  #acquisition rifampicin resistance
  paste0("T","s0",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), "t1 [ ",
         "T","s0",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t1", 
         " * rho1000",
         "] -> Ts1",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), "t1" ),
  
  #acquisistion isoniazid
  paste0("T","s",c("0","1"),"0",rep(c("0","1"),each=2),rep(c("0","1"),each=4), "t1", " [ ",
         "T","s",c("0","1"),"0",rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t1 ",
         "* rho0100",
         "] -> T","s",c("0","1"),"1",rep(c("0","1"),each=2),rep(c("0","1"),each=4), "t1"),
  
  
  #acquisistion short course old 2nd line
  paste0("T","s",c("0","1"),rep(c("0","1"),each=2),"0", rep(c("0","1"),each=4),"t2", " [ ",
         "T","s",c("0","1"),rep(c("0","1"),each=2),"0",rep(c("0","1"),each=4),"t2 ",
         "* rho0010 ] -> Ts",c("0","1"),rep(c("0","1"),each=2),"1",rep(c("0","1"),each=4),"t2"),
 
   #acquisistion short course new 2nd line
  paste0("T","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"0", "t3", " [ ",
         "T","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"0","t3 ",
         "* rho0001 ] -> Ts",c("0","1"),rep(c("0","1"),each=2),"1",rep(c("0","1"),each=4),"t3"),
  
  #Treatment shift
  #not treated with the first line but sensible 
  paste0(c("Ts0010", "Ts0001", "Ts0011"), rep(c("t2","t3", "t4"), each=3), " [ ",
         c("Ts0010", "Ts0001", "Ts0011"), rep(c("t2","t3", "t4"), each=3), " * ", 
         "ita", 
         " ] -> ", c("Ts0010", "Ts0001", "Ts0011"), "t1"), 
  
  #resistant to the the 1st line but not to the 2nd and not treated with the second line old one
  paste0(c("Ts1000", "Ts0100", "Ts1100"), rep(c("t1","t4"), each=3), " [ ",
         c("Ts1000", "Ts0100", "Ts1100"), rep(c("t1","t4"), each=3), " * ", 
         "ita_old", 
         " ] -> ", c("Ts1000", "Ts0100", "Ts1100"), "t2"), 
  
  
  #resistant to the the 1st line but not to the 2nd and not treated with the second line old one
  paste0(c("Ts10", "Ts01", "Ts11"), rep(c("0","1"),each=3),"0", rep(c("t1","t4"), each=6), " [ ",
         c("Ts10", "Ts01", "Ts11"), rep(c("0","1"),each=3),"0", rep(c("t1","t4"), each=6), " * ", 
         "ita_new", 
         " ] -> ", c("Ts10", "Ts01", "Ts11"), rep(c("0","1"),each=3),"0t3"), 
  
  #resistant to the the 1st line and the 2nd one
  paste0(c("Ts10", "Ts01", "Ts11"), rep(c("10","01", "11"),each=3), rep(c("t1","t2", "t3"), each=9), " [ ",
         c("Ts10", "Ts01", "Ts11"), rep(c("10","01", "11"),each=3), rep(c("t1","t2", "t3"), each=9), " * ", 
         "ita", 
         " ] -> ", c("Ts10", "Ts01", "Ts11"), rep(c("10","01", "11"),each=3),"t4"), 
  
  ####death due to TB
  paste0(c(paste0("I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8)),
           paste0("T","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8), 
                  "t", rep(c("1","2", "3", "4"), each=16))), " [ ",
         "theta * ",
         c(paste0("I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8)),
           paste0("T","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8), 
                  "t", rep(c("1","2", "3", "4"), each=16))),
         " ] -> S"),
  
  ####Relapse from cured to active stage
  paste0(paste0("R","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),rep(c("0","1"),each=8))," [ ",
         "phi * ",
         paste0("R","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),rep(c("0","1"),each=8)), 
         " ] -> ",paste0("I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),rep(c("0","1"),each=8))),
  
  ####natural death
  paste0(c("S",
           paste0("L","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8)),
           paste0("E","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8)),
           paste0("I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8)), 
           paste0("T","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), rep(c("0","1"),each=8), 
                  "t", rep(c("1","2", "3", "4"), each=16)),
           paste0("R","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),rep(c("0","1"),each=8))),
         " [  mu * N ] -> S")
)




#####



###Build the simulator that will allow to run multiple trajectories

simplification<-c(

paste0("N=S+Es0000+Es1000+Es0100+Es1100+Es0010+Es1010+Es0110+Es1110+Es0001+Es1001+Es0101+Es1101+Es0011+Es1011+Es0111+Es1111+",
"Ls0000+Ls1000+Ls0100+Ls1100+Ls0010+Ls1010+Ls0110+Ls1110+Ls0001+Ls1001+Ls0101+Ls1101+Ls0011+Ls1011+Ls0111+Ls1111+",
"Is0000+Is1000+Is0100+Is1100+Is0010+Is1010+Is0110+Is1110+Is0001+Is1001+Is0101+Is1101+Is0011+Is1011+Is0111+Is1111+",
"Ts0000t1+Ts1000t1+Ts0100t1+Ts1100t1+Ts0010t1+Ts1010t1+Ts0110t1+Ts1110t1+Ts0001t1+Ts1001t1+Ts0101t1+Ts1101t1+Ts0011t1+Ts1011t1+Ts0111t1+Ts1111t1+",
"Ts0000t2+Ts1000t2+Ts0100t2+Ts1100t2+Ts0010t2+Ts1010t2+Ts0110t2+Ts1110t2+Ts0001t2+Ts1001t2+Ts0101t2+Ts1101t2+Ts0011t2+Ts1011t2+Ts0111t2+Ts1111t2+", 
"Ts0000t3+Ts1000t3+Ts0100t3+Ts1100t3+Ts0010t3+Ts1010t3+Ts0110t3+Ts1110t3+Ts0001t3+Ts1001t3+Ts0101t3+Ts1101t3+Ts0011t3+Ts1011t3+Ts0111t3+Ts1111t3+",
"Ts0000t4+Ts1000t4+Ts0100t4+Ts1100t4+Ts0010t4+Ts1010t4+Ts0110t4+Ts1110t4+Ts0001t4+Ts1001t4+Ts0101t4+Ts1101t4+Ts0011t4+Ts1011t4+Ts0111t4+Ts1111t4+",
"Rs0000+Rs1000+Rs0100+Rs1100+Rs0010+Rs1010+Rs0110+Rs1110+Rs0001+Rs1001+Rs0101+Rs1101+Rs0011+Rs1011+Rs0111+Rs1111"), 

"Ts0000=Ts0000t1+Ts0000t2+Ts0000t3+Ts0000t4",

"Ts1000=Ts1000t1+Ts1000t2+Ts1000t3+Ts1000t4",

"Ts0100=Ts0100t1+Ts0100t2+Ts0100t3+Ts0100t4",

"Ts0010=Ts0010t1+Ts0010t2+Ts0010t3+Ts0010t4",

"Ts0001=Ts0001t1+Ts0001t2+Ts0001t3+Ts0001t4",

"Ts1100=Ts1100t1+Ts1100t2+Ts1100t3+Ts1100t4",

"Ts1010=Ts1010t1+Ts1010t2+Ts1010t3+Ts1010t4",

"Ts1001=Ts1001t1+Ts1001t2+Ts1001t3+Ts1001t4",

"Ts0110=Ts0110t1+Ts0110t2+Ts0110t3+Ts0110t4",

"Ts0101=Ts0101t1+Ts0101t2+Ts0101t3+Ts0101t4",

"Ts0011=Ts0011t1+Ts0011t2+Ts0011t3+Ts0011t4",

"Ts1101=Ts1101t1+Ts1101t2+Ts1101t3+Ts1101t4",

"Ts1110=Ts1110t1+Ts1110t2+Ts1110t3+Ts1110t4",

"Ts1011=Ts1011t1+Ts1011t2+Ts1011t3+Ts1011t4",

"Ts0111=Ts0111t1+Ts0111t2+Ts0111t3+Ts0111t4",

"Ts1111=Ts1111t1+Ts1111t2+Ts1111t3+Ts1111t4"

)


sir_simu <- build_simulator(reactions, functions=simplification)

#####
###Initial States


initialStates <- c(S =70000, 
                   
                   Es0000 = 20000, Es1000 = 0, Es0100 = 0, Es1100 = 0, Es0010 = 0, Es1010 = 0,Es0110 = 0,Es1110 = 0,
                   Es0001 = 20000, Es1001 = 0, Es0101 = 0, Es1101 = 0, Es0011 = 0, Es1011 = 0,Es0111 = 0,Es1111 = 0,
                   
                   
                   Ls0000 = 20000, Ls1000 = 0, Ls0100 = 0, Ls1100 = 0, Ls0010 = 0, Ls1010 = 0,Ls0110 = 0,Ls1110 = 0,
                   Ls0001 = 20000, Ls1001 = 0, Ls0101 = 0, Ls1101 = 0, Ls0011 = 0, Ls1011 = 0,Ls0111 = 0,Ls1111 = 0,
                  
                   
                   Is0000 = 20000, Is1000 = 0, Is0100 = 0, Is1100 = 0, Is0010 = 0, Is1010 = 0,Is0110 = 0,Is1110 = 0,
                   Is0001 = 20000, Is1001 = 0, Is0101 = 0, Is1101 = 0, Is0011 = 0, Is1011 = 0,Is0111 = 0,Is1111 = 0,
                  
                   
                   Ts0000t1 = 0, Ts1000t1 = 0, Ts0100t1 = 0, Ts1100t1 = 0, Ts0010t1 = 0, Ts1010t1 = 0,Ts0110t1 = 0,Ts1110t1 = 0,
                   Ts0001t1 = 0, Ts1001t1 = 0, Ts0101t1 = 0, Ts1101t1 = 0, Ts0011t1 = 0, Ts1011t1 = 0,Ts0111t1 = 0,Ts1111t1 = 0,
                   Ts0000t2 = 0, Ts1000t2 = 0, Ts0100t2 = 0, Ts1100t2 = 0, Ts0010t2 = 0, Ts1010t2 = 0,Ts0110t2 = 0,Ts1110t2 = 0,
                   Ts0001t2 = 0, Ts1001t2 = 0, Ts0101t2 = 0, Ts1101t2 = 0, Ts0011t2 = 0, Ts1011t2 = 0,Ts0111t2 = 0,Ts1111t2 = 0,
                   Ts0000t3 = 0, Ts1000t3 = 0, Ts0100t3 = 0, Ts1100t3 = 0, Ts0010t3 = 0, Ts1010t3 = 0,Ts0110t3 = 0,Ts1110t3 = 0,
                   Ts0001t3 = 0, Ts1001t3 = 0, Ts0101t3 = 0, Ts1101t3 = 0, Ts0011t3 = 0, Ts1011t3 = 0,Ts0111t3 = 0,Ts1111t3 = 0,
                   Ts0000t4 = 0, Ts1000t4 = 0, Ts0100t4 = 0, Ts1100t4 = 0, Ts0010t4 = 0, Ts1010t4 = 0,Ts0110t4 = 0,Ts1110t4 = 0,
                   Ts0001t4 = 0, Ts1001t4 = 0, Ts0101t4 = 0, Ts1101t4 = 0, Ts0011t4 = 0, Ts1011t4 = 0,Ts0111t4 = 0,Ts1111t4 = 0,
                 
                   Rs0000 = 20000, Rs1000 = 0, Rs0100 = 0, Rs1100 = 0, Rs0010 = 0, Rs1010 = 0,Rs0110 = 0,Rs1110 = 0,
                   Rs0001 = 20000, Rs1001 = 0, Rs0101 = 0, Rs1101 = 0, Rs0011 = 0, Rs1011 = 0,Rs0111 = 0,Rs1111 = 0
                   
)
#####

###Duration of the simuI_isorn = 0, I_isoriso=0, I_isorr=0, I_isore=0, I_isorp=0, I_isorb=0, I_isord=0ation
time <- c(0,20,40)
dT <- round(1/365, 4) # use a daily time step

#####
###Parameters value
theta <- list(pi = 0.01856/365.25, #birth
             
              alphas0000 = 1.779603e-06, alphas1000 = 1.779603e-05, alphas0100 = 1.779603e-05, alphas0010 = 1.779603e-05, 
              alphas1100 = 1.779603e-05, alphas1010 = 1.779603e-05, alphas0110 = 1.779603e-05, alphas1110 = 1.779603e-05,
              alphas0001 = 1.779603e-06, alphas1001 = 1.779603e-05, alphas0101 = 1.779603e-05, alphas0011 = 1.779603e-05, 
              alphas1101 = 1.779603e-05, alphas1011 = 1.779603e-05, alphas0111 = 1.779603e-05, alphas1111 = 1.779603e-05,#infection 
              
              
              epsilon = 0.2/365.25, #progression rate from early latent to late latent
             
              delta1 = 0.4/365.25, delta2 = 0.0008/365.25, #progression to active disease
              
              lambda1=0.21/365.25, lambda2=0.21/365.25, #reinfection to early latent stage
              
              psi=0, #detected rate
              
              sigma=1, #tested to know which strain is it
              
              omega=1, #test is right
              
              iota1=1/4, iota2=1/4, iota3=1/4, iota4=1/4,#rate wrong result for the 0,1,2 treatment
              
              ####treated when the result is good
              betas0000t1=1,betas1000t1=0,betas0100t1=0,betas1100t1=0,betas0010t1=1,betas1010t1=0,betas0110t1=0,betas1110t1=0,
              betas0001t1=1,betas1001t1=0,betas0101t1=0,betas1101t1=0,betas0011t1=1,betas1011t1=0,betas0111t1=0,betas1111t1=0,
              
              betas0000t2=0,betas1000t2=1,betas0100t2=1,betas1100t2=1,betas0010t2=0,betas1010t2=0,betas0110t2=0,betas1110t2=0,
              betas0001t2=0,betas1001t2=1,betas0101t2=1,betas1101t2=1,betas0011t2=0,betas1011t2=0,betas0111t2=0,betas1111t2=0,
              
              betas0000t3=0,betas1000t3=1,betas0100t3=1,betas1100t3=1,betas0010t3=0,betas1010t3=1,betas0110t3=1,betas1110t3=1,
              betas0001t3=0,betas1001t3=0,betas0101t3=0,betas1101t3=0,betas0011t3=0,betas1011t3=0,betas0111t3=0,betas1111t3=0,
              
              betas0000t4=0,betas1000t4=0,betas0100t4=0,betas1100t4=0,betas0010t4=0,betas1010t4=1,betas0110t4=1,betas1110t4=1,
              betas0001t4=0,betas1001t4=1,betas0101t4=1,betas1101t4=1,betas0011t4=0,betas1011t4=1,betas0111t4=1,betas1111t4=1,
            
              
              
              ####reinfection 
              LambdaL0000=0,LambdaL1000=0,LambdaL0100=0,LambdaL1100=0,LambdaL0010=0,LambdaL1010=0,LambdaL0110=0,LambdaL1110=0,
              LambdaL0001=0,LambdaL1001=0,LambdaL0101=0,LambdaL1101=0,LambdaL0011=0,LambdaL1011=0,LambdaL0111=0,LambdaL1111=0,
              
              LambdaR0000=0,LambdaR1000=0,LambdaR0100=0,LambdaR1100=0,LambdaR0010=0,LambdaR1010=0,LambdaR0110=0,LambdaR1110=0,
              LambdaR0001=0,LambdaR1001=0,LambdaR0101=0,LambdaR1101=0,LambdaR0011=0,LambdaR1011=0,LambdaR0111=0,LambdaR1111=0,
              
              
              
              
              omicron=0.02/365.25, #self cured
            
              rho1000= 30/365.25, rho0100=30/365.25, rho0010=30/365.25, rho0001=30/365.25, #acquisition of a new resistance
              #rho100=0.8/365.25, rho010=0.8/365.25, rho001=0.08/365.25, #acquisition of a new resistance
              
              ita=0.03, ita_old= 0, ita_new=0,  #treatment shift due to resistance
              
              theta= 0.06/365.25, #death due to tuberculosis 
              #(1. Figures of the dead: a decade of tuberculosis mortality registrations in South Africa. https://journals.co.za/doi/epdf/10.7196/SAMJ.2019.v109i10.14073 doi:10.7196/SAMJ.2019.v109i10.14073.)
              
              phi=0.01/365.25, #re-activation from recovered to infected active stage
              
              mu= 0.00926) #natural death

#####
###In some simulations, the population size of a deme compartment may be zero before the upper time limit is reached, because of stochasticity or parameter values. 
#In this case, the simulation is considered to have failed and is halted.
#To bypass these failures, we can define the following wrapper:
safe_run <- function(f, ...) {
  out <- list()
  while(! length(out)) {out <- f(...)}
  out
}


#A safe version of our simulator sir_simu() is then:
safe_sir_simu <- function(...) safe_run(sir_simu, ...)




###Direct method
traj_dm <- safe_sir_simu(
  paramValues = theta,
  initialStates = initialStates,
  times = time,
  tau = dT,
  method = "mixed",
  nTrials=5)
#####

####Visualisation
colnames(traj_dm$traj)
tail(traj_dm$traj)


layout(matrix(1:1,1,1))
plot(traj_dm$traj[,1],log(traj_dm$traj[,4]), type = "l", ylim=c(0,12) , xlab="time (years)", ylab="number of individuals")
lines(traj_dm$traj[,1], log(traj_dm$traj[,21]), type= "l", col = "blue")
lines(traj_dm$traj[,1], log(traj_dm$traj[,13]), type= "l", col = "purple")
lines(traj_dm$traj[,1], log(traj_dm$traj[,5]), type= "l", col = "red")
lines(traj_dm$traj[,1], log(traj_dm$traj[,29]), type= "l", col = "green")
lines(traj_dm$traj[,1], log(traj_dm$traj[,53]), type= "l", col = "green")




plot(traj_dm$traj[,1], traj_dm$traj[,49], type= "l", col = "green")


layout(matrix(1:4,2,2))
plot(traj_dm$traj[,1],traj_dm$traj[,4], type = "l")
plot(traj_dm$traj[,1], traj_dm$traj[,21], type= "l", col = "blue")
plot(traj_dm$traj[,1], traj_dm$traj[,13], type= "l", col = "purple")
plot(traj_dm$traj[,1], traj_dm$traj[,5], type= "l", col = "red")



Is000<-as.vector(rowSums(traj_dm$traj[,c(21,29,37)]))
Is100<-as.vector(rowSums(traj_dm$traj[,c(22,30,38)]))
Is010<-as.vector(rowSums(traj_dm$traj[,c(23,31,39)]))
Is110<-as.vector(rowSums(traj_dm$traj[,c(24,32,40)]))
Is001<-as.vector(rowSums(traj_dm$traj[,c(25,33,41)]))
Is101<-as.vector(rowSums(traj_dm$traj[,c(26,34,42)]))
Is011<-as.vector(rowSums(traj_dm$traj[,c(27,35,43)]))
Is111<-as.vector(rowSums(traj_dm$traj[,c(28,36,44)]))

plot(traj_dm$traj[,1], Is000, type= "l", col = "blue")
lines(traj_dm$traj[,1], Is100, type= "l", col = "red")
lines(traj_dm$traj[,1], Is010, type= "l", col = "green")
lines(traj_dm$traj[,1], Is110, type= "l", col = "black")
lines(traj_dm$traj[,1], Is001, type= "l", col = "purple")
lines(traj_dm$traj[,1], Is101, type= "l", col = "grey")
lines(traj_dm$traj[,1], Is011, type= "l", col = "orange")
lines(traj_dm$traj[,1], Is111, type= "l", col = "yellow")


ncol(traj_dm$traj)


Susceptible<-as.vector(traj_dm$traj[,4])
EarlyLatent<-as.vector(rowSums(traj_dm$traj[,5:12]))
LateLatent<-as.vector(rowSums(traj_dm$traj[,13:20]))
Infected<-as.vector(rowSums(traj_dm$traj[,21:44]))
Recovered<-as.vector(rowSums(traj_dm$traj[,44:68]))
Infected_Resistant<-as.vector(rowSums(traj_dm$traj[,c(22,23,24,25,26,27,28,30,31,32,33,34,35,36,38,39,40,41,42,43,44)]))

layout(matrix(1:1,1,1))
plot(traj_dm$traj[,1], log(Susceptible),ylim=c(0,12), type= "l", col = "black", ylab="log number of individuals", xlab="time (years)")
lines(traj_dm$traj[,1], log(EarlyLatent), type= "l", col = "red")
lines(traj_dm$traj[,1], log(LateLatent), type= "l", col = "purple")
lines(traj_dm$traj[,1], log(Infected), type= "l", col = "blue")
lines(traj_dm$traj[,1], log(Recovered), type= "l", col = "green")



layout(matrix(1:4,2,2))
plot(traj_dm$traj[,1],Susceptible, type = "l", xlab="time (years)")
plot(traj_dm$traj[,1], Recovered, type= "l", col = "green",xlab="time (years)")
plot(traj_dm$traj[,1], Infected, type= "l", col = "blue", xlab="time (years)")
plot(traj_dm$traj[,1], Infected_Resistant, type="l", col= "darkblue", xlab="time (years)")
