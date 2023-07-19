#This script varies the half-sat constant of phyto uptake of nutrients.

library(tidyverse)
library(ggpubr)
library(grid)

#time
dt = 0.01 # time step [1/100 of a day]
months = 600
timesteps = 3000 * months

#Area
fishpond_area_sequences <- list(seq(0.1, 0.1, by = 0.1))
patch_area_grid <- do.call(expand.grid, fishpond_area_sequences)
patch_area_grid$Var2 <- 0.34 - patch_area_grid$Var1
patch_area_grid$Var3 <- 0.66
#change order of columns
patch_area_grid <- patch_area_grid[,c(3,2,1)]
colnames(patch_area_grid) <- c("kbay", "estuary", "fishpond")
patch_area_list <- split(patch_area_grid, 1:nrow(patch_area_grid))
number_patches <- ncol(patch_area_grid)

#for later
patch_area_grid_m2 <- patch_area_grid * 4.14e+7 
colnames(patch_area_grid_m2) <- c("kbay", "estuary", "fishpond")

ks <- list(seq(0, 0.3, by =0.05), seq(0, 0.3, by =0.05), seq(0, 0.3, by =0.05))
ks_grid <- do.call(expand.grid, ks)
ks_grid <- ks_grid %>% filter(Var1==Var2) #only want patch 1 and 2 to have same ks
ks_list <- split(ks_grid, 1:nrow(ks_grid))

#parameter combinations
parameter_grid <- expand.grid(patch_area = patch_area_list,
                              ks = ks_list)


#rest of parameters
#outputs
outcome_nutrients <- matrix(0, ncol = number_patches, nrow = nrow(parameter_grid))
outcome_phyto <- matrix(0, ncol = number_patches, nrow = nrow(parameter_grid))
outcome_fish <- matrix(0, ncol = number_patches, nrow = nrow(parameter_grid))
outcome_harvest <- matrix(0, ncol = number_patches, nrow = nrow(parameter_grid))

#Rates
Vm = 2     # max phtyoplankton growth rate
Rm = 1.1     # max fish ingestion rate, 1.1
#ks = 0.31 #.31   # half-saturation constant for phytoplankton N uptake 
m = 0.1      # phytoplankton mortality
kg = 12      # grazing half saturation or (1/Ivlev constant)
lmbda = 0.1  # unassimilated phytoplankton fraction
g = 0.001      # fish mortality
fishing_effort = c(0, 0, 0)

#DISPERSAL
dispersal <- function(patch_area, number_patches, S) {
  D_common_pool <- matrix(rep(patch_area, number_patches), ncol = number_patches, byrow = TRUE) # common pool dispersal such that the movement is proportional to patch relative area
  require(ramify) # need this for matrix functions (e.g., triu)
  if(S > 0){
    tmp1 <- D_common_pool
    tmp1_triu <- triu(tmp1, 1) # set values in the upper triangular part of the matrix (in the 2x2 matrix this mean the cross patch movement)
    tmp1_tril <- tril(tmp1,-1) # set values in the lower triangular part of the matrix (in the 2x2 matrix this mean the cross patch movement)
    tmp1_offdiag <- tmp1_triu + tmp1_tril # Original off diagonal, or cross patch movement (used below)
  } else {
    return(D_common_pool)
  }
  for(D_row in 1:length(diag(tmp1))){ # Repeat for each patch, taking the length of the diagonal is equivalent to the total number of patches
    if(tmp1[D_row,D_row] > 0){ # If that site has positive habitat area  
      if((sum(tmp1[D_row, ]) - tmp1[D_row, D_row]) > 0){ # and the cross recruit sites have positive habitat area, this should always be true for our case
        # calculate the difference between larval pool recruitment and 100% self-recruitment
        diff_self_pool <- 1 - tmp1[D_row, D_row] # the 1 being 100% self-recruitment
        # add enhanced site-fidelity to common pool
        enhanced_site_fidelity <- tmp1[D_row, D_row] + (diff_self_pool * S) # tmp1[D_row, D_row] indicates self recruitment to a patch 
        # Standardize the proportion of the off diagonals such that they sum to 1 (see offdiag above)
        tmp1_offdiag_row_standardized <- tmp1_offdiag[D_row, ] / sum(tmp1_offdiag[D_row, ])      
        # Subtract from the common pool cross-recruit probability the standardized diff*proportion enhancement
        reduced_cross_recruit <- tmp1_offdiag[D_row, ] - (tmp1_offdiag_row_standardized * (diff_self_pool*S))
        # Insert the new values
        tmp1[D_row, ] <- reduced_cross_recruit 
        tmp1[D_row, D_row] <- enhanced_site_fidelity
      }
    } 
  }
  return(tmp1) 
}

S <- 0.5
S_recruit <- 0.5

ingestion <- function(predation_type, kg, P){
  if (predation_type == 'quadratic'){
    I = P^2/(kg^2+P^2)
  }   
  else if (predation_type == 'Ivlev'){
    I = 1-exp(-kg^P) 
  }
  else if (predation_type == 'M-P'){
    I = kg^-1*P*(1-exp(-kg^-1*P))
  }
  else {
    I = P/(kg+P) 
  }
  return(I)
}
predation_type = 'M-P'


phyto_growth_rate <- function(N,Vm,ks){
  mu = Vm*N/(ks+N) # growth rate [day^-1]
  return(mu)
}

#starting values
##100 timesteps per day, 30 days per month, so every 3000 timesteps starts a new month

phytoplankton <- array(NA, dim = c(timesteps, number_patches))

nutrients <- array(NA, dim = c(timesteps, number_patches))

#nutrient_inputs <- c(0.07, 0.13, 0.13)
fish <- array(NA, dim = c(timesteps, number_patches))

fraction_harvested <- array(NA, dim = c(timesteps, number_patches))
harvest <- array(NA, dim = c(timesteps, number_patches))

prod <- array(NA, dim = c(timesteps, number_patches))
pmort <- array(NA, dim = c(timesteps, number_patches))
grazing <- array(NA, dim = c(timesteps, number_patches))
fmort <- array(NA, dim = c(timesteps, number_patches))

#nutrient inputs
l2 <- c(rep(0.1377, times = 6), rep(0.1583, times = 6))
l2 <- c(rep(l2, times = months/12)) #l2 <- c(rep(0.1377, times = monthly_timesteps)) #this is the max from the NERR 2021 Wai 2 dataset
l3 <- l2
l1 <- l2 * 0.06 #6% of what goes into the fishpond, calculated from NERR 2021 R9 data where the max was 6% of the max of Wai 2




for (iter in 1:nrow(parameter_grid)) {
  starting_nutrients <- c(0.1377*0.06,0.1377,0.1377)
  starting_phyto <- c(0.65,1.23,1.23)
  starting_fish <- c(10,10,10)
  fecundity_rate <-  c(0.001, 0.00025, 0.00025) #0.001 for p1
  nutrient_inputs <- cbind(l1, l2, l3)
  
  for (i in 1:number_patches) {
    if (as.numeric(parameter_grid[['patch_area']][[iter]])[i] == 0) {
      starting_nutrients[i] <- 0
      nutrient_inputs[,i] <- 0 
      starting_fish[i] <- 0
      starting_phyto[i] <- 0
      fecundity_rate[i] <- 0
    }
  }
  phytoplankton[1,] <- starting_phyto
  nutrients[1,] <- starting_nutrients
  fish[1,] <- starting_fish
  
  #scaled areas for dispersal
  kbay_area_scaled <- as.numeric(parameter_grid[['patch_area']][[iter]])[1] / (as.numeric(parameter_grid[['patch_area']][[iter]])[1] + as.numeric(parameter_grid[['patch_area']][[iter]])[2])
  estuary_area_scaled = as.numeric(parameter_grid[['patch_area']][[iter]])[2] / (as.numeric(parameter_grid[['patch_area']][[iter]])[1] + as.numeric(parameter_grid[['patch_area']][[iter]])[2])
  patch_area_scaled = c(kbay_area_scaled, estuary_area_scaled)
  number_patches_scaled = length(patch_area_scaled)
  patch_area_m2 <- as.numeric(parameter_grid[['patch_area']][[iter]]) * 4.14e+7 
  
  #dispersal rates
  tmp1 <- dispersal(patch_area_scaled, number_patches_scaled, S)
  id_matrix_1 <- matrix(0, number_patches, number_patches_scaled)
  diag(id_matrix_1) <- 1
  id_matrix_2 <- matrix(0, number_patches_scaled, number_patches)
  diag(id_matrix_2) <- 1
  adult_disp <- id_matrix_1 %*% tmp1 %*% id_matrix_2
  for (i in number_patches_scaled:number_patches) {
    if (i > 2)
      adult_disp[i, i] <- 1
  }
  #set diagonal to be zeros - new matrix
  adult_immigration <- adult_disp
  for (i in 1:number_patches) {
    adult_immigration[i, i] <- 0
  }
  adult_emigration <- c((1-adult_disp[1,1]), (1-adult_disp[2,2]), (1-adult_disp[3,3]))
  recruit_migration <- dispersal(as.numeric(parameter_grid[['patch_area']][[iter]]), number_patches, S_recruit)
  
  #ACTUAL LOOP
  for(t in 2:timesteps){
    if ((t %% 3000) == 0) {
      nutrients[t-1,] <- nutrients[t-1,] + nutrient_inputs[t/3000,]
    }
    prod[t-1,] = phytoplankton[t-1,] * ((Vm * nutrients[t-1,])/(as.numeric(parameter_grid[['ks']][[iter]]) + nutrients[t-1,]))
    pmort[t-1,] = phytoplankton[t-1,] * m
    grazing[t-1,] = Rm * ingestion(predation_type, kg, phytoplankton[t-1,]) * fish[t-1,]
    fmort[t-1,] = g*fish[t-1,]
    
    phytoplankton[t,] = phytoplankton[t-1,] + (prod[t-1,] - pmort[t-1,] - grazing[t-1,]) *dt
    nutrients[t,] = nutrients[t-1,] + ((-prod[t-1,] + pmort[t-1,] + (fmort[t-1,]*.75) + lmbda*grazing[t-1,]) *dt)
    fish[t-1,] = fish[t-1,] + (
      (fish[t-1,] %*% (adult_immigration/365)) - 
        (fish[t-1,] * (adult_emigration/365)) +
        ((fish[t-1,] * fecundity_rate) %*% (recruit_migration/365)) +
        ((1-lmbda) * grazing[t-1,]) -
        fmort[t-1,]
    ) * dt
    
    fraction_harvested[t-1,] = (1 - exp(-fishing_effort)) * dt
    harvest[t-1,] = fish[t-1,] * fraction_harvested[t-1,] * patch_area_m2
    fish[t,] = fish[t-1,] * (1 - fraction_harvested[t-1,])
  }
  
  outcome_nutrients[iter, ] <- colMeans(nutrients[t-360000:t,]) #nutrients[t,]
  outcome_phyto[iter, ] <- colMeans(phytoplankton[t-360000:t,]) #phytoplankton[t,]
  outcome_fish[iter, ] <- colMeans(fish[t-360000:t,]) #fish[t,]
  outcome_harvest[iter, ] <-colMeans(harvest[t-360000:(t-1),]) #harvest[t-1,]
}


#plots
colnames(outcome_nutrients) <- c("p1_nutrients", "p2_nutrients", "p3_nutrients")
colnames(outcome_phyto) <- c("p1_phyto", "p2_phyto", "p3_phyto")
colnames(outcome_fish) <- c("p1_fish", "p2_fish", "p3_fish")
colnames(outcome_harvest) <- c("p1_harvest", "p2_harvest", "p3_harvest")

outcome <- cbind(parameter_grid, outcome_nutrients, outcome_phyto, outcome_fish, outcome_harvest)
outcome$area_p1 <- map_dbl(outcome$patch_area, 1)
outcome$area_p2 <- map_dbl(outcome$patch_area, 2)
outcome$area_p3 <- map_dbl(outcome$patch_area, 3)
outcome$ks_p1 <- map_dbl(outcome$ks, 1)
outcome$ks_p2 <- map_dbl(outcome$ks, 2)
outcome$ks_p3 <- map_dbl(outcome$ks, 3)

outcome$area_p1_m <- outcome$area_p1 * 4.14e+7
outcome$area_p2_m <- outcome$area_p2 * 4.14e+7
outcome$area_p3_m <- outcome$area_p3 * 4.14e+7

#population x area to get total g in the patch. convert to kg.
outcome$p1_fish_kg <- outcome$p1_fish * outcome$area_p1_m / 1000
outcome$p2_fish_kg <- outcome$p2_fish * outcome$area_p2_m / 1000
outcome$p3_fish_kg <- outcome$p3_fish * outcome$area_p3_m / 1000


#weighted avg by area?? idk
outcome$est_fish <- ((outcome$p1_fish * outcome$area_p1) + 
                       (outcome$p2_fish * outcome$area_p2)) / (outcome$area_p1 + outcome$area_p2)

outcome$est_phyto <- ((outcome$p1_phyto * outcome$area_p1) + 
                        (outcome$p2_phyto * outcome$area_p2)) / (outcome$area_p1 + outcome$area_p2)

outcome$est_nutrients <- ((outcome$p1_nutrients * outcome$area_p1) + 
                            (outcome$p2_nutrients * outcome$area_p2)) / (outcome$area_p1 + outcome$area_p2)

outcome <- replace(outcome, is.na(outcome), 0)

df <- apply(outcome,2,as.character)
write.csv(df, file='ks_sens.csv')

outcome <- read.csv("ks_sens.csv")


p1<-ggplot(data=outcome %>% filter(ks_p3==0.1 & ks_p1>0), aes(x=ks_p1, y=est_nutrients)) +
  geom_line(size=2) +
  theme_minimal() +
  theme(plot.title = element_text(size = 16), 
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position = "right",
        plot.margin = margin(1,1,1,1, "cm"),
        strip.text.x = element_text(size = 20), 
        axis.text.x = element_text(vjust = 0.5, hjust=1)) +
  labs(x = "Bay half-saturation constant", y= bquote('Bay nutrients density'~(g/m^3)))

p2<-ggplot(data=outcome %>% filter(ks_p3==0.1 & ks_p1>0), aes(x=ks_p1, y=est_phyto)) +
  geom_line(size=2) +
  theme_minimal() +
  theme(plot.title = element_text(size = 16), 
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position = "right",
        plot.margin = margin(1,1,1,1, "cm"),
        strip.text.x = element_text(size = 20), 
        axis.text.x = element_text(vjust = 0.5, hjust=1)) +
  labs(x = "Bay half-saturation constant", y= bquote('Bay phytoplankton density'~(g/m^2)))

p3<-ggplot(data=outcome %>% filter(ks_p3==0.1 & ks_p1>0), aes(x=ks_p1, y=est_fish)) +
  geom_line(size=2) +
  theme_minimal() +
  theme(plot.title = element_text(size = 16), 
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position = "right",
        plot.margin = margin(1,1,1,1, "cm"),
        strip.text.x = element_text(size = 20), 
        axis.text.x = element_text(vjust = 0.5, hjust=1)) +
  labs(x = "Bay half-saturation constant", y= bquote('Bay fish density'~(g/m^2)))


#loko ia
p4<-ggplot(data=outcome %>% filter(ks_p1==0.1 & ks_p3>0), aes(x=ks_p3, y=p3_nutrients)) +
  geom_line(size=2) +
  theme_minimal() +
  theme(plot.title = element_text(size = 16), 
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position = "right",
        plot.margin = margin(1,1,1,1, "cm"),
        strip.text.x = element_text(size = 20), 
        axis.text.x = element_text(vjust = 0.5, hjust=1)) +
  labs(x = "Loko i'a half-saturation constant", y= bquote("Loko i'a nutrients density"~(g/m^3)))

p5<-ggplot(data=outcome %>% filter(ks_p1==0.1 & ks_p3>0), aes(x=ks_p3, y=p3_phyto)) +
  geom_line(size=2) +
  theme_minimal() +
  theme(plot.title = element_text(size = 16), 
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position = "right",
        plot.margin = margin(1,1,1,1, "cm"),
        strip.text.x = element_text(size = 20), 
        axis.text.x = element_text(vjust = 0.5, hjust=1)) +
  labs(x = "Loko i'a half-saturation constant", y= bquote("Loko i'a phytoplankton density"~(g/m^3)))

p6<-ggplot(data=outcome %>% filter(ks_p1==0.1 & ks_p3>0), aes(x=ks_p3, y=p3_fish)) +
  geom_line(size=2) +
  theme_minimal() +
  theme(plot.title = element_text(size = 16), 
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position = "right",
        plot.margin = margin(1,1,1,1, "cm"),
        strip.text.x = element_text(size = 20), 
        axis.text.x = element_text(vjust = 0.5, hjust=1)) +
  labs(x = "Loko i'a half-saturation constant", y= bquote("Loko i'a fish density"~(g/m^3)))

plot<-ggarrange(p1+rremove("xlab"),p2+rremove("xlab"),p3+rremove("xlab"),
                p4+rremove("xlab"),p5+rremove("xlab"),p6+rremove("xlab"),nrow=2,ncol=3)

annotate_figure(plot, bottom = text_grob(bquote("Half saturation constant"~(k[s])), size=20))


