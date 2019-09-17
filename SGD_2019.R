library("readr")
library("dplyr")
library("tidyr")
library("mvtnorm")
library("numDeriv")
library("nFactors")
library("psych")
library("corrplot")
library("formattable")
library("kableExtra")
library("ggplot2")
library("ggrepel")
library("directlabels")
library("tidyverse")
library("gridExtra")
library("grid")
library("ggpubr")

Sys.setenv('R_MAX_VSIZE'=16000000000)
options(scipen=999) #deafult = 0

path <- "/Users/Liang/Dropbox/Yuliang Weng - Dynamic Factor Models/"
data_path <- paste(path, "data/nba_data", sep = "/")

year_directories <- dir(data_path, full.names = TRUE)
year_directories <- year_directories[5]


spec <- cols(
  Rk = col_double(),
  G = col_double(),
  Date = col_date(format = ""),
  Age = col_character(),
  Tm = col_character(),
  `Unnamed: 5` = col_character(),
  Opp = col_character(),
  `Unnamed: 7` = col_character(),
  GS = col_character(),
  MP = col_character(),
  FG = col_double(),
  FGA = col_double(),
  `FG%` = col_double(),
  `3P` = col_double(),
  `3PA` = col_double(),
  `3P%` = col_double(),
  FT = col_double(),
  FTA = col_double(),
  `FT%` = col_double(),
  ORB = col_double(),
  DRB = col_double(),
  TRB = col_double(),
  AST = col_double(),
  STL = col_double(),
  BLK = col_double(),
  TOV = col_double(),
  PF = col_double(),
  PTS = col_double(),
  GmSc = col_double(),
  `+/-` = col_double()
)
vars <- names(spec$cols)


player_data_2019 <- NULL
for (y in year_directories) {
  player_file <- dir(y)
  player_name <- gsub("_|regular|playoff|\\.csv", "", player_file)
  regular_season <- grepl("regular", player_file)
  for (j in seq_along(player_file)) {
    c_data <- readr::read_csv(paste(y, player_file[j], sep = "/"), col_type = spec)
    ## If any variables are missing from the csv enrich c_data with those columns and give NA values
    ind <- vars %in% names(c_data)
    if (!all(ind)) {
      c_data[vars[!ind]] <- NA
    }
    ## Just making sure that the variables appear in the right columns
    c_data <- c_data[vars]
    c_data$regular_season <- regular_season[j]
    c_data$player <- player_name[j]
    player_data_2019 <- rbind(player_data_2019, c_data)
  }
}

toAge <- function(x) {
  x <- as.numeric(strsplit(x, "-")[[1]])
  x[1] + x[2]/365
}

player_data_2019 <- player_data_2019 %>%
  dplyr::rename(`Home` = `Unnamed: 5`, `game_result` = `Unnamed: 7`) %>%
  dplyr::mutate(Home = (Home == "@") %in% TRUE) %>%
  dplyr::mutate(Age = toAge(Age))
## Summaries

## Number of unique players in the database
player_data_2019 %>% count(player) %>% nrow
## 874 players
nrow(player_data_2019)


select_summaries <- c("Date","GS","FGA","FG%","3PA","3P%","FTA","FT%","AST","STL","ORB","DRB","PTS","+/-","player")
len_sum <- length(select_summaries) - 3
## Select relevant summaries and only keep those matches that players played
player_data_active_2019 <- player_data_2019 %>% dplyr::select(select_summaries) %>% filter(GS==1 | GS==0) %>% mutate_if(is.numeric, funs(replace_na(., 0)))  %>% arrange(Date)
# replace(player_data_active, is.na("FG%"), 0)
player_data_active_2019 %>% count(player) %>% nrow

player_stat_2019 <- player_data_active_2019 %>% count(player)
no_player_2019 <- length(player_stat_2019$player)

# Arrange 
time_point_2019 <- player_data_active_2019 %>% group_split(Date)
time_point_all_2019 <- do.call("rbind", time_point_2019)


############################################################################################################################################
#####################################################################################################################################
########################################## Maximum likelihood estimation ##############################################################
#####################################################################################################################################
############################################################################################################################################
options(digits.secs = 4) 

start.time <- Sys.time()

start <- c(rep(0, len_sum), rep(1, len_sum), log(2))

MLE <- function(parameter, data1){
  parameters <- c(0, parameter)
  alpha <- parameters[1:len_sum] 
  beta <- parameters[(len_sum+1):(2*len_sum)]
  psi <- parameters[(2*len_sum+1)]
  vc <- tcrossprod(beta)
  diag(vc) <- diag(vc) + exp(psi)
  ivc <- 0
  ivc_try <- try(solve(vc))
  
  if("try-error" %in% class(ivc_try)){
    return(-Inf)
  } else{
    ivc = ivc_try
  }
  
  # no_player should be p_t 
  n_player <- sapply(data1, nrow)
  part0 <- sum(n_player) * (len_sum-1) * log(2 * pi)
  part1 <- sum(n_player) * (len_sum-1) * psi
  part2 <- sum(n_player) * log(sum(beta^2) + exp(psi))
  part3 <- 0
  for(tp in data1){
    data_mat <- data.matrix(tp[,3:(2+len_sum)])
    x_a <- sweep(data_mat, 2, alpha, "-")
    part3 <- part3 + sum(diag(x_a %*% ivc %*% t(x_a)))
  }
  loglik <- -0.5 * (part0 + part1 + part2 + part3)
  return(-1 * loglik)    
  
}

MLE_results_2019 <- optim(par = as.vector(start[-1]), MLE, data1 = time_point_2019, method = "BFGS" )
MLE_para_2019 <- as.matrix(MLE_results_2019$par)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

#########################################################################################################
########### Compare learning rate of MLE & SGD ###################################
#########################################################################################################

start <- c(0.05, 0.1, 0.0003, 1, 1)
#start <- c(1,1,0.03,1,2)

MLE_learningRate_2019 <- optim(par = as.vector(start), MLE_learningRate, data1 = time_point_2019,results=MLE_results_2019$par)
MLE_LRpara_2019 <- as.matrix(MLE_learningRate_2019$par)
MLE_LRpara_2019


############################################################################################################################################
#####################################################################################################################################
########################################## Stochastic Gradient Descent ##############################################################
#####################################################################################################################################
############################################################################################################################################
start.time <- Sys.time()
## Function for intialise dataframe for store each player's abilities over time
IntialiseList <- function(player){
  df = data.frame("Date"=as.Date(time_point_all_2019[1,1][[1]],"%y-%m-%d"),"Abilities"=0)
  #df
  #df= df %>% mutate("Date" = col_date(format = "")) %>% mutate("Abilities" = col_double())
  df$Date = as.Date(df$Date,"%y-%m-%d")
  #typeof(df$Date)
  df[-1,]
}

returnZero <- function(x){
  if(is.na(x) == FALSE){
    return(x)
  }
  else{
    return(0)
  }
}

# Store team names for each player
team_name_2019 <- data.frame(Player_Name = player_stat_2019$player, Team_Name = rep(" ", no_player_2019), stringsAsFactors=FALSE )

for(i in 1:nrow(team_name_2019)){
  name = team_name_2019[i,1]
  team_name_2019[i,2] = player_data_2019[which(player_data_2019$player == name),5][1,1]
}


# Make a list of player names and their abilities 
player_ability_2019 <- rep(0,no_player_2019)
names(player_ability_2019) <- player_stat_2019$player 

abilityList_2019 <- lapply(player_stat_2019$player, IntialiseList)
names(abilityList_2019) <- player_stat_2019$player


# Initialisation of parameters and learning rate
alpha = matrix(rep(0, len_sum), nrow = len_sum, ncol = 1)
beta = matrix(rep(1, len_sum), nrow = len_sum, ncol = 1)
sigmaS = 2
psi = log(sigmaS)


alpha_lRate_0 = MLE_LRpara_2019[1,1]  
beta_lRate_0 = MLE_LRpara_2019[2,1]
sigmaS_lRate_0 = MLE_LRpara_2019[3,1]
lrate_a = MLE_LRpara_2019[4,1]
lrate_c = MLE_LRpara_2019[5,1]

#0.054288121
#0.029594859
#0.002879629
#1.075588335
#0.955340221

# Player's ability is calculated as: current updated ability at time t * factor + (1-factor) * ability at time t-1.
player_update_factor = 0.05

lrate_alpha <- c(alpha_lRate_0)
lrate_beta <- c(beta_lRate_0)
lrate_sigma <- c(sigmaS_lRate_0)

j <- 0
# For every time point, iteratively update parameters alpha, beta and sigma^2, and update player's abilities
for(tp in time_point_2019){
  
  alpha_lRate = 1
  beta_lRate = 1
  sigmaS_lRate = 1
  # Update learning rate
  if(j == 0){
    alpha_lRate = alpha_lRate_0
    beta_lRate = beta_lRate_0
    sigmaS_lRate = sigmaS_lRate_0
  } else{
    alpha_lRate = alpha_lRate_0 * (1 + lrate_a * alpha_lRate_0 * j)^(-lrate_c)
    lrate_alpha <- c(lrate_alpha,alpha_lRate)
    beta_lRate = beta_lRate_0 * (1 + lrate_a * beta_lRate_0 * j)^(-lrate_c)
    lrate_beta <- c(lrate_beta,beta_lRate)
    sigmaS_lRate = sigmaS_lRate_0 * (1 + lrate_a * sigmaS_lRate_0 * j)^(-lrate_c)
    lrate_sigma <- c(lrate_sigma, sigmaS_lRate)
  }
  
  j <- j+1
  # Retrieve columns relevant to selected summaries
  data_mat <- data.matrix(tp[,3:(2+len_sum)])
  
  len_tp <- nrow(tp)
  ########################### Update parameters ###################################
  
  # Variance-covariance matrix & its Inverse
  covar = (beta %*% t(beta))
  diag(covar) = diag(covar) + exp(psi)
  icovar = solve(covar) 
  
  ######## Update ALPHA #########
  
  sum_matrix_a = icovar %*% t(sweep(data_mat, 2, alpha, "-"))   # X-alpha
  alpha = alpha + (alpha_lRate * rowSums(sum_matrix_a))
  alpha[1] <- 0
  ########   Update BETA  #########
  x_a = sweep(data_mat, 2, alpha, "-")     ## Subtract alpha from summaries
  beta_part1 = -len_tp * beta / (sum(beta^2) + exp(psi))  ## First part of beta's derivative
  beta_part2 = matrix(c(rep(0,len_sum)), nrow=len_sum, ncol=1)      ## Initialise a empty vector at time t to store beta
  
  for(i in 1:len_sum){
    # Calculate the derivative of var-covar matrix w.r.t each beta
    cov_deri = matrix(0L, nrow = len_sum, ncol = len_sum)
    cov_deri[i,] = cov_deri[i,] + beta
    cov_deri[,i] = cov_deri[,i] + beta
    # Part 2 of beta's derivative
    beta_ts = 0.5 * ( x_a %*% icovar %*% cov_deri %*% icovar %*% t(x_a) )
    beta_part2[i] = sum(diag(beta_ts))
  }
  beta_deri = beta_part1 + beta_part2
  beta = beta + beta_lRate * beta_deri
  
  ########   Update SIGMA  #########
  sig_part1 = -len_tp * (len_sum-1) / (2 * exp(psi))
  sig_part2 = -0.5 * len_tp / (sum(beta^2) + exp(psi))
  sig_part3 = 0.5 * sum(diag(x_a %*% icovar %*% icovar %*% t(x_a)))
  
  sigmaS_t = exp(psi) * (sig_part1 + sig_part2 + sig_part3)    
  psi = psi + sigmaS_lRate*sigmaS_t
  
  ######## Update Player's Abilities  #############
  pa_2019 = colSums(t(sweep(data_mat,2,alpha,"-"))  *as.vector(beta)) / as.numeric(t(beta) %*% beta + exp(psi))
  names(pa_2019) = tp$player
  for(player in tp$player){
    new_abilities = returnZero((1-player_update_factor)*player_ability_2019[player][[1]]) + returnZero(player_update_factor*pa_2019[player][[1]])
    player_ability_2019[player][[1]] = new_abilities
    abilityList_2019[player][[1]] <- rbind(abilityList_2019[player][[1]], data.frame("Date"=as.Date(tp[1,1][[1]],"%y-%m-%d"),"Abilities"=new_abilities))
  }
}

################ Create a dataframe for the result ###########################3
df_2019 <- data.frame(player_ability_2019) 
df_2019 <- data.frame("Player_Name"=rownames(df_2019),"Player_Abilities"=df_2019$player_ability)
df_2019 <- df_2019[order(-df_2019$Player_Abilities),]

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

######################################################################
############### Team abilities #################################
###############################################################
Team_Name_2019 <- c()

# match player's team
for(i in 1:nrow(df_2019)){
  temp = team_name_2019[which(team_name_2019$Player_Name==df_2019[i,1]),][1,2]
  Team_Name_2019 <- c(Team_Name_2019, temp)
}
# assign team name to each player
df_2019$Team_Name <- Team_Name_2019
conference <- rep("Eastern", 30)
conference[c(7,8,10,11,13,14,15,18,19,21,24,25,26,27,29)] <- "Western"
# Group by team name, and sum all player's abilities for each team
Team_Abilities_2019 <- setNames(aggregate(df_2019$Player_Abilities, by=list(Category=df_2019$Team_Name), FUN=sum), c("Team", "Abilities"))
Team_Abilities_2019$Conference <- conference


ggbarplot(Team_Abilities_2019, x = "Team", y = "Abilities", fill = "Conference", color = "white", palette = "jco", 
          sort.val = "desc", sort.by.groups = FALSE, x.text.angle = 90, ggtheme = theme_pubclean(), ) + 
          font("x.text", size = 12, vjust = 0.5) + labs(title = "Team abilities for the season 2018-19") +
  theme(legend.position = "bottom",text=element_text(size=17),legend.title=element_text(size=26))

#Abilities of each conference
aggregate(Team_Abilities_2019$Abilities, by=list(Category=Team_Abilities_2019$Conference), FUN=sum)

########################################################################
############## Plot Learning rate curves ############################
################################################################################

rownames(df_2019) <- 1:nrow(df_2019)
par(mfrow=c(1,3)) 
plot(lrate_alpha, main="Learning rate curve of alpha", ylab = "Learning rate", 
     xlab = "Timepoint", cex.lab=2, cex.main=2, cex.axis=2)
plot(lrate_beta, main="Learning rate curve of beta", ylab = "Learning rate", 
     xlab = "Timepoint", cex.lab=2, cex.main=2, cex.axis=2)
plot(lrate_sigma, main="Learning rate curve of psi", ylab = "Learning rate", 
     xlab = "Timepoint", cex.lab=2, cex.main=2, cex.axis=2)

#Plot estimators
par(mfrow=c(1,1)) 
ggplot(data.frame(MLE=MLE_results_2019[["par"]],SGD=c(alpha[-1], beta,psi), parameters=c(rep("alpha",11),rep("beta",12),"psi")), aes(x=MLE, y=SGD)) + 
  geom_point( aes(shape=parameters, color=parameters), size=5) + geom_smooth(method=lm, se=TRUE,color = "black",size=0.5) + 
  labs(title="Comparison between estimators from MLE and SGD for the season 2018-2019") +
  scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) + scale_shape_manual(values=c(17, 16, 18)) +
  theme(legend.position = "bottom", text=element_text(size=17), legend.title=element_text(size=22) ) 


################################################################################################
########################## Add abilities on Last day ############################################
######################################################################################################
Last_Day_2019 <- as.Date(tail(time_point_all_2019,1)[1,1][[1]])

for(player in player_stat_2019$player){
  #if player's last game is not on final
  player_last_day = as.Date(tail(abilityList_2019[player][[1]][[1]],1))
  if(player_last_day != Last_Day_2019) {
    #retrive ability 
    ability = tail(abilityList_2019[player][[1]][[2]],1)
    #add ability on last day
    abilityList_2019[player][[1]] <- rbind(abilityList_2019[player][[1]], data.frame("Date"=Last_Day_2019,"Abilities"=ability))
  }
}
data.frame(df_2019[1:20,],df_2018[1:20,])

################################################################################################
########################## Player abilities trend ############################################
######################################################################################################
#c(seq(2,0.2,length.out=nrow(bind_rows(abilityList[df[1:20,1]]))))   family="serif"
ggplot(bind_rows(abilityList_2019[df_2019[1:20,1]], .id="df"), aes(Date, Abilities, group=df, colour=df )) +
  geom_line(size = .6, alpha = 1) + 
  theme_minimal() + ylim(0.2, 2.4) + 
  theme(legend.position = "bottom",text=element_text(size=17),legend.title=element_text(size=26)) + 
  #xlim(as.Date(c("2018-10-27","2019-07-31"))) +
  labs(title = "Abilities trends of 20 players with top comprehensive abilities for the NBA season 2018-2019.") + 
  geom_vline(aes(xintercept = as.Date("2019-04-13")), linetype = "longdash", alpha=.7) + 
  scale_linetype_manual(values=c("longdash","dotdash")) +
  scale_x_date(date_labels = ("%b %Y"), date_breaks='1 month', limits = as.Date(c("2018-10-16","2019-07-18"))) +
  scale_color_discrete(name="Players") +
  #geom_dl(aes(label=df), method=list(dl.combine("angled.boxes")))+ 
  geom_label_repel(aes(label = if_else(Date == max(Date), as.character(df), NA_character_)),
                   nudge_x = 30, direction = "y", na.rm = FALSE,segment.size = 0.2,hjust = 0.5, size = 5,
                   arrow = arrow(length = unit(0.01, "npc"), type = "closed", ends = "first"),force = 3,)

########10##############3
ggplot(bind_rows(abilityList_2019[df_2019[c(1:10,46),1]], .id="df"), aes(Date, Abilities, group=df, colour=df )) +
  geom_line(size = .6, alpha = 1) + 
  theme_minimal() + ylim(0.6, 2.2) + 
  theme(legend.position = "bottom",text=element_text(size=17),legend.title=element_text(size=26)) + 
  #xlim(as.Date(c("2018-10-16","2019-07-31"))) +
  labs(title = "Abilities trends of 10 players with top comprehensive abilities for the NBA season 2018-2019.") + 
  geom_vline(aes(xintercept = as.Date("2019-04-13")), linetype = "longdash", alpha=.7) + 
  scale_linetype_manual(values=c("longdash","dotdash")) +
  scale_x_date(date_labels = ("%b %Y"), date_breaks='1 month', limits = as.Date(c("2018-10-16","2019-07-18"))) +
  scale_color_discrete(name="Players") +
  geom_dl(aes(label=df), method=list(dl.combine("angled.boxes")))+ #+  scale_linetype_manual(values=rep(1,10))  #values = c(1,5,6,1,5,6,1,5,6,1))
  geom_label_repel(aes(label = if_else(Date == max(Date), as.character(df), NA_character_)),
                   nudge_x = 30, direction = "y", na.rm = FALSE,segment.size = 0.2,hjust = 0.5, size = 5,
                   arrow = arrow(length = unit(0.01, "npc"), type = "closed", ends = "first"),force = 3,)


###### Compare abilities from last year ################################################
diff_1819 <- merge(df_2018, df_2019, by.x="Player_Name", by.y = "Player_Name")
diff_1819$Improvement <- diff_1819$Player_Abilities.y - diff_1819$Player_Abilities.x
names(diff_1819) <- c("Player_Name", "Abilities_2018", "Abilities_2019", "Improvement")
diff_1819 <- diff_1819[order(-diff_1819$Improvement),]
#grid.table(head(diff_1819,20))


######################################################################################################
################## Compare PER and SGD rankings ########################################################
############################################################################################################
SplitName <- function(name){
  split = strsplit(as.character(name), "\\\\")[[1]][1]
  return(split)
}
# Read in player statistics from basketball reference
PER_2019 <- read.csv("Dropbox/Yuliang Weng - Dynamic Factor Models/data/PER_2019.csv",header=TRUE)
nrow(PER_2019)
# Select only player names and PER 
PER_2019 <- PER_2019[,c("Player","PER")]
# Correct players names
PER_2019$Player <- lapply(PER_2019[,1], FUN = SplitName)
PER_2019$Player <- as.character(PER_2019$Player)
###########
PER_2019_rank <- PER_2019[order(-PER_2019$PER),]
PER_2019_rank$PER_Rank <- 1:nrow(PER_2019)

# Get SGD Rankings 
df_2019_rank <- df_2019
df_2019_rank$SGD_Rank <- 1:nrow(df_2019)
df_2019_rank$Player_Name <- as.character(df_2019_rank$Player_Name)

# Merge two rankings 
player_Ranking_2019 <- merge(df_2019_rank, PER_2019_rank, by.x = "Player_Name", by.y = "Player")
player_Ranking_2019 <- player_Ranking_2019 %>% dplyr::select("Player_Name", "SGD_Rank", "PER_Rank")

######### Reorder rankings ################
rankdata_2019 <-  data.frame(Player=player_Ranking_2019$Player_Name, SGD=player_Ranking_2019$SGD_Rank, PER=player_Ranking_2019$PER_Rank)
rankdata_2019 <- rankdata_2019[order(rankdata_2019$SGD),]
rankdata_2019$SGD <- 1:nrow(rankdata_2019)

############## Rank Correlation ####################################
wilcox <- wilcox.test(rankdata_2019$SGD, rankdata_2019$PER,paired = TRUE)
kendall <- cor.test(rankdata_2019$SGD, rankdata_2019$PER, method = c("kendall"))
spearman <- cor.test(rankdata_2019$SGD, rankdata_2019$PER, method = c("spearman"))

corr_test <- data.frame("Correlation_Method" = c("Wilcox","Kendall","Spearman"), 
           "Sample_estimates" = c("", kendall$estimate[[1]], spearman$estimate[[1]]),
           "P_value" = c(wilcox$p.value, kendall$p.value, spearman$p.value),
           "Alternative_hypothesis" = c("true location shift is not equal to 0",
                                        "true tau is not equal to 0", 
                                        "true rho is not equal to 0"),
           "Statstics" = c(wilcox$statistic, kendall$statistic[[1]], spearman$statistic[[1]]))

grid.table(corr_test)


length(which(rankdata_2019$PER > rankdata_2019$SGD + tolerance | rankdata_2019$SGD > rankdata_2019$PER + tolerance))

###### Scatter plot for rankings ##################
tolerance <- 30
tolerance_label <- 100

ggplot(rankdata_2019, aes(x=SGD, y=PER, label=Player)) + 
  geom_point(color = 
      dplyr::case_when( 
          rankdata_2019$PER > rankdata_2019$SGD + tolerance | rankdata_2019$SGD > rankdata_2019$PER + tolerance ~ "red",
          rankdata_2019$PER < rankdata_2019$SGD + tolerance | rankdata_2019$SGD < rankdata_2019$PER + tolerance ~ "forestgreen",
          rankdata_2019$PER > rankdata_2019$SGD + tolerance_label | rankdata_2019$SGD > rankdata_2019$PER + tolerance_label ~ "blue",
          TRUE ~ "yellow"), alpha=.7, size=1.3)  + 
  labs(title = "The scatter plot of players rankings in Player Efficiency Rating and in estimations from SGD (for the NBA season 2018-19)") + 
  labs(subtitle = paste0("The blue solid line has the equation PER = SGD, the black dash line shown above and below of the solid line has the equation PER = 30 + SGD and PER = -30 + SGD respectively. 
  A player marked green indicates that the difference between his SGD ranking and PER ranking is less than ", tolerance, ". ", "If the player marked red then indicates that the difference is more than ", tolerance, "." ) ) +
  theme(legend.position = "bottom", text=element_text(size=17),legend.title=element_text(size = 28),plot.subtitle=element_text(size=13) ) +
  scale_x_continuous(position = "top") + scale_y_reverse() + xlab("SGD ranking") + ylab("PER ranking") + 
  geom_abline(slope = -1, intercept = tolerance, linetype="longdash") + 
  geom_abline(slope = -1, intercept = -tolerance,linetype="longdash") +
  geom_abline(slope = -1, intercept = 0, colour="deepskyblue2", size=1.2) +
  #label players with PER > SGD + 50
  geom_label_repel(data = subset(rankdata_2019, rankdata_2019$PER > rankdata_2019$SGD + tolerance_label), 
                   aes(label = Player), nudge_x = -20, direction = "y", na.rm = FALSE, segment.size = 0.3,
                   hjust = 0.5, force = 4, size = 5) +
  #label players with SGD > PER + 50
  geom_label_repel(data = subset(rankdata_2019, rankdata_2019$SGD > rankdata_2019$PER + tolerance_label), 
                   aes(label = Player), nudge_x = -20, direction = "y", na.rm = FALSE, segment.size = 0.3,
                   hjust = 0.5, force = 4, size = 5) +
  # label top players
  geom_label_repel(data = subset(rankdata_2019, rankdata_2019$SGD <= 5 | rankdata_2019$PER <= 5), 
                   aes(label = Player), nudge_y = 10, direction = "x", na.rm = FALSE, segment.size = 0.3,
                   hjust = 0.7, force = 10, size = 5) 



###########################################################################################
######################### MLE learning rate ##############################################
########################################################################################

MLE_learningRate <- function(parameters, data1,results){
  ## 2018
  result_MLE <- results
  
  len_sum = 12
  # Initialisation of parameters and learning rate
  alpha = matrix(rep(0, len_sum), nrow = len_sum, ncol = 1)
  beta = matrix(rep(1, len_sum), nrow = len_sum, ncol = 1)
  sigmaS = 2
  psi = log(sigmaS)
  
  alpha_lRate_0 = parameters[1]
  beta_lRate_0 = parameters[2]
  sigmaS_lRate_0 = parameters[3]
  
  # Facor a and c for updating learning rates
  lrate_a = parameters[4]
  lrate_c = parameters[5]#1
  
  j <- 0
  # For every time point, iteratively update parameters alpha, beta and sigma^2, and update player's abilities
  for(tp in data1){
    
    alpha_lRate = 1
    beta_lRate = 1
    sigmaS_lRate = 1
    # Update learning rate
    if(j == 0){
      alpha_lRate = alpha_lRate_0
      beta_lRate = beta_lRate_0
      sigmaS_lRate = sigmaS_lRate_0
    } else{
      alpha_lRate = alpha_lRate_0 * (1 + lrate_a * alpha_lRate_0 * j)^(-lrate_c)
      lrate_alpha <- c(lrate_alpha,alpha_lRate)
      beta_lRate = beta_lRate_0 * (1 + lrate_a * beta_lRate_0 * j)^(-lrate_c)
      sigmaS_lRate = sigmaS_lRate_0 * (1 + lrate_a * sigmaS_lRate_0 * j)^(-lrate_c)
      lrate_sigma <- c(lrate_sigma, sigmaS_lRate)
    }
    
    j <- j+1
    # Retrieve columns relevant to selected summaries
    data_mat <- data.matrix(tp[,3:(2+len_sum)])
    
    len_tp <- nrow(tp)
    ########################### Update parameters ###################################
    
    # Variance-covariance matrix & its Inverse
    covar <- (beta %*% t(beta))
    diag(covar) <- diag(covar) + exp(psi)
    icovar <- 0
    
    ivc_try <- try(solve(covar))
    
    if("try-error" %in% class(ivc_try)){
      return(-Inf)
    } else{
      icovar = ivc_try
    }
    
    ######## Update ALPHA #########
    
    sum_matrix_a <- icovar %*% t(sweep(data_mat, 2, alpha, "-"))   # X-alpha
    alpha = alpha + (alpha_lRate * rowSums(sum_matrix_a))
    alpha[1] <- 0
    ########   Update BETA  #########
    x_a <- sweep(data_mat, 2, alpha, "-")     ## Subtract alpha from summaries
    beta_part1 <- -len_tp * beta / (sum(beta^2) + exp(psi))  ## First part of beta's derivative
    beta_part2 <- matrix(c(rep(0,len_sum)), nrow=len_sum, ncol=1)      ## Initialise a empty vector at time t to store beta
    
    for(i in 1:len_sum){
      # Calculate the derivative of var-covar matrix w.r.t each beta
      cov_deri <- matrix(0L, nrow = len_sum, ncol = len_sum)
      cov_deri[i,] <- cov_deri[i,] + beta
      cov_deri[,i] <- cov_deri[,i] + beta
      # Part 2 of beta's derivative
      beta_ts <- 0.5 * ( x_a %*% icovar %*% cov_deri %*% icovar %*% t(x_a) )
      beta_part2[i] <- sum(diag(beta_ts))
    }
    beta_deri <- beta_part1 + beta_part2
    beta = beta + beta_lRate * beta_deri
    
    ########   Update SIGMA  #########
    sig_part1 <- -len_tp * (len_sum-1) / (2 * exp(psi))
    sig_part2 <- -0.5 * len_tp / (sum(beta^2) + exp(psi))
    sig_part3 <- 0.5 * sum(diag(x_a %*% icovar %*% icovar %*% t(x_a)))
    
    sigmaS_t <- exp(psi) * (sig_part1 + sig_part2 + sig_part3)    
    psi = psi + sigmaS_lRate*sigmaS_t
    
  }
  optimal_result <- sum((alpha[2:(len_sum)] - result_MLE[1:(len_sum-1)])^2) + sum((beta - result_MLE[len_sum:(2*len_sum-1)])^2) + sum((psi - result_MLE[2*len_sum])^2)
  return(optimal_result)
}

