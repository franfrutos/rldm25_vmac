# Load packages ----
if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(
  dplyr,
  Rmisc
)

source("Scripts/functions.R")

dStan <- list()
dPlots <- list()

# Load VMAC ----
d <- read.csv("Input/le2023_e2.csv") %>%
  filter(phase == 2) %>%
  mutate(ID = as.integer(as.factor(ID))) %>%
  mutate(reward = trialPay/500)

subj <- d$ID
j <- length(unique(subj)) # 40 participants
d$ID <- rep(1:j, each = 384)

dat <-  list(
  N = nrow(d),
  J = j,
  subj = subj,
  reward = d$reward,
  Singleton = d$distractType,
  valid =  ifelse(d$timeout != 1 & d$trialPropGoodSamples > .25, 1, 0),
  trial = (d$trial - mean(d$trial)) / (max(d$trial) - min(d$trial))
)

dat[["N_valid"]] <- sum(dat[["valid"]])
dat[["Omission"]] <- d$omissionTrial[which(dat[["valid"]] == 1)]

dStan[["VMAC"]] <- dat
dPlots[["VMAC"]] <- d %>%
  mutate(Epoch = create_epochs(trial, 32)) %>%
  filter(timeout != 1 & trialPropGoodSamples > .25) %>%
  dplyr::summarise(P = mean(omissionTrial), 
                   .by = c(ID, Epoch, distractType))%>%
  summarySEwithin(measurevar = "P", withinvars = c("Epoch", "distractType"),
                  idvar = "ID") %>%
  dplyr::rename("Singleton" = distractType) %>%
  mutate(Singleton = ifelse(Singleton == 1, "High", "Low"),
         Data = "Le et al. (2024)")
# Load UMAC ----
d <- read.csv("Input/df_chow.csv") %>%
  filter(phase == 2) %>%
  mutate(ID = as.integer(as.factor(ID))) %>%
  mutate(reward = trialPay/40)

subj <- d$ID
j <- length(unique(subj)) # 40 participants
d$ID <- rep(1:j, each = 512)

dat <-  list(
  N = nrow(d),
  J = j,
  subj = subj,
  reward = d$reward,
  Singleton = d$distractType,
  valid =  ifelse(d$searchTimeout != 1 & d$trialPropGoodSamples > .25, 1, 0),
  trial = (d$trial - mean(d$trial)) / (max(d$trial) - min(d$trial))
)

dat[["N_valid"]] <- sum(dat[["valid"]])
dat[["Omission"]] <- d$omission1Trial[which(dat[["valid"]] == 1)]

dStan[["UMAC"]] <- dat
dPlots[["UMAC"]] <- d %>%
  mutate(Epoch = create_epochs(trial, 32)) %>%
  filter(searchTimeout != 1 & trialPropGoodSamples > .25) %>%
  dplyr::summarise(P = mean(omission1Trial), 
            .by = c(ID, Epoch, distractType))%>%
  summarySEwithin(measurevar = "P", withinvars = c("Epoch", "distractType"),
                  idvar = "ID") %>%
  dplyr::rename("Singleton" = distractType) %>%
  mutate(Singleton = ifelse(Singleton == 2, "High", "Low"),
         Data = "Chow et al. (2024)")

rm(d, subj, j)

