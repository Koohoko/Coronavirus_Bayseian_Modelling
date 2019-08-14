library(Biostrings)
library(tidyverse)
library(ggplot2)
library(SynMut)

cds_omsn <- readDNAStringSet("../data/cleaned_cds_osmn.fasta")
meta_omsn <- read_csv("../data/meta_data_osmn.csv")

# codon usage -------------------------------------------------------------

rscu_cds <- get_rscu(cds_omsn)
rscu_cds <- as_tibble(rscu_cds)
index_tmp <- !names(rscu_cds) %in% c("ATG", "TGG", "TAA", "TAG", "TGA")
rscu_cds <- rscu_cds[,index_tmp]  # remove stop codons, atg and tgg

## prepare data
df_tmp <- meta_omsn %>% select(id, protein_name, `Virus Genus`, Host)
df_tmp <- bind_cols(df_tmp, rscu_cds)
df_tmp <- df_tmp %>% gather(key = "codon", value = "RSCU", AAA:TTT)
df_tmp$AA <- as.character(translate(DNAStringSet(df_tmp$codon)))

## codon usage between virus genus
ggplot(df_tmp) + 
    geom_boxplot(aes(x = codon, y = RSCU, fill = AA, color = AA),
                 outlier.size = 0.1, outlier.alpha = 0.2)+
    facet_grid(vars(`Virus Genus`),vars(AA),scales = "free")+
    theme(axis.text.x = element_text(angle = 90))


# jags model --------------------------------------------------------------

library(rjags)
df_jags <- subset(df_tmp, codon == "TTT")
df_jags$group_values <- as.numeric(factor(df_jags$`Virus Genus`))

mod_string = " model {
for (i in 1:length(RSCU)) {
  RSCU[i] ~ dgamma(alpha[group_values[i]], beta[group_values[i]])
}

for (j in 1:max(group_values)) {
  alpha[j] ~ dnorm(theta1, 1/sigma1)
  beta[j] ~ dnorm(theta2, 1/sigma2)
}

theta1 ~ dnorm(1.0, 1/100)
theta2 ~ dnorm(1.0, 1/100)

sigma1 ~ dexp(1.0)
sigma2 ~ dexp(1.0)

} "

set.seed(2019)

data_jags = as.list(df_jags)

params = c("alpha", "beta", "theta1", "theta2", "sigma1", "sigma2")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e5)
Sys.time()
mod_csim = as.mcmc(do.call(rbind, mod_sim))

## convergence diagnostics
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)

## compute DIC
# dic = dic.samples(mod, n.iter=1e3)

## model checking
(pm_params = colMeans(mod_csim))

alpha_s <- pm_params[1:4]
beta_s <- pm_params[5:8]
mean_s <- pm_params[1:4]/pm_params[5:8]
yhat = mean_s[data_jags$group_values]
resid = data_jags$RSCU - yhat
plot(resid)

plot(jitter(yhat), resid)

#### we can check the posterior probability that the RSCU values for the 
# Alphacoronavirus is greater than the RSCU values for the Betacoronavirus
n_sim <- nrow(mod_csim)
set.seed(2019)
alpha_rscu_pred <- rgamma(n_sim,mod_csim[,1], mod_csim[,5])
beta_rscu_pred <- rgamma(n_sim,mod_csim[,2], mod_csim[,6])

mean(alpha_rscu_pred > beta_rscu_pred)

# 
# n_sim <- nrow(mod_csim)
# rscu__pred <- rgamma(n_sim, mod_csim[,1:4], mod_csim[,5:8])



# render rmd --------------------------------------------------------------
rmarkdown::render("./index.rmd")
