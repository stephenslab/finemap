library(data.table)
library(susieR)
set.seed(1)

# Load the summary data.
dat1 <- readRDS("small_data_11.rds")
dat2 <- readRDS("small_data_11_sim_gaussian_pve_n_8.rds")
dat3 <- readRDS("small_data_11_sim_gaussian_pve_n_8_get_sumstats_n_1.rds")
b    <- drop(dat2$meta$true_coef)
maf  <- dat1$maf$in_sample
bhat <- dat3$sumstats$bhat
shat <- dat3$sumstats$shat
z    <- bhat/shat
cat("True causal SNPs:\n")
print(which(b != 0))

# Run susie.
n      <- 800
ldfile <- "small_data_11_sim_gaussian_pve_n_8_get_sumstats_n_1.ld_sample_n_file.in_n.ld"
R      <- as.matrix(fread(ldfile))
fit    <- susie_rss(z,R,n = 800,min_abs_corr = 0.1,refine = FALSE,
                    verbose = TRUE)
cat("SuSiE CSs:\n")
print(fit$sets[c("cs","purity")])

# Run FINEMAP.
p   <- length(b)
dat <- data.frame(rsid       = 1:p,
                  chromosome = rep(1,p),
                  position   = rep(1,p),
                  allele1    = rep("A",p),
                  allele2    = rep("C",p),
                  maf        = round(maf,digits = 6),
                  beta       = round(bhat,digits = 6),
                  se         = round(shat,digits = 6))
write.table(dat,"sim.z",quote = FALSE,col.names = TRUE,row.names = FALSE)
system(paste("./finemap_v1.4.1_x86_64 --sss --log --in-files sim.master",
             "--n-causal-snps 5"))
