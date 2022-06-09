# Run susie.
n <- 800
ldfile <- "small_data_11_sim_gaussian_pve_n_8_get_sumstats_n_1.ld_sample_n_file.in_n.ld"
# ldfile <- "small_data_11.ld_refout_file.refout.ld"
R <- as.matrix(fread(ldfile))
fit <- susie_rss(z,R,n = 800,min_abs_corr = 0.1,refine = FALSE,
                 verbose = TRUE)

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
