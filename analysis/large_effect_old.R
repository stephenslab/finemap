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
