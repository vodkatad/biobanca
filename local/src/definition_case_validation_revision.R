t2 <- read.table("/scratch/trcanmed/biobanca/local/share/data/TABLE_2_SImo_201222.tsv",
                 quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# Quindi, se vi va facciamo in parallelo un controllo con Table 2 alla mano: una volta filtrati via i QC "Not passed" in colonna B, rimangono un tot di 242 casi di partenza.
# Sempre a partire da questi 242 casi, per fare un check che vi tornino i casi dei vari sottogruppi, vi dico i filtri che ho applicato nella Table 2 per individuare i vari sottogruppi che metterei nel circles di Fig 1E:
# - i Validation successful: in verde (per continuità con Fig 1D), tot 132 casi (tutti inclusi i 4 LMH) (li ottieni filtrando via dalla colonna E i "Failed", i "Not Performed" e i "(Blancs)" che sarebbero le celle rimaste bianche in corrispondenza dei not established)
# - Validation failed: in nero, tot 24 (sono i "failed" in colonna E)
# - Validation not performed: in grigio, include tutti gli established che non hanno fatto la validazione (sono i "not performed" in colonna E)
# - Not established: potremmo farli in blu (invece che grigio scuro come dicevo prima, sempre per coerenza con Fig 1 B e C), tot 57 (a ripartire dai 242 con QC passed, sono i "Failed" nella colonna D della Derivation)

t2 <- t2 %>% filter(!QC..passed.not.passed. == "Not passed" ) # 242 samples

v_suc <- t2 %>% filter(!VALIDATION..successful.failed.not.performed. == "Failed")
v_suc <- v_suc %>% filter(!VALIDATION..successful.failed.not.performed. == "Not performed")
v_suc$VALIDATION..successful.failed.not.performed.[v_suc$VALIDATION..successful.failed.not.performed.==""]<-"blank"
v_suc <- v_suc %>% filter(!VALIDATION..successful.failed.not.performed. == "blank")

## validation_successful 132 confirmed

v_failed <- t2 %>% filter(VALIDATION..successful.failed.not.performed. == "Failed")

## validation_failed 24 confirmed

v_notp <- t2 %>% filter(VALIDATION..successful.failed.not.performed. == "Not performed")

not_est <- t2 %>% filter(DERIVATION..successful.failed. == "Failed")

### not est confirmed

v_suc$type <- "Validation successful"
v_suc <- v_suc[,c(1,6)]
## aggiungo nuovamente il 1241
v_suc <- rbind(v_suc, c("CRC1241", "Validation successful"))
v_failed$type <- "Validation failed"
v_failed <- v_failed[,c(1,6)]
v_notp$type <- "Validation not performed"
v_notp <- v_notp[,c(1,6)]
not_est$type <- "Not established"
not_est <- not_est[,c(1,6)]

lmh <- t2 %>% filter(ORIGIN..patient.sample.PDX.fresh.frozen.. == "patient sample")
lmh2 <- t2 %>% filter(ORIGIN..patient.sample.PDX.fresh.frozen.. == "patient sample ")
lmh3 <- t2 %>% filter(ORIGIN..patient.sample.PDX.fresh.frozen.. == "PDX (fresh)+patient sample")

lmh <- rbind(lmh, lmh2, lmh3)

whoiswho <- as.data.frame(rbind(v_suc, v_failed, v_notp, not_est))

for (i in seq(rownames(whoiswho))) {
  if (whoiswho[i, "type"] == "Not established")  {
    whoiswho[i, "derivation_type"] <- "Derivation failed"
  } else {
    whoiswho[i, "derivation_type"] <- "Derivation successful"
  }
}

write.table(whoiswho, file = "/scratch/trcanmed/biobanca/local/share/data/whoiswho_validation_xen_revision_derivation.tsv", quote = FALSE,
            sep = "\t", col.names = TRUE, row.names = FALSE)

#write.table(lmh, file = "/scratch/trcanmed/biobanca/local/share/data/lmh_detele_multivariate.tsv", quote = FALSE,
#            sep = "\t", col.names = TRUE, row.names = FALSE)

#lmh <- lmh$CASE
#whoiswho <- whoiswho %>% filter(!CASE %in% lmh)

#write.table(whoiswho, file = "/scratch/trcanmed/biobanca/local/share/data/whoiswho_validation_xen_nolmh.tsv", quote = FALSE,
#            sep = "\t", col.names = TRUE, row.names = FALSE)
