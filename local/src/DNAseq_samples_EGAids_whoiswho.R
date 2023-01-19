library(tidyverse)
library(readxl)
library(openxlsx)

samples_f <- snakemake@input[["dataset"]]
buoni_f <- snakemake@input[["pdo"]]
table_dna <- snakemake@output[["res"]]

## targeted
#samples_f <- "/scratch/trcanmed/biobanca/local/share/data/DNAseq_samples_EGAids.xlsx"
samples <- read_excel("/scratch/trcanmed/biobanca/local/share/data/DNAseq_samples_EGAids.xlsx", sheet = "targeted")
samples$type <- substr(samples$genalogy, 8, 10)
xh <- samples %>% filter(!type == "LMO")
xh$validation_status <- NA
xh$type <- NULL
o <- samples %>% filter(type == "LMO")

#buoni_f <- "/scratch/trcanmed/biobanca/local/share/data/whoiswho_validation_xen.tsv"
buoni <- read.table(buoni_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(buoni) <- c("smodel", "validation_status")
merged <- merge(o, buoni, by = "smodel")
## si perde il 2041
merged$validation_status.x <- NULL
merged$type <- NULL
names(merged)[names(merged)=="validation_status.y"] <- "validation_status"

samplestot <- rbind(xh, merged)

## mettere gli lmo derivati da lmh con una colonna "note" -> derived from human tissue
lmh <- c("CRC1337", "CRC1342", "CRC1563", "CRC1566")
lmh <- samplestot %>% filter(smodel %in% lmh)
lmh$type <- substr(lmh$genalogy, 8, 10)
lmh <- lmh %>% filter(type == "LMO")
lmh <- lmh$genalogy

samplestot <- as.data.frame(samplestot)
samplestot$note <- ""
rownames(samplestot) <- samplestot$genalogy
samplestot["CRC1337LMO0A04008001D02000", "note"] <- "derived from lmh"
samplestot["CRC1342LMO0A04007001D02000", "note"] <- "derived from lmh"
samplestot["CRC1563LMO0A04008001D02000", "note"] <- "derived from lmh"
samplestot["CRC1566LMO0A02005003D01000", "note"] <- "derived from lmh"

samplestot <- samplestot[order(samplestot$smodel),]
names(samplestot)[names(samplestot) == "genalogy"] <- "genealogy"

targeted <- samplestot
row_names_remove <- c("CRC0282-01-A-1")
targeted <- targeted[!(row.names(targeted) %in% row_names_remove),]
### shallowseq

samples <- read_excel("/scratch/trcanmed/biobanca/local/share/data/DNAseq_samples_EGAids.xlsx", sheet = "shallowseq")
samples$type <- substr(samples$genalogy, 8, 10)
samples <- samples %>% filter(!type == " CL")
samples <- samples %>% filter(!type == "LM ")
xh <- samples %>% filter(!type == "LMO")
xh$validation_status <- NA
xh$type <- NULL
o <- samples %>% filter(type == "LMO")

#buoni_f <- "/scratch/trcanmed/biobanca/local/share/data/whoiswho_validation_xen.tsv"
buoni <- read.table(buoni_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(buoni) <- c("smodel", "validation_status")
merged <- merge(o, buoni, by = "smodel")
## si perde il 2041
merged$validation_status.x <- NULL
merged$type <- NULL
names(merged)[names(merged)=="validation_status.y"] <- "validation_status"

samplestot <- rbind(xh, merged)

## mettere gli lmo derivati da lmh con una colonna "note" -> derived from human tissue
lmh <- c("CRC1337", "CRC1342", "CRC1563", "CRC1566")
lmh <- samplestot %>% filter(smodel %in% lmh)
lmh$type <- substr(lmh$genalogy, 8, 10)
lmh <- lmh %>% filter(type == "LMO")
lmh <- lmh$genalogy

samplestot <- as.data.frame(samplestot)
samplestot$note <- ""
rownames(samplestot) <- samplestot$genalogy
samplestot["CRC1337LMO0A04008001D02000", "note"] <- "derived from lmh"
samplestot["CRC1342LMO0A04007001D02000", "note"] <- "derived from lmh"
samplestot["CRC1563LMO0A04008001D02000", "note"] <- "derived from lmh"
samplestot["CRC1566LMO0A04010001D01000", "note"] <- "derived from lmh"

samplestot <- samplestot[order(samplestot$smodel),]
names(samplestot)[names(samplestot) == "genalogy"] <- "genealogy"

shallowseq <- samplestot

# Create a blank workbook
OUT <- createWorkbook()

# Add some sheets to the workbook
addWorksheet(OUT, "targeted")
addWorksheet(OUT, "shallowseq")

# Write the data to the sheets
writeData(OUT, sheet = "targeted", x = targeted)
writeData(OUT, sheet = "shallowseq", x = shallowseq)

# Export the file
saveWorkbook(OUT, table_dna)



