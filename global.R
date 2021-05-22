library(survival)
library(survminer)
library(riskRegression)
library(rms)
library(caret)
library(glmnet)
library(BeSS)
library(shiny)
library(shinyWidgets)
library(readr)
library(rintrojs)
library(Cairo)
library(ggrepel)
####################################################################
stomach.surv <- read_delim("STAD/STAD_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
liver.surv <- read_delim("LIHC/LIHC_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
breast.surv <- read_delim("BRCA/BRCA_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
lung.surv <- read_delim("LUAD/LUAD_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
brain.surv <- read_delim("GBM/GBM_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
kidney.surv <- read_delim("KIRC/KIRC_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
aculeuk.surv <- read_delim("LAML/LAML_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
adreno.surv <- read_delim("ACC/ACC_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
bladder.surv <- read_delim("BLCA/BLCA_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
brain.lowgrad.surv <- read_delim("LGG/LGG_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
cervix.surv <- read_delim("CESC/CESC_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
chol.surv <- read_delim("CHOL/CHOL_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
colon.surv <- read_delim("COAD/COAD_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
esophagus.surv <- read_delim("ESCA/ESCA_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
headneck.surv <- read_delim("HNSC/HNSC_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
kdchrom.surv <- read_delim("KICH/KICH_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
kdrnlp.surv <- read_delim("KIRP/KIRP_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
lngsq.surv <- read_delim("LUSC/LUSC_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
blymph.surv <- read_delim("DLBC/DLBC_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
mesothel.surv <- read_delim("MESO/MESO_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
ovary.surv <- read_delim("OV/OV_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
pancr.surv <- read_delim("PAAD/PAAD_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
prost.surv <- read_delim("PRAD/PRAD_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
paragang.surv <- read_delim("PCPG/PCPG_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rect.surv <- read_delim("READ/READ_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
sarc.surv <- read_delim("SARC/SARC_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
skin.surv <- read_delim("SKCM/SKCM_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
tstc.surv <- read_delim("TGCT/TGCT_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
thymus.surv <- read_delim("THYM/THYM_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
thyroid.surv <- read_delim("THCA/THCA_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
uterus.surv <- read_delim("UCS/UCS_survival.txt" , "\t", escape_double = FALSE, trim_ws = TRUE)
utcendo.surv <- read_delim("UCEC/UCEC_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
uva.surv <- read_delim("UVM/UVM_survival.txt", "\t", escape_double = FALSE, trim_ws = TRUE)


datasets <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", 
              "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", 
              "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

#################################################################
set.surv <- function(dataset) {
  if (dataset == "STAD") {
    surv <- stomach.surv
  } else if (dataset == "GBM") {
    surv <- brain.surv
  } else if (dataset == "LIHC") {
    surv <- liver.surv
  } else if (dataset == "LUAD") {
    surv <- lung.surv
  } else if (dataset == "ACC") {
    surv <- adreno.surv
  } else if (dataset == "BLCA") {
    surv <- bladder.surv
  } else if (dataset == "KIRC") {
    surv <- kidney.surv
  } else if (dataset == "BRCA") {
    surv <- breast.surv
  } else if (dataset == "CESC") {
    surv <- cervix.surv
  } else if (dataset == "CHOL") {
    surv <- chol.surv
  } else if (dataset == "COAD") {
    surv <- colon.surv
  } else if (dataset == "DLBC") {
    surv <- blymph.surv
  } else if (dataset == "ESCA") {
    surv <- esophagus.surv
  } else if (dataset == "HNSC") {
    surv <- headneck.surv
  } else if (dataset == "KICH") {
    surv <- kdchrom.surv 
  } else if (dataset == "KIRP") {
    surv <- kdrnlp.surv
  } else if (dataset == "LAML") {
    surv <- aculeuk.surv
  } else if (dataset == "LGG") {
    surv <- brain.lowgrad.surv
  } else if (dataset == "LUSC") {
    surv <- lngsq.surv
    colnames(surv)[1] <- gsub("^xena_", "", colnames(surv)[1])
  } else if (dataset == "MESO") {
    surv <- mesothel.surv
  } else if (dataset == "OV") {
    surv <- ovary.surv
  } else if (dataset == "PAAD") {
    surv <- pancr.surv
  } else if (dataset == "PCPG") {
    surv <- paragang.surv
  } else if (dataset == "PRAD") {
    surv <- prost.surv
  } else if (dataset == "READ") {
    surv <- rect.surv
  } else if (dataset == "SARC") {
    surv <- sarc.surv
  } else if (dataset == "SKCM") {
    surv <- skin.surv
  } else if (dataset == "TGCT") {
    surv <- tstc.surv
  } else if (dataset == "THCA") {
    surv <- thyroid.surv
  } else if (dataset == "THYM") {
    surv <- thymus.surv
  } else if (dataset == "UCEC") {
    surv <- utcendo.surv
  } else if (dataset == "UCS") {
    surv <- uterus.surv
  } else if (dataset == "UVM") {
    surv <- uva.surv
  }
  return(surv)
}

######################################################

sub.pancan <- function(pancan.file, dataset) {
  subpan <- pancan.file[which(pancan.file$`cancer type abbreviation` == dataset & (pancan.file$sample_type == "Primary Tumor" | 
                                pancan.file$sample_type == "Primary Blood Derived Cancer - Peripheral Blood")), ]
  return(subpan)
}

set.genes <- function(subpan) {
  genes <- colnames(subpan[3:(ncol(subpan)-2)])
  return(genes)
}

set.exp <- function(subpan, genes) {
  xena.rpkm.exp <- subpan[, c("sample", genes)]
  return(xena.rpkm.exp)
}

set.foldchange <- function(subpan, genes, fold.cols) {
   xena.rpkm.fold <- merge(subpan, fold.cols, by = "sample")
   xena.rpkm.fold <- xena.rpkm.fold[, c(1, (ncol(xena.rpkm.fold) - (ncol(fold.cols)-2)):ncol(xena.rpkm.fold))]
   colnames(xena.rpkm.fold)[2:ncol(xena.rpkm.fold)] <- genes
   return(xena.rpkm.fold)
}

glmnet.os <- function(glm.surv.data, nfold, alph) {
  glm.surv.data <- na.omit(glm.surv.data)
  glm.surv.data$OS.time <- glm.surv.data$OS.time + 0.00001
  surv <- Surv(glm.surv.data$OS.time, glm.surv.data$OS)
  if (nfold == "sample_size") {
    cv.glmfit <- cv.glmnet(as.matrix(glm.surv.data[4:ncol(glm.surv.data)]), surv, alpha = alph, family = "cox", nfolds = nrow(glm.surv.data))
  } else {
    set.seed(123)
    folds <- createFolds(surv, 10, list = TRUE, returnTrain = FALSE)
    foldids <- rep(1, length(surv))
    foldids[folds$Fold02] <- 2
    foldids[folds$Fold03] <- 3
    foldids[folds$Fold04] <- 4
    foldids[folds$Fold05] <- 5
    foldids[folds$Fold06] <- 6
    foldids[folds$Fold07] <- 7
    foldids[folds$Fold08] <- 8
    foldids[folds$Fold09] <- 9
    foldids[folds$Fold10] <- 10
    cv.glmfit <- cv.glmnet(as.matrix(glm.surv.data[4:ncol(glm.surv.data)]), surv, foldid = foldids, alpha = alph, family = "cox")
  }
  coefs <- as.matrix(coef(cv.glmfit))
  coefs <- cbind(rownames(coefs), coefs)
  colnames(coefs) <- c("genes", "coeffs")
  coefs <- coefs[coefs[, 2] != 0]
  ncoef <- length(coefs)/2
  coef.vec <- as.numeric(coefs[(ncoef + 1):length(coefs)])
  names(coef.vec) <- coefs[1:ncoef]
  prognostic.index <- function(glm.surv.data, coef.vec) {
    prog.ind <- numeric(nrow(glm.surv.data))
    for (i in 1:length(coef.vec)) {
      prog.ind <- prog.ind + glm.surv.data[names(coef.vec)[i]]*coef.vec[i]
    }
    names(prog.ind) <- "PI"
    return(prog.ind)
  }
  if (ncoef == 0) {
    glm.data.km <- NULL
  } else {
    glm.data.km <- glm.surv.data %>% mutate(prog_ind = prognostic.index(glm.surv.data, coef.vec),
                                            Risk = ifelse(prog_ind > median(prognostic.index(glm.surv.data, coef.vec)$PI), "High", "Low"))
  }  
  return(list(coef.vec, cv.glmfit, glm.data.km))
}

glmnet.dss <- function(glm.surv.data, nfold, alph) {
  glm.surv.data <- na.omit(glm.surv.data)
  glm.surv.data$DSS.time <- glm.surv.data$DSS.time + 0.00001
  surv <- Surv(glm.surv.data$DSS.time, glm.surv.data$DSS)
  if (nfold == "sample_size") {
    cv.glmfit <- cv.glmnet(as.matrix(glm.surv.data[4:ncol(glm.surv.data)]), surv, alpha = alph, family = "cox", nfolds = nrow(glm.surv.data))
  } else {
    set.seed(123)
    folds <- createFolds(surv, 10, list = TRUE, returnTrain = FALSE)
    foldids <- rep(1, length(surv))
    foldids[folds$Fold02] <- 2
    foldids[folds$Fold03] <- 3
    foldids[folds$Fold04] <- 4
    foldids[folds$Fold05] <- 5
    foldids[folds$Fold06] <- 6
    foldids[folds$Fold07] <- 7
    foldids[folds$Fold08] <- 8
    foldids[folds$Fold09] <- 9
    foldids[folds$Fold10] <- 10
    cv.glmfit <- cv.glmnet(as.matrix(glm.surv.data[4:ncol(glm.surv.data)]), surv, foldid = foldids, alpha = alph, family = "cox")
  }
  coefs <- as.matrix(coef(cv.glmfit))
  coefs <- cbind(rownames(coefs), coefs)
  colnames(coefs) <- c("genes", "coeffs")
  coefs <- coefs[coefs[, 2] != 0]
  ncoef <- length(coefs)/2
  coef.vec <- as.numeric(coefs[(ncoef + 1):length(coefs)])
  names(coef.vec) <- coefs[1:ncoef]
  prognostic.index <- function(glm.surv.data, coef.vec) {
    prog.ind <- numeric(nrow(glm.surv.data))
    for (i in 1:length(coef.vec)) {
      prog.ind <- prog.ind + glm.surv.data[names(coef.vec)[i]]*coef.vec[i]
    }
    names(prog.ind) <- "PI"
    return(prog.ind)
  }
  if (ncoef == 0) {
    glm.data.km <- NULL
  } else {
    glm.data.km <- glm.surv.data %>% mutate(prog_ind = prognostic.index(glm.surv.data, coef.vec),
                                            Risk = ifelse(prog_ind > median(prognostic.index(glm.surv.data, coef.vec)$PI), "High", "Low"))
  }  
  return(list(coef.vec, cv.glmfit, glm.data.km))
}

glmnet.dfi <- function(glm.surv.data, nfold, alph) {
  glm.surv.data <- na.omit(glm.surv.data)
  glm.surv.data$DFI.time <- glm.surv.data$DFI.time + 0.00001
  surv <- Surv(glm.surv.data$DFI.time, glm.surv.data$DFI)
  if (nfold == "sample_size") {
    cv.glmfit <- cv.glmnet(as.matrix(glm.surv.data[4:ncol(glm.surv.data)]), surv, alpha = alph, family = "cox", nfolds = nrow(glm.surv.data))
  } else {
    set.seed(123)
    folds <- createFolds(surv, 10, list = TRUE, returnTrain = FALSE)
    foldids <- rep(1, length(surv))
    foldids[folds$Fold02] <- 2
    foldids[folds$Fold03] <- 3
    foldids[folds$Fold04] <- 4
    foldids[folds$Fold05] <- 5
    foldids[folds$Fold06] <- 6
    foldids[folds$Fold07] <- 7
    foldids[folds$Fold08] <- 8
    foldids[folds$Fold09] <- 9
    foldids[folds$Fold10] <- 10
    cv.glmfit <- cv.glmnet(as.matrix(glm.surv.data[4:ncol(glm.surv.data)]), surv, foldid = foldids, alpha = alph, family = "cox")
  }  
  coefs <- as.matrix(coef(cv.glmfit))
  coefs <- cbind(rownames(coefs), coefs)
  colnames(coefs) <- c("genes", "coeffs")
  coefs <- coefs[coefs[, 2] != 0]
  ncoef <- length(coefs)/2
  coef.vec <- as.numeric(coefs[(ncoef + 1):length(coefs)])
  names(coef.vec) <- coefs[1:ncoef]
  prognostic.index <- function(glm.surv.data, coef.vec) {
    prog.ind <- numeric(nrow(glm.surv.data))
    for (i in 1:length(coef.vec)) {
      prog.ind <- prog.ind + glm.surv.data[names(coef.vec)[i]]*coef.vec[i]
    }
    names(prog.ind) <- "PI"
    return(prog.ind)
  }
  if (ncoef == 0) {
    glm.data.km <- NULL
  } else {
    glm.data.km <- glm.surv.data %>% mutate(prog_ind = prognostic.index(glm.surv.data, coef.vec),
                                            Risk = ifelse(prog_ind > median(prognostic.index(glm.surv.data, coef.vec)$PI), "High", "Low"))
  }  
  return(list(coef.vec, cv.glmfit, glm.data.km))
}

glmnet.pfi <- function(glm.surv.data, nfold, alph) {
  glm.surv.data <- na.omit(glm.surv.data)
  glm.surv.data$PFI.time <- glm.surv.data$PFI.time + 0.00001
  surv <- Surv(glm.surv.data$PFI.time, glm.surv.data$PFI)
  if (nfold == "sample_size") {
    cv.glmfit <- cv.glmnet(as.matrix(glm.surv.data[4:ncol(glm.surv.data)]), surv, alpha = alph, family = "cox", nfolds = nrow(glm.surv.data))
  } else {
    set.seed(123)
    folds <- createFolds(surv, 10, list = TRUE, returnTrain = FALSE)
    foldids <- rep(1, length(surv))
    foldids[folds$Fold02] <- 2
    foldids[folds$Fold03] <- 3
    foldids[folds$Fold04] <- 4
    foldids[folds$Fold05] <- 5
    foldids[folds$Fold06] <- 6
    foldids[folds$Fold07] <- 7
    foldids[folds$Fold08] <- 8
    foldids[folds$Fold09] <- 9
    foldids[folds$Fold10] <- 10
    cv.glmfit <- cv.glmnet(as.matrix(glm.surv.data[4:ncol(glm.surv.data)]), surv, foldid = foldids, alpha = alph, family = "cox")
  }
  coefs <- as.matrix(coef(cv.glmfit))
  coefs <- cbind(rownames(coefs), coefs)
  colnames(coefs) <- c("genes", "coeffs")
  coefs <- coefs[coefs[, 2] != 0]
  ncoef <- length(coefs)/2
  coef.vec <- as.numeric(coefs[(ncoef + 1):length(coefs)])
  names(coef.vec) <- coefs[1:ncoef]
  prognostic.index <- function(glm.surv.data, coef.vec) {
    prog.ind <- numeric(nrow(glm.surv.data))
    for (i in 1:length(coef.vec)) {
      prog.ind <- prog.ind + glm.surv.data[names(coef.vec)[i]]*coef.vec[i]
    }
    names(prog.ind) <- "PI"
    return(prog.ind)
  }
  if (ncoef == 0) {
    glm.data.km <- NULL
  } else {
    glm.data.km <- glm.surv.data %>% mutate(prog_ind = prognostic.index(glm.surv.data, coef.vec),
                                            Risk = ifelse(prog_ind > median(prognostic.index(glm.surv.data, coef.vec)$PI), "High", "Low"))
  }
  return(list(coef.vec, cv.glmfit, glm.data.km))
}

#####################################################
bess.os <- function(bess.result, bess.data) {
  data <- na.omit(bess.data)
  df.os <- as.data.frame(bess.result$bestmodel[[1]])
  rownames(df.os) <- gsub("xbest", "", rownames(df.os))
  df.os <- cbind(rownames(df.os), df.os)
  colnames(df.os) <- c("genes", "coeffs")
  rownames(df.os) <- NULL
  prognostic.index <- function(data, df.os) {
    prog.ind <- numeric(nrow(data))
    for (i in 1:nrow(df.os)) {
      prog.ind <- prog.ind + data[df.os$genes[i]]*df.os$coeffs[i]
    }
    names(prog.ind) <- "PI"
    return(prog.ind)
  }
  if (nrow(df.os) == 0) {
    data.km <- NULL
  } else {
    data.km <- data %>% mutate(prog_ind = prognostic.index(data, df.os),
                               Risk = ifelse(prog_ind > median(prognostic.index(data, df.os)$PI), "High", "Low"))
  }
  return(list(df.os, data.km))
}

bess.dss <- function(bess.result, bess.data) {
  data <- na.omit(bess.data)
  df.dss <- as.data.frame(bess.result$bestmodel[[1]])
  rownames(df.dss) <- gsub("xbest", "", rownames(df.dss))
  df.dss <- cbind(rownames(df.dss), df.dss)
  colnames(df.dss) <- c("genes", "coeffs")
  rownames(df.dss) <- NULL
  prognostic.index <- function(data, df.dss) {
    prog.ind <- numeric(nrow(data))
    for (i in 1:nrow(df.dss)) {
      prog.ind <- prog.ind + data[df.dss$genes[i]]*df.dss$coeffs[i]
    }
    names(prog.ind) <- "PI"
    return(prog.ind)
  }
  if (nrow(df.dss) == 0) {
    data.km <- NULL
  } else {
    data.km <- data %>% mutate(prog_ind = prognostic.index(data, df.dss),
                               Risk = ifelse(prog_ind > median(prognostic.index(data, df.dss)$PI), "High", "Low"))
  }
  return(list(df.dss, data.km))
}

bess.dfi <- function(bess.result, bess.data) {
  data <- na.omit(bess.data)
  df.dfi <- as.data.frame(bess.result$bestmodel[[1]])
  rownames(df.dfi) <- gsub("xbest", "", rownames(df.dfi))
  df.dfi <- cbind(rownames(df.dfi), df.dfi)
  colnames(df.dfi) <- c("genes", "coeffs")
  rownames(df.dfi) <- NULL
  prognostic.index <- function(data, df.dfi) {
    prog.ind <- numeric(nrow(data))
    for (i in 1:nrow(df.dfi)) {
      prog.ind <- prog.ind + data[df.dfi$genes[i]]*df.dfi$coeffs[i]
    }
    names(prog.ind) <- "PI"
    return(prog.ind)
  }
  if (nrow(df.dfi) == 0) {
    data.km <- NULL
  } else {
    data.km <- data %>% mutate(prog_ind = prognostic.index(data, df.dfi),
                               Risk = ifelse(prog_ind > median(prognostic.index(data, df.dfi)$PI), "High", "Low"))
  }
  return(list(df.dfi, data.km))
}

bess.pfi <- function(bess.result, bess.data) {
  data <- na.omit(bess.data)
  df.pfi <- as.data.frame(bess.result$bestmodel[[1]])
  rownames(df.pfi) <- gsub("xbest", "", rownames(df.pfi))
  df.pfi <- cbind(rownames(df.pfi), df.pfi)
  colnames(df.pfi) <- c("genes", "coeffs")
  rownames(df.pfi) <- NULL
  prognostic.index <- function(data, df.pfi) {
    prog.ind <- numeric(nrow(data))
    for (i in 1:nrow(df.pfi)) {
      prog.ind <- prog.ind + data[df.pfi$genes[i]]*df.pfi$coeffs[i]
    }
    names(prog.ind) <- "PI"
    return(prog.ind)
  }
  if (nrow(df.pfi) == 0) {
    data.km <- NULL
  } else {
    data.km <- data %>% mutate(prog_ind = prognostic.index(data, df.pfi),
                               Risk = ifelse(prog_ind > median(prognostic.index(data, df.pfi)$PI), "High", "Low"))
  }
  return(list(df.pfi, data.km))
}

##################################
multi.cox.model.os <- function(surv.data, exprs) {
  colnames(exprs) <- gsub("-", "_", colnames(exprs))
  exprs.surv.data <- merge(surv.data, exprs)
  exprs.surv.data <- exprs.surv.data[, c(1, 3:4, 12:ncol(exprs.surv.data))]
  exprs.covars <- paste(colnames(exprs)[2:ncol(exprs)], collapse = " + ")
  exprs.surv.form <- as.formula(paste("Surv(OS.time, OS)", "~", "(", exprs.covars, ")"))
  return(coxph(exprs.surv.form, data = exprs.surv.data, model = TRUE, x = TRUE, y = TRUE))
}
########################################################
multi.cox.model.dss <- function(surv.data, exprs) {
  colnames(exprs) <- gsub("-", "_", colnames(exprs))
  exprs.surv.data <- merge(surv.data, exprs)
  exprs.surv.data <- exprs.surv.data[, c(1, 5:6, 12:ncol(exprs.surv.data))]
  exprs.covars <- paste(colnames(exprs)[2:ncol(exprs)], collapse = " + ")
  exprs.surv.form <- as.formula(paste("Surv(DSS.time, DSS)", "~", "(", exprs.covars, ")"))
  return(coxph(exprs.surv.form, data = exprs.surv.data, model = TRUE, x = TRUE, y = TRUE))
}

##############################################################
multi.cox.model.dfi <- function(surv.data, exprs) {
  colnames(exprs) <- gsub("-", "_", colnames(exprs))
  exprs.surv.data <- merge(surv.data, exprs)
  exprs.surv.data <- exprs.surv.data[, c(1, 7:8, 12:ncol(exprs.surv.data))]
  exprs.covars <- paste(colnames(exprs)[2:ncol(exprs)], collapse = " + ")
  exprs.surv.form <- as.formula(paste("Surv(DFI.time, DFI)", "~", "(", exprs.covars, ")"))
  return(coxph(exprs.surv.form, data = exprs.surv.data, model = TRUE, x = TRUE, y = TRUE))
}

############################################################
multi.cox.model.pfi <- function(surv.data, exprs) {
  colnames(exprs) <- gsub("-", "_", colnames(exprs))
  exprs.surv.data <- merge(surv.data, exprs)
  exprs.surv.data <- exprs.surv.data[, c(1, 9:10, 12:ncol(exprs.surv.data))]
  exprs.covars <- paste(colnames(exprs)[2:ncol(exprs)], collapse = " + ")
  exprs.surv.form <- as.formula(paste("Surv(PFI.time, PFI)", "~", "(", exprs.covars, ")"))
  return(coxph(exprs.surv.form, data = exprs.surv.data, model = TRUE, x = TRUE, y = TRUE))
}
################################
km.os <- function(cox.result, cox.data) {
  cox.data <- na.omit(cox.data)
  cox.os <- summary(cox.result)$coefficients
  cox.os <- cox.os[cox.os[, 5] < 0.05, , drop = FALSE]
  prognostic.index <- function(data, cox.os) {
    prog.ind <- numeric(nrow(data))
    if (nrow(cox.os) == 1) {
      prog.ind <- prog.ind + data[rownames(cox.os)]*cox.os[,1] 
    } else  {
      for (i in 1:nrow(cox.os)) {
        prog.ind <- prog.ind + data[names(cox.os[,1])[i]]*cox.os[,1][i]
      }
    }
    names(prog.ind) <- "PI"
    return(prog.ind)
  }
  if (nrow(cox.os) == 0) {
    data.km <- NULL
  } else {
    data.km <- cox.data %>% mutate(prog_ind = prognostic.index(cox.data, cox.os),
                               Risk = ifelse(prog_ind > median(prognostic.index(cox.data, cox.os)$PI), "High", "Low"))
  }
  return(list(cox.os, data.km))
}

km.dss <- function(cox.result, cox.data) {
  cox.data <- na.omit(cox.data)
  cox.dss <- summary(cox.result)$coefficients
  cox.dss <- cox.dss[cox.dss[, 5] < 0.05, , drop = FALSE]
  prognostic.index <- function(data, cox.dss) {
    prog.ind <- numeric(nrow(data))
    if (nrow(cox.dss) == 1) {
      prog.ind <- prog.ind + data[rownames(cox.dss)]*cox.dss[,1] 
    } else  {
      for (i in 1:nrow(cox.dss)) {
        prog.ind <- prog.ind + data[names(cox.dss[,1])[i]]*cox.dss[,1][i]
      }
    }
    names(prog.ind) <- "PI"
    return(prog.ind)
  }
  if (nrow(cox.dss) == 0) {
    data.km <- NULL
  } else {
    data.km <- cox.data %>% mutate(prog_ind = prognostic.index(cox.data, cox.dss),
                                   Risk = ifelse(prog_ind > median(prognostic.index(cox.data, cox.dss)$PI), "High", "Low"))
  }
  return(list(cox.dss, data.km))
}

km.dfi <- function(cox.result, cox.data) {
  cox.data <- na.omit(cox.data)
  cox.dfi <- summary(cox.result)$coefficients
  cox.dfi <- cox.dfi[cox.dfi[, 5] < 0.05, , drop = FALSE]
  prognostic.index <- function(data, cox.dfi) {
    prog.ind <- numeric(nrow(data))
    if (nrow(cox.dfi) == 1) {
      prog.ind <- prog.ind + data[rownames(cox.dfi)]*cox.dfi[,1] 
    } else  {
      for (i in 1:nrow(cox.dfi)) {
        prog.ind <- prog.ind + data[names(cox.dfi[,1])[i]]*cox.dfi[,1][i]
      }
    }
    names(prog.ind) <- "PI"
    return(prog.ind)
  }
  if (nrow(cox.dfi) == 0) {
    data.km <- NULL
  } else {
    data.km <- cox.data %>% mutate(prog_ind = prognostic.index(cox.data, cox.dfi),
                                   Risk = ifelse(prog_ind > median(prognostic.index(cox.data, cox.dfi)$PI), "High", "Low"))
  }
  return(list(cox.dfi, data.km))
}

km.pfi <- function(cox.result, cox.data) {
  cox.data <- na.omit(cox.data)
  cox.pfi <- summary(cox.result)$coefficients
  cox.pfi <- cox.pfi[cox.pfi[, 5] < 0.05, , drop = FALSE]
  prognostic.index <- function(data, cox.pfi) {
    prog.ind <- numeric(nrow(data))
    if (nrow(cox.pfi) == 1) {
      prog.ind <- prog.ind + data[rownames(cox.pfi)]*cox.pfi[,1] 
    } else  {
      for (i in 1:nrow(cox.pfi)) {
        prog.ind <- prog.ind + data[names(cox.pfi[,1])[i]]*cox.pfi[,1][i]
      }
    }
    names(prog.ind) <- "PI"
    return(prog.ind)
  }
  if (nrow(cox.pfi) == 0) {
    data.km <- NULL
  } else {
    data.km <- cox.data %>% mutate(prog_ind = prognostic.index(cox.data, cox.pfi),
                                   Risk = ifelse(prog_ind > median(prognostic.index(cox.data, cox.pfi)$PI), "High", "Low"))
  }
  return(list(cox.pfi, data.km))
}









