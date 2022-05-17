library(pedigreemm)
library(MASS)
library(Rcpp)

data_manipulation <- function(data) {
    id <- c(as.character(data$ID),
            as.character(data$DamID),
            as.character(data$SireID))
    id <- unique(id)
    id <- id[id != "0"]
    new_data <- data.frame("ID" = id)
    new_data <- merge(new_data, data, by = "ID", all.x = TRUE)
    new_data <- new_data[, c(1,2,3,5)]
    new_data <- new_data[order(new_data$Byear),]
    n <- nrow(new_data)
    new_data$newID <- 1:n
    data_tmp <- new_data[, c(1, 5)]
    colnames(data_tmp) <- c("SireID", "newSireID")
    new_data <- merge(new_data, data_tmp, by = "SireID", all.x = TRUE)
    data_tmp <- new_data[, c(2, 5)]
    colnames(data_tmp) <- c("DamID", "newDamID")
    new_data <- merge(new_data, data_tmp, by = "DamID", all.x = TRUE)
    new_data <- new_data[order(new_data$Byear),]
    new_data <- new_data[, c(-1, -2, -4)]
    new_data$newSireID[new_data$newSireID > new_data$newID &
                        !is.na(new_data$newSireID)] <- NA
    new_data$newDamID[new_data$newDamID > new_data$newID &
                        !is.na(new_data$newDamID)] <- NA
    return(new_data)
}
add_information_columns <- function(data, new_data) {
    new_data$y1 <- merge(new_data, data[data$ID %in% new_data$ID, ][, c("Y1", "ID")], by="ID")$Y1
    new_data$y2 <- merge(new_data, data[data$ID %in% new_data$ID, ][, c("Y2", "ID")], by="ID")$Y2
    new_data$sex <- merge(new_data, data[data$ID %in% new_data$ID, ][, c("Sex", "ID")], by="ID")$Sex
    new_data$byear <- merge(new_data, data[data$ID %in% new_data$ID, ][, c("Byear", "ID")], by="ID")$Byear
    new_dataSNP <- merge(new_data, data[data$ID %in% new_data$ID, -c(2,3,4,5,6,7)], by = "ID")
    new_dataSNP <- new_dataSNP[order(new_dataSNP$newID), ]
    return(new_dataSNP)
}

prepare_data <- function(data_path, size) {
    data <- read.table(data_path, header = TRUE, sep = ";")
    new_data <- data_manipulation(data)
    new_dataSNP <- add_information_columns(data, new_data)
    new_dataSNP <- new_dataSNP[0:size, ]
    return(new_dataSNP)
}

get_freq_A <- function(newData) {
    first_marker <- which(colnames(newData) == "SNP1")
    data_snp <- newData[first_marker:ncol(newData)]
    freq_A <- numeric(ncol(data_snp))
    for (i in seq_len(ncol(data_snp))) {
        n <- length(data_snp[, i]) * 2
        A <- 0
        for (j in data_snp[, i]) {
            if (j == -1) {
                A <- A + 2
            }
            else if (j == 0) {
                A <- A + 1
            }
        }
        freq_A[i] <- A/n
    }
    return(freq_A)
}

remove_low_MAF <- function(newData, freq_A) {
    first_marker <- which(colnames(newData) == "SNP1")
    for (i in first_marker:ncol(newData)){
        if (freq_A[i-first_marker + 1] < 0.05 | 1-freq_A[i-first_marker + 1] < 0.05) {
            newData <- newData[, -i]
        }
    }
    return(newData)
}

mme = function(y, X, Z, A, sigma_a, sigma_e) {
    cat("#Starting mme\n")
    alpha = sigma_e / sigma_a
    print(alpha)
    cat("#Calculating ginvA\n")
    invA = ginv(A)
    cat("#Calculationg C\n")
    C = rbind(cbind(t(X)%*%X, t(X)%*%Z),
              cbind(t(Z)%*%X, t(Z)%*%Z+invA*c(alpha)))
    cat("#Calculating rhs\n")
    rhs = rbind(t(X)%*%y, t(Z)%*%y)
    cat("#invC\n")
    invC = ginv(C)
    cat("#estimators\n")
    estimators = invC%*%rhs
    list(C = C, est = estimators, invA = invA, invC = invC)
}

EM = function(y, X, Z, A, sigma_a, sigma_e, output) {
  n = nrow(X)
  p = ncol(X) 
  q = nrow(A) 
  t = 1 #iteration number 1
  tmp = 0.1 #test for convergance
  while (tmp > 0.01) {
    cat("Loop :")
    cat(t)
    cat("\n")
    cat("##Start mme\n")
    mme_new = mme(y, X, Z, A, sigma_a, sigma_e)
    cat("##Ginv C\n")
    C_new = mme_new$invC
    Ck = C_new[(p+1):(p+q), (p+1):(p+q)]
    mme2 = mme_new$est
    a = as.matrix(mme2[(p+1):(p+q)])
    invA = mme_new$invA
    cat("##Sigma a new\n")
    sigma_a_new = (t(a)%*%invA%*%a + sum(diag(invA%*%Ck))*c(sigma_e))/q
    res = as.matrix(y-X%*%as.matrix(mme2[1:p]) - Z%*%as.matrix(mme2[(p+1):(p+q)]))
    X.tmp1 = cbind(X,Z) %*% C_new
    X.tmp2 = t(cbind(X,Z))
    cat("##Sigma e new\n")
    sigma_e_new = (t(res)%*%res + sum(diag(X.tmp1%*%X.tmp2))*c(sigma_e))/n
    tmp = max(abs(sigma_a - sigma_a_new), abs(sigma_e - sigma_e_new))
    sigma_a = sigma_a_new
    sigma_e = sigma_e_new
    write.csv(c(sigma_a,sigma_e,t, tmp), output, col.names = FALSE, row.names = FALSE, append = TRUE, quote = FALSE)
    t = t + 1
  }
  list(t = t, sigma_a = sigma_a, sigma_e = sigma_e)
}

mme2 = function(y, X, Z1, Z2, A, G, sigma_a, sigma_g, sigma_e) {
    alpha1 = sigma_e / sigma_a
    alpha2 = sigma_e / sigma_g
    invA = ginv(A)
    invG = ginv(G)
    C = rbind(cbind(t(X)%*%X, t(X)%*%Z1, t(X)%*%Z2),
              cbind(t(Z1)%*%X, t(Z1)%*%Z1+invA*c(alpha1), t(Z1)%*%Z2),
              cbind(t(Z2)%*%X, t(Z2)%*%Z1, t(Z2)%*%Z2 + invG*c(alpha2)))
    rhs = rbind(t(X)%*%y, t(Z1)%*%y, t(Z2)%*%y)
    invC = ginv(C)
    estimators = invC%*%rhs
    list(C = C, est = estimators)
}

ssSNPBLUP <- function(data_path, size = 10000, phenotype = c("y1", "y2")) {
    cat("###Preparing file...\n")
    data <- prepare_data(data_path, size)
    cat("###Cleaning markers...\n")
    freq_A <- get_freq_A(data)
    newData <- remove_low_MAF(data, freq_A)
    cat("###Calculating pedigree matrix...\n")
    ped <- pedigree(sire=newData$newSireID, dam=newData$newDamID, label=newData$newID)
    A <- as.matrix(getA(ped))
    cat("###Setting variables...\n")
    first_marker <- which(colnames(newData) == "SNP1")
    snp_data <- newData[first_marker:ncol(newData)]
    print(dim(snp_data))
    if (phenotype == "y1"){
        y <- newData$y1
    }
    else if(phenotype == "y2"){
        y <- newData$y2
    }
    else{
        stop("Wrong phenotype! Try y1 or y2.")
    }
    byear <- newData$byear
    sex <- newData$sex
    X <- matrix(0,size,3)
    X[, 1] <- sex
    X[, 2] <- -sex
    X[, 3] <- byear
    Z <- diag(size)
    p2 <- get_freq_A(newData)
    cat("###Estimating sigma_a and sigma_e...\n")
    sigma_file_name <- paste("output/Sigma_estiamtes", phenotype, ".csv", sep = "")
    estimate_sigma = EM(y, X, Z, A, sigma_a = 0.3*var(y), sigma_e = 0.7*var(y), sigma_file_name)
    cat("###Setting new variables to mme2...\n")
    Z2 <- data.matrix(snp_data)
    sigma_a = estimate_sigma$sigma_a
    sigma_e = estimate_sigma$sigma_e
    sigma_g = sigma_a / length(p2)
    G = diag(length(p2))
    cat("###Calculating mme2...\n")
    results = mme2(y, X, Z, Z2, A, G, sigma_a, sigma_g, sigma_e)
    C_file_name = paste("output/C", phenotype, ".csv", sep = "")
    est_file_name = paste("output/estimates", phenotype, ".csv", sep = "")
    print(dim(results$C))
    write.csv(results$C, C_file_name, col.names = FALSE, row.names = FALSE, append = TRUE, quote = FALSE)
    write.csv(results$est, est_file_name, col.names = FALSE, row.names = FALSE, append = TRUE, quote = FALSE)
    list(C = results$C, est = results$est, sigma_a = sigma_a, sigme_e = sigma_e, sigma_g = sigma_g)
}

Wald_test <- function(est, Z, A, sigma_a, sigma_e, size) {
    G = A*c(sigma_a)
    R = diag(size)*c(sigma_e)
    V = Z%*%G%*%t(Z) + R
    varB = ginv(t(X)%*%ginv(V)%*%X)
    seB = sqrt(diag(varB))
    testWalda = est[1:3] / seB
    p_value = 2*pnorm(abs(testWalda), lower.tail = FALSE)
    return(p_value)
}

estimate_acc <- function(y, X, Z, A, sigma_a, sigma_e, C, size) {
    invC = ginv(C)
    invC22 = invC[3:size+3, 3:size+3]
    alpha = sigma_e / sigma_a
    r2 = diag(1 - invC22*c(alpha))
    r = sqrt(r2)
    list(r = r, r2 = r2)
}

snp_effect <- function(estimates, size, size_perm_effects) {
    n = ncol(estimates)
    est_snp = estimates[size_perm_effects + size : n]
    cat("SNPs length: \n")
    cat(length((est_snp)))
    s = sd(est_snp)
    m = mean(est_snp)
    est_snp = est_snp / s
    cat("SD for SNPs:\n")
    cat(sd(est_snp))
    W = est_snp
    p.value = 2*pnorm(abs(W), lower.tail = FALSE)
    snp_s <- which(p.value < 0.05)
    list(snp = snp_s, m = m, s = s)
}


#########################################################
dir.create("output", showWarnings = FALSE)
sourceCpp("MatrixInverse.cpp")
#First phenotype
y1_results <- ssSNPBLUP("data.csv", size = 1000, phenotype = "y1")
#Wald Test
w_test <- Wald_test(y1_results$est)
y2_results <- ssSNPBLUP("data.csv", size = 1000, phenotype = "y2")