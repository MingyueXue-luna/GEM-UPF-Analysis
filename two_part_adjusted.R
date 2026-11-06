# ============================================================
# Two-part GEE model for microbiome
# Adjusted for:
# Age + Sex + total_reads + Ethnicity + Income + Total.fiber.intake
# ============================================================

rm(list = ls())

setwd("C:/Users/Administrator/Desktop/Paper/NOVA")

# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------

library(geepack)
library(Matrix)
library(aod)
library(MESS)
library(doParallel)

# ------------------------------------------------------------
# Read data
# ------------------------------------------------------------

sheet.genus <- read.csv(
  "Nova_genus_interaction.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

str(sheet.genus[, 1:18])
colnames(sheet.genus)

# ------------------------------------------------------------
# Basic data preparation
# ------------------------------------------------------------

# family ID
sheet.genus$fam <- as.factor(
  gsub("-\\d{2}$", "", sheet.genus$SubjectID)
)

# variable type
sheet.genus$Age <- as.numeric(sheet.genus$Age)
sheet.genus$Sex <- as.factor(sheet.genus$Sex)
sheet.genus$Ethnicity <- as.factor(sheet.genus$Ethnicity)
sheet.genus$Income <- as.factor(sheet.genus$Income)

sheet.genus$total_reads <- as.numeric(sheet.genus$total_reads)
sheet.genus$Total.fiber.intake <- as.numeric(sheet.genus$Total.fiber.intake)
sheet.genus$kathy_kcal_4.0 <- as.numeric(sheet.genus$kathy_kcal_4.0)


# ============================================================
# QIC function
# ============================================================

QIC_TP <- function(logistic_obj, cont_obj) {
  
  if (
    !(
      "geeglm" %in% class(logistic_obj) &&
      "geeglm" %in% class(cont_obj)
    )
  ) {
    stop("QIC_TP requires geeglm objects as input")
  }
  
  # -------------------------
  # Logistic part
  # -------------------------
  
  mu <- logistic_obj$fitted.values
  y <- logistic_obj$y
  
  type <- family(logistic_obj)$family
  
  if (type != "binomial") {
    stop("Error: logistic part is not binomial")
  }
  
  # 避免 log(0)
  eps <- 1e-8
  mu <- pmin(pmax(mu, eps), 1 - eps)
  
  quasi_logistic_i <-
    y * log(mu / (1 - mu)) +
    log(1 - mu)
  
  pi_i <- predict(
    logistic_obj,
    type = "response"
  )
  
  # -------------------------
  # Continuous part
  # -------------------------
  
  mu1 <- cont_obj$fitted.values
  y1 <- cont_obj$y
  
  type1 <- family(cont_obj)$family
  
  quasi_cont_i0 <- switch(
    
    type1,
    
    poisson =
      (y1 * log(mu1)) - mu1,
    
    gaussian =
      ((y1 - mu1)^2) / -2,
    
    Gamma =
      -(y1 / mu1 - log(mu1)),
    
    stop("Error: distribution not recognized")
  )
  
  # -------------------------
  # Merge
  # -------------------------
  
  quasi_0 <- merge(
    data.frame(
      id = rownames(quasi_logistic_i),
      quasi_logistic = quasi_logistic_i
    ),
    data.frame(
      id = rownames(quasi_cont_i0),
      quasi_cont = quasi_cont_i0
    ),
    all.x = TRUE,
    sort = FALSE
  )
  
  quasi_0$quasi_cont[
    is.na(quasi_0$quasi_cont)
  ] <- 0
  
  quasi_0 <- merge(
    quasi_0,
    data.frame(
      id = rownames(pi_i),
      pi = pi_i
    ),
    sort = FALSE
  )
  
  # -------------------------
  # Calculations
  # -------------------------
  
  quasi <- sum(
    quasi_0$quasi_logistic * quasi_0$pi +
      quasi_0$quasi_cont * (1 - quasi_0$pi)
  )
  
  n <- nrow(quasi_0)
  
  params <-
    length(coef(logistic_obj)) +
    length(coef(cont_obj))
  
  kpm <-
    params +
    length(logistic_obj$geese$alpha) +
    length(cont_obj$geese$alpha)
  
  V <- as.matrix(
    bdiag(
      logistic_obj$geese$vbeta.naiv,
      cont_obj$geese$vbeta.naiv
    )
  )
  
  AIinverse <- solve(V)
  
  Vr <- as.matrix(
    bdiag(
      logistic_obj$geese$vbeta,
      cont_obj$geese$vbeta
    )
  )
  
  trace <- sum(
    diag(
      AIinverse %*% Vr
    )
  )
  
  QIC <- -2 * (quasi - trace)
  QICu <- -2 * (quasi - params)
  
  if ((n - kpm - 1) > 0) {
    
    QICC <-
      QIC +
      (2 * kpm * (kpm + 1)) /
      (n - kpm - 1)
    
  } else {
    
    QICC <- NA
  }
  
  output <- c(
    QIC,
    QICu,
    quasi,
    trace,
    params,
    QICC
  )
  
  names(output) <- c(
    "QIC",
    "QICu",
    "Quasi Lik",
    "CIC",
    "params",
    "QICC"
  )
  
  return(output)
}


# ============================================================
# Wald test
# ============================================================

wald.test.lowert <- function(
    Sigma,
    b,
    Terms = NULL,
    L = NULL,
    H0 = NULL,
    df = NULL,
    verbose = FALSE
) {
  
  if (is.null(Terms) & is.null(L)) {
    stop("One of Terms or L must be used.")
  }
  
  if (!is.null(Terms) & !is.null(L)) {
    stop("Only one of Terms or L must be used.")
  }
  
  if (is.null(Terms)) {
    
    w <- nrow(L)
    
    Terms <- seq(length(b))[
      colSums(L) > 0
    ]
    
  } else {
    
    w <- length(Terms)
  }
  
  if (is.null(H0)) {
    H0 <- rep(0, w)
  }
  
  if (w != length(H0)) {
    stop(
      "Vectors of tested coefficients and null hypothesis differ in length"
    )
  }
  
  if (is.null(L)) {
    
    L <- matrix(
      0,
      nrow = w,
      ncol = length(b)
    )
    
    for (i in 1:w) {
      L[i, Terms[i]] <- 1
    }
  }
  
  dimnames(L) <- list(
    paste0(
      "L",
      seq_len(NROW(L))
    ),
    names(b)
  )
  
  f <- L %*% b
  
  V <- Sigma
  
  mat <- qr.solve(
    L %*% V %*% t(L)
  )
  
  stat <-
    t(f - H0) %*%
    mat %*%
    (f - H0)
  
  p <- pchisq(
    stat,
    df = w,
    lower.tail = FALSE
  )
  
  if (is.null(df)) {
    
    res <- list(
      chi2 = c(
        chi2 = stat,
        df = w,
        P = p
      )
    )
    
  } else {
    
    fstat <- stat / nrow(L)
    
    df1 <- nrow(L)
    df2 <- df
    
    res <- list(
      
      chi2 = c(
        chi2 = stat,
        df = w,
        P = p
      ),
      
      Ftest = c(
        Fstat = fstat,
        df1 = df1,
        df2 = df2,
        P = 1 - pf(
          fstat,
          df1,
          df2
        )
      )
    )
  }
  
  structure(
    
    list(
      Sigma = Sigma,
      b = b,
      Terms = Terms,
      H0 = H0,
      L = L,
      result = res,
      verbose = verbose,
      df = df
    ),
    
    class = "wald.test"
  )
}


# ============================================================
# Main two-part function
# ============================================================

TP_v2 <- function(data) {
  
  # ----------------------------------------------------------
  # 1. Basic preparation
  # ----------------------------------------------------------
  
  data$taxa <- suppressWarnings(
    as.numeric(data$taxa)
  )
  
  # presence / absence
  data$Y_bin <- ifelse(
    is.na(data$taxa),
    NA,
    ifelse(
      data$taxa == 0,
      0,
      1
    )
  )
  
  # 不直接 log(0)
  data$ln_relab <- NA_real_
  
  idx_positive <-
    !is.na(data$taxa) &
    data$taxa > 0
  
  data$ln_relab[
    idx_positive
  ] <- log(
    data$taxa[
      idx_positive
    ]
  )
  
  # factors
  data$fam <- as.factor(data$fam)
  data$Sex <- as.factor(data$Sex)
  data$Ethnicity <- as.factor(data$Ethnicity)
  data$Income <- as.factor(data$Income)
  
  
  # ----------------------------------------------------------
  # 2. Binary data
  # ----------------------------------------------------------
  
  vars_bin <- c(
    "Y_bin",
    "kathy_kcal_4.0",
    "Age",
    "Sex",
    "total_reads",
    "Ethnicity",
    "Income",
    "Total.fiber.intake",
    "fam"
  )
  
  data_bin <- data[
    complete.cases(
      data[, vars_bin]
    ),
  ]
  
  # numeric finite check
  data_bin <- data_bin[
    is.finite(data_bin$kathy_kcal_4.0) &
      is.finite(data_bin$Age) &
      is.finite(data_bin$total_reads) &
      is.finite(data_bin$Total.fiber.intake),
  ]
  
  data_bin <- droplevels(
    data_bin
  )
  
  # family必须相邻
  data_bin <- data_bin[
    order(data_bin$fam),
  ]
  
  
  # ----------------------------------------------------------
  # 3. Continuous positive-abundance data
  # ----------------------------------------------------------
  
  vars_logn <- c(
    "ln_relab",
    "kathy_kcal_4.0",
    "Age",
    "Sex",
    "total_reads",
    "Ethnicity",
    "Income",
    "Total.fiber.intake",
    "fam"
  )
  
  data_logn <- data[
    !is.na(data$Y_bin) &
      data$Y_bin == 1,
  ]
  
  data_logn <- data_logn[
    complete.cases(
      data_logn[, vars_logn]
    ),
  ]
  
  data_logn <- data_logn[
    is.finite(data_logn$ln_relab) &
      is.finite(data_logn$kathy_kcal_4.0) &
      is.finite(data_logn$Age) &
      is.finite(data_logn$total_reads) &
      is.finite(data_logn$Total.fiber.intake),
  ]
  
  data_logn <- droplevels(
    data_logn
  )
  
  data_logn <- data_logn[
    order(data_logn$fam),
  ]
  
  
  # ----------------------------------------------------------
  # 4. Zero proportion
  # ----------------------------------------------------------
  
  if (nrow(data_bin) == 0) {
    
    return(
      list(
        zp = NA,
        R1 = rep(NA, 3),
        R2 = rep(NA, 3),
        subglobalpvalue0 = NA,
        subglobalpvalue1 = NA,
        globalpvalue = NA,
        R3 = c(
          LOGNORMAL = NA,
          TP = NA
        ),
        err_1 = NA,
        err_2 = NA
      )
    )
  }
  
  zp <- 1 - mean(
    data_bin$Y_bin
  )
  
  
  # ----------------------------------------------------------
  # 5. Formulas
  # ----------------------------------------------------------
  
  form.bin <- as.formula(
    paste(
      "Y_bin ~",
      "kathy_kcal_4.0 +",
      "Age + Sex + total_reads +",
      "Ethnicity + Income + Total.fiber.intake"
    )
  )
  
  form.logn <- as.formula(
    paste(
      "ln_relab ~",
      "kathy_kcal_4.0 +",
      "Age + Sex + total_reads +",
      "Ethnicity + Income + Total.fiber.intake"
    )
  )
  
  
  # ==========================================================
  # 6. Binary part
  # ==========================================================
  
  if (
    zp < 0.01 ||
    zp > 0.80
  ) {
    
    TPMa <- NULL
    
  } else {
    
    m_test <- tryCatch(
      
      glm(
        Y_bin ~ kathy_kcal_4.0,
        data = data_bin,
        family = binomial(
          link = "logit"
        )
      ),
      
      error = function(e) {
        NULL
      }
    )
    
    if (is.null(m_test)) {
      
      TPMa <- NULL
      
    } else {
      
      beta_test <- coef(m_test)[
        grepl(
          "kathy_kcal_4.0",
          names(coef(m_test))
        )
      ]
      
      if (
        length(beta_test) == 0 ||
        any(!is.finite(beta_test))
      ) {
        
        TPMa <- NULL
        
      } else {
        
        val <- max(
          abs(beta_test),
          na.rm = TRUE
        )
        
        if (val > 10) {
          
          TPMa <- NULL
          
        } else {
          
          TPMa <- tryCatch(
            
            geeglm(
              form.bin,
              data = data_bin,
              id = fam,
              corstr = "independence",
              family = binomial(
                link = "logit"
              ),
              std.err = "san.se"
            ),
            
            error = function(e) {
              NULL
            }
          )
        }
      }
    }
  }
  
  
  # ==========================================================
  # 7. Continuous log-normal part
  # ==========================================================
  
  if (
    nrow(data_logn) < 10
  ) {
    
    TPMb_0 <- NULL
    
  } else {
    
    # 关键：model.matrix 必须和实际模型一致
    modelmat <- tryCatch(
      
      model.matrix(
        ~ kathy_kcal_4.0 +
          Age +
          Sex +
          total_reads +
          Ethnicity +
          Income +
          Total.fiber.intake,
        data = data_logn
      ),
      
      error = function(e) {
        NULL
      }
    )
    
    if (
      is.null(modelmat)
    ) {
      
      TPMb_0 <- NULL
      
    } else if (
      anyNA(modelmat) ||
      any(!is.finite(modelmat))
    ) {
      
      TPMb_0 <- NULL
      
    } else {
      
      qqrr <- tryCatch(
        qr(modelmat),
        error = function(e) NULL
      )
      
      if (
        is.null(qqrr)
      ) {
        
        TPMb_0 <- NULL
        
      } else if (
        qqrr$rank < ncol(modelmat)
      ) {
        
        TPMb_0 <- NULL
        
      } else {
        
        TPMb_0 <- tryCatch(
          
          geeglm(
            form.logn,
            data = data_logn,
            id = fam,
            corstr = "exchangeable",
            family = gaussian(),
            std.err = "san.se"
          ),
          
          error = function(e) {
            NULL
          }
        )
      }
    }
  }
  
  
  # ==========================================================
  # 8. Extract results
  # ==========================================================
  
  if (
    !is.null(TPMa) &&
    !is.null(TPMb_0)
  ) {
    
    Wald0 <- summary(TPMa)$coefficients[
      ,
      c(
        "Estimate",
        "Std.err",
        "Pr(>|W|)"
      ),
      drop = FALSE
    ]
    
    Wald1 <- summary(TPMb_0)$coefficients[
      ,
      c(
        "Estimate",
        "Std.err",
        "Pr(>|W|)"
      ),
      drop = FALSE
    ]
    
    C0 <- TPMa$coefficients
    C1 <- TPMb_0$coefficients
    
    C <- c(
      C0,
      C1
    )
    
    V0 <- summary(TPMa)$cov.unscaled
    V1 <- summary(TPMb_0)$cov.unscaled
    
    V <- as.matrix(
      bdiag(
        V0,
        V1
      )
    )
    
    globalpvalue <- tryCatch(
      
      wald.test.lowert(
        b = C,
        Sigma = V,
        Terms = which(
          grepl(
            "kathy_kcal_4.0",
            names(C)
          )
        )
      )$result$chi2[3],
      
      error = function(e) {
        NA
      }
    )
    
    subglobalpvalue0 <- tryCatch(
      
      wald.test.lowert(
        b = C0,
        Sigma = as.matrix(
          bdiag(V0)
        ),
        Terms = which(
          grepl(
            "kathy_kcal_4.0",
            names(C0)
          )
        )
      )$result$chi2[3],
      
      error = function(e) {
        NA
      }
    )
    
    subglobalpvalue1 <- tryCatch(
      
      wald.test.lowert(
        b = C1,
        Sigma = as.matrix(
          bdiag(V1)
        ),
        Terms = which(
          grepl(
            "kathy_kcal_4.0",
            names(C1)
          )
        )
      )$result$chi2[3],
      
      error = function(e) {
        NA
      }
    )
    
    R1 <- Wald0[
      grepl(
        "kathy_kcal_4.0",
        rownames(Wald0)
      ),
      ,
      drop = FALSE
    ]
    
    R2 <- Wald1[
      grepl(
        "kathy_kcal_4.0",
        rownames(Wald1)
      ),
      ,
      drop = FALSE
    ]
    
    R1 <- as.numeric(
      R1[1, ]
    )
    
    R2 <- as.numeric(
      R2[1, ]
    )
    
    err_1 <- as.numeric(
      summary(TPMa)$error == 0
    )
    
    err_2 <- as.numeric(
      summary(TPMb_0)$error == 0
    )
    
    R3 <- tryCatch(
      
      c(
        LOGNORMAL =
          MESS:::QIC(
            TPMb_0
          )[1],
        
        TP =
          QIC_TP(
            TPMa,
            TPMb_0
          )[1]
      ),
      
      error = function(e) {
        
        c(
          LOGNORMAL = NA,
          TP = NA
        )
      }
    )
    
    
  } else if (
    !is.null(TPMa) &&
    is.null(TPMb_0)
  ) {
    
    Wald0 <- summary(TPMa)$coefficients[
      ,
      c(
        "Estimate",
        "Std.err",
        "Pr(>|W|)"
      ),
      drop = FALSE
    ]
    
    C0 <- TPMa$coefficients
    V0 <- summary(TPMa)$cov.unscaled
    
    subglobalpvalue0 <- tryCatch(
      
      wald.test.lowert(
        b = C0,
        Sigma = as.matrix(
          bdiag(V0)
        ),
        Terms = which(
          grepl(
            "kathy_kcal_4.0",
            names(C0)
          )
        )
      )$result$chi2[3],
      
      error = function(e) {
        NA
      }
    )
    
    subglobalpvalue1 <- NA
    
    R1 <- Wald0[
      grepl(
        "kathy_kcal_4.0",
        rownames(Wald0)
      ),
      ,
      drop = FALSE
    ]
    
    R1 <- as.numeric(
      R1[1, ]
    )
    
    R2 <- rep(
      NA,
      3
    )
    
    err_1 <- as.numeric(
      summary(TPMa)$error == 0
    )
    
    err_2 <- NA
    globalpvalue <- NA
    
    R3 <- c(
      LOGNORMAL = NA,
      TP = NA
    )
    
    
  } else if (
    is.null(TPMa) &&
    !is.null(TPMb_0)
  ) {
    
    Wald1 <- summary(TPMb_0)$coefficients[
      ,
      c(
        "Estimate",
        "Std.err",
        "Pr(>|W|)"
      ),
      drop = FALSE
    ]
    
    C1 <- TPMb_0$coefficients
    V1 <- summary(TPMb_0)$cov.unscaled
    
    subglobalpvalue0 <- NA
    
    subglobalpvalue1 <- tryCatch(
      
      wald.test.lowert(
        b = C1,
        Sigma = as.matrix(
          bdiag(V1)
        ),
        Terms = which(
          grepl(
            "kathy_kcal_4.0",
            names(C1)
          )
        )
      )$result$chi2[3],
      
      error = function(e) {
        NA
      }
    )
    
    R1 <- rep(
      NA,
      3
    )
    
    R2 <- Wald1[
      grepl(
        "kathy_kcal_4.0",
        rownames(Wald1)
      ),
      ,
      drop = FALSE
    ]
    
    R2 <- as.numeric(
      R2[1, ]
    )
    
    err_2 <- as.numeric(
      summary(TPMb_0)$error == 0
    )
    
    err_1 <- NA
    globalpvalue <- NA
    
    R3 <- tryCatch(
      
      c(
        LOGNORMAL =
          MESS:::QIC(
            TPMb_0
          )[1],
        
        TP = NA
      ),
      
      error = function(e) {
        
        c(
          LOGNORMAL = NA,
          TP = NA
        )
      }
    )
    
    
  } else {
    
    R1 <- rep(
      NA,
      3
    )
    
    R2 <- rep(
      NA,
      3
    )
    
    err_1 <- NA
    err_2 <- NA
    
    globalpvalue <- NA
    subglobalpvalue0 <- NA
    subglobalpvalue1 <- NA
    
    R3 <- c(
      LOGNORMAL = NA,
      TP = NA
    )
  }
  
  
  # ==========================================================
  # 9. Return fixed-length vector
  # ==========================================================
  
  r <- c(
    
    ZP = zp,
    
    coefbin = R1[1],
    std_bin = R1[2],
    p_bin = R1[3],
    
    coef_cont = R2[1],
    std_cont = R2[2],
    p_cont = R2[3],
    
    pvalue.bin =
      subglobalpvalue0,
    
    pvalue.cont =
      subglobalpvalue1,
    
    global.p =
      globalpvalue,
    
    QIC_LOGNORMAL =
      R3["LOGNORMAL"],
    
    QIC_TP =
      R3["TP"],
    
    err_1 =
      err_1,
    
    err_2 =
      err_2
  )
  
  return(
    as.numeric(r)
  )
}


# ============================================================
# Batch analysis
# ============================================================

# 你原来是 i = 1:102
n_taxa <- 103

# 14列结果
result_TP <- data.frame(
  matrix(
    NA_real_,
    nrow = n_taxa,
    ncol = 14
  )
)

# 保存 taxa 名称
taxa_names <- colnames(
  sheet.genus
)[17:(16 + n_taxa)]

for (i in seq_len(n_taxa)) {
  
  cat(
    "Running",
    i,
    "/",
    n_taxa,
    ":",
    taxa_names[i],
    "\n"
  )
  
  data <- data.frame(
    taxa = sheet.genus[[i + 17]],
    Age = sheet.genus$Age,
    Sex = sheet.genus$Sex,
    total_reads = sheet.genus$total_reads,
    Ethnicity = sheet.genus$Ethnicity,
    Income = sheet.genus$Income,
    Total.fiber.intake = sheet.genus$Total.fiber.intake,
    kathy_kcal_4.0 = sheet.genus$kathy_kcal_4.0,
    fam = sheet.genus$fam
  )
  
  # family连续排列
  data <- data[
    order(data$fam),
  ]
  
  tmp <- tryCatch(
    TP_v2(data),
    
    error = function(e) {
      
      message(
        "Taxa ",
        taxa_names[i],
        " failed: ",
        conditionMessage(e)
      )
      
      rep(
        NA_real_,
        14
      )
    }
  )
  
  result_TP[i, ] <- tmp
}


# ============================================================
# Result names
# ============================================================

names(result_TP) <- c(
  "ZP",
  "coefbin",
  "std_bin",
  "Pr_bin",
  "coef_cont",
  "std_cont",
  "Pr_cont",
  "pvalue.bin",
  "pvalue.cont",
  "global.p",
  "QIC_LOGNORMAL",
  "QIC_TP",
  "err_1",
  "err_2"
)

# taxa名称放第一列
result_TP$Taxa <- taxa_names

result_TP <- result_TP[
  ,
  c(
    "Taxa",
    setdiff(
      names(result_TP),
      "Taxa"
    )
  )
]


# ============================================================
# QIC-based P selection
# ============================================================

result_TP$pvalue.QIC <- mapply(
  
  function(
    zp,
    p_cont,
    p_global,
    qic_logn,
    qic_tp
  ) {
    
    if (
      is.na(zp)
    ) {
      return(NA_real_)
    }
    
    # 很少0：主要使用continuous part
    if (
      zp <= 0.05
    ) {
      return(
        p_cont
      )
    }
    
    # 较多0：使用two-part global P
    if (
      zp >= 0.15
    ) {
      return(
        p_global
      )
    }
    
    # 中间区域比较QIC
    if (
      zp > 0.05 &&
      zp < 0.15
    ) {
      
      # 两个QIC都存在
      if (
        !is.na(qic_logn) &&
        !is.na(qic_tp)
      ) {
        
        if (
          qic_logn <= qic_tp
        ) {
          
          return(
            p_cont
          )
          
        } else {
          
          return(
            p_global
          )
        }
        
      }
      
      # 只有lognormal QIC
      if (
        !is.na(qic_logn) &&
        is.na(qic_tp)
      ) {
        
        return(
          p_cont
        )
      }
      
      # 只有TP QIC
      if (
        is.na(qic_logn) &&
        !is.na(qic_tp)
      ) {
        
        return(
          p_global
        )
      }
      
      return(
        NA_real_
      )
    }
    
    return(
      NA_real_
    )
  },
  
  zp =
    result_TP$ZP,
  
  p_cont =
    result_TP$pvalue.cont,
  
  p_global =
    result_TP$global.p,
  
  qic_logn =
    result_TP$QIC_LOGNORMAL,
  
  qic_tp =
    result_TP$QIC_TP
)


# ============================================================
# BH-FDR
# ============================================================

result_TP$qvalue <- p.adjust(
  result_TP$pvalue.QIC,
  method = "BH"
)


# 按q-value排序
result_TP <- result_TP[
  order(
    result_TP$qvalue
  ),
]


# ============================================================
# Save
# ============================================================

write.csv(
  result_TP,
  "NOVA_kathy_kcal_4.0_genus_twopart_adjusted-updated.csv",
  row.names = FALSE
)

# 看前20个
head(
  result_TP,
  20
)
