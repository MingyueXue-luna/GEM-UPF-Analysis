rm(list = ls())
library(survival)
library(stringr)

setwd("C:/Users/Administrator/Desktop/Paper/NOVA")
aa <- read.csv("Nova_genus_interaction.csv")

# 家系/聚类与结局变量
aa$fam <- str_sub(aa$SubjectID, 1, -4)
aa$CD1 <- ifelse(aa$CD == 1, 1, 0)
colnames(aa)
# 变量列表（按你原来的区间）
variables <- names(aa)[17:120]

# 1) 先对 kathy_kcal_4.0 做 scale 标准化（转成数值向量）
aa$`kathy_kcal_4.0_sc` <- as.numeric(scale(aa$`kathy_kcal_4.0`))

# 用于汇总所有模型的结果
all_results <- data.frame()

for (var in variables) {
  x <- aa[[var]]
  
  # 尽量转数值
  if (!is.numeric(x)) {
    suppressWarnings(x_num <- as.numeric(x))
  } else {
    x_num <- x
  }
  
  # 无法转成数值（全 NA）就直接跳过
  if (all(is.na(x_num))) {
    message("变量 ", var, " 转成数值后全是 NA，跳过。")
    next
  }
  
  # 计算零值比例
  zero_ratio <- tryCatch(mean(x_num == 0, na.rm = TRUE), error = function(e) NA_real_)
  
  # 连续变量 vs 二元化
  if (!is.na(zero_ratio) && zero_ratio > 0.15) {
    new_var_name <- paste0(var, "_bin")
    aa[[new_var_name]] <- ifelse(is.na(x_num), NA, ifelse(x_num == 0, 0, 1))
  } else {
    new_var_name <- paste0(var, "_sc")
    aa[[new_var_name]] <- as.numeric(scale(x_num))
  }
  
  # 👉 再检查一次：这个新变量是不是全 NA 或者没有有效观测
  if (all(is.na(aa[[new_var_name]]))) {
    message("变量 ", new_var_name, " 全是 NA，跳过。")
    next
  }
  
  # 构造本次模型要用到的数据子集
  dat_tmp <- aa[, c("Fellow_up_durations", "CD1", "kathy_kcal_4.0_sc",
                    new_var_name, "Sex", "Age", "fam")]
  rows_use <- complete.cases(dat_tmp)
  
  if (!any(rows_use)) {
    message("变量 ", new_var_name, " 没有完整观测（全有 NA），跳过。")
    next
  }
  
  # 构造公式
  formula_string <- paste0(
    "Surv(Fellow_up_durations, CD1) ~ ",
    "`kathy_kcal_4.0_sc` * `", new_var_name, "`",
    " + Sex + Age + cluster(fam)"
  )
  formula_obj <- as.formula(formula_string)
  
  # 为了方便排查，再包一层 tryCatch
  model <- tryCatch(
    coxph(formula_obj, data = aa[rows_use, ]),
    error = function(e) {
      message("拟合变量 ", var, " 时出错：", e$message)
      return(NULL)
    }
  )
  
  # 如果这个变量拟合失败，就直接下一个
  if (is.null(model)) next
  
  model_summary <- summary(model)
  coef_summary  <- model_summary$coefficients
  coef_names    <- rownames(coef_summary)
  
  exp_coef     <- coef_summary[, "exp(coef)"]
  confint_exp  <- exp(confint(model))
  
  for (i in seq_along(coef_names)) {
    temp_results <- data.frame(
      Model_Variable = var,
      Used_Feature   = new_var_name,
      Coefficient    = coef_names[i],
      Estimate       = coef_summary[i, "coef"],
      Exp_Coef       = exp_coef[i],
      Pr_z           = coef_summary[i, "Pr(>|z|)"],
      Lower_95       = confint_exp[i, 1],
      Upper_95       = confint_exp[i, 2],
      stringsAsFactors = FALSE
    )
    all_results <- rbind(all_results, temp_results)
  }
}

# 写出结果
write.csv(
  all_results,
  file = "interaction_kathy_kcal4.0_scaled_and_taxa_scaled_or_binarized8-25.csv",
  row.names = FALSE
)


grep("Marvin", names(aa), value = TRUE)
marvin_bin <- "d__Bacteria.Firmicutes.Clostridia.Lachnospirales.Lachnospiraceae.Marvinbryantia"
table(
  aa[[marvin_bin]],
  useNA="always"
)

carrier_data <- aa[
  aa$Marvinbryantia_bin == 1,
]

cox_carrier <- coxph(
  Surv(Fellow_up_durations, CD1)~
    kathy_kcal_4.0_sc+
    Sex+
    Age+
    cluster(fam),
  data=carrier_data
)

summary(cox_carrier)

###########
noncarrier_data <- aa[
  aa$Marvinbryantia_bin == 0,
]


cox_noncarrier <- coxph(
  Surv(Fellow_up_durations, CD1)~
    kathy_kcal_4.0_sc+
    Sex+
    Age+
    cluster(fam),
  data=noncarrier_data
)

summary(cox_noncarrier)


interactio <- coxph(
  Surv(Fellow_up_durations,CD1)~
    kathy_kcal_4.0_sc*
    Marvinbryantia_bin+
    Sex+
    Age+
    cluster(fam),
  data=aa
)
summary(interactio)$coefficients
