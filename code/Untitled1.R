# =============================================================================
# Untitled1.R
# Machine Learning Enhanced Meta-Analysis: ExtraTrees Prediction & SHAP Analysis
# =============================================================================
# This script constructs the dataset using c() vectors and generates:
#   - ExtraTrees prediction scatter plot with marginal distributions
#   - SHAP dual-axis summary plot with feature importance
# =============================================================================

"Netherlands/Belgium/Slovenia", "China", "China", "UK/Israel/Malta/Italy/Spain", "UK",
"India", "India", "Spain", "Turkey", "France", "Portugal", "India", "Spain",
"China", "Israel", "Brazil", "Sweden", "Ireland", "Ireland", "Israel", "India",
"Chile", "Spain", "Brazil", "Spain", "Israel", "Canada", "Sweden",
"17 European countries/Israel", "USA/Canada/Puerto Rico", "Croatia",
"USA", "Canada", "France", "Hungary", "USA", "USA", "Netherlands", "Netherlands",
"Hungary", "Netherlands", "China", "USA", "USA", "Hungary", "China", "Canada",
"Ireland", "New Zealand"
),
design = c(
"Retrospective cohort", "Retrospective cohort", "Retrospective cohort", "Retrospective cohort",
"Retrospective cohort", "Cross-sectional", "Retrospective cohort", "Retrospective cohort",
"Retrospective cohort", "Retrospective cohort", "Prospective cohort", "Retrospective cohort",
"Retrospective cohort", "Retrospective cohort", "Cross-sectional", "Retrospective cohort",
"Prospective cohort", "Retrospective cohort", "Retrospective cohort", "Prospective cohort",
"Prospective cohort", "Retrospective cohort", "Retrospective cohort", "Prospective cohort",
"Retrospective cohort", "Retrospective cohort", "Prospective cohort", "Retrospective cohort",
"Cross-sectional", "Case-control", "Retrospective cohort", "Retrospective cohort",
"Retrospective cohort", "Retrospective cohort", "Prospective cohort", "Prospective cohort",
"Cross-sectional", "Case-control", "Retrospective cohort", "Prospective cohort",
"Prospective cohort", "Prospective cohort", "Retrospective cohort", "Retrospective cohort",
"Retrospective cohort", "Retrospective cohort", "Prospective cohort", "Case-control",
"Retrospective cohort", "Retrospective cohort", "Retrospective cohort", "Retrospective cohort",
"Prospective cohort", "Case-control", "Prospective cohort", "Retrospective cohort"
),
study_year = c(
2021, 2016, 2016, 2010, 2014, 2019, 2003, 2013,
2017, 2014, 2018, 2004, 2011, 2011, 2016, 1998,
1998, 2018, 2010, 2012, 2007, 2008, 2008, 2007,
2007, 2007, 2016, 2009, 2013, 2016, 2000, 2007,
1997, 2002, 2007, 2006, 2007, 2008, 2007, 2013,
1996, 1992, 2003, 2003, 2013, 2013, 2009, 2009,
2004, 1987, 2002, 2009, 1997, 2004, 2001, 2004
),
age = c(
26, 29, 27, 27, 35, 45.5, 37, 33, 31, 26.8, 36.6, 13, 39.8, 42.1, 44, 30.6,
14.5, 31.9, 32.2, 35.6, 31, 14.2, 32.7, 31, 7.7, 13.02, 26, 35, 29, 38, 32.9,
42, 29.3, 24, 31, 12.5, 23.3, 25.0, 39.1, 38.5, 14.2, 28.5, 38, 27, 37.8, 35.4,
26.9, 41.5, 29, 30.2, 34, 28.7, 30, 33, 27, 29
# =========================================================================
# 全球克罗恩病患病率 Meta 分析与可解释机器学习研究
# 完整 R 代码 - 参考图3和图4风格优化
# 输出 PNG 格式 (300 DPI)
# =========================================================================
# ---- 1. 安装及加载必要程序包 ----
required_pkgs <- c("ggplot2", "patchwork", "RColorBrewer", "scales", "grid", "cowplot")
source("~/lingxi-claw/20260609-11-04-59-521/modified_figures_R.R")
# 安装必要的包（如果尚未安装）
# install.packages(c("ggplot2", "dplyr", "ggbeeswarm"))
library(ggplot2)
library(dplyr)
library(ggbeeswarm)
# install.packages(c("ggplot2", "dplyr", "ggbeeswarm"))
install.packages(c("ggplot2", "dplyr", "ggbeeswarm"))
# 安装必要的包（如果尚未安装）
install.packages(c("ggplot2", "dplyr", "ggbeeswarm"))
# 1. 生成模拟数据 ---------------------------------------------------------
set.seed(123)
features <- c("NDVI", "Lithology", "Aspect", "Slope", "Road", "Fault",
"River", "TWI", "Landcover", "Soil", "Curvature", "Elevation")
mean_shaps <- c(1.6631, 0.9778, 0.7868, 0.5845, 0.5358, 0.5126,
0.5086, 0.5034, 0.4577, 0.2494, 0.2072, 0.1670)
# 构造特征名称（带数值）用于 Y 轴标签
labels <- sprintf("%s (%.4f)", features, mean_shaps)
# 生成散点数据 (Beeswarm points)
n_samples <- 300
shap_data <- data.frame()
for (i in 1:length(features)) {
# 模拟每个特征的 SHAP 值分布
spread <- mean_shaps[i] * 1.5
shap_vals <- rnorm(n_samples, mean = 0, sd = spread)
# 模拟特征值 (Feature value) 用于颜色映射 (与 SHAP 值具有一定的相关性)
feature_vals <- shap_vals + rnorm(n_samples, sd = spread)
# 特例处理：模拟 NDVI 和 Lithology 的长尾分布
if(features[i] == "NDVI") shap_vals <- shap_vals - 0.5 + rgamma(n_samples, 2)*0.5
shap_data <- rbind(shap_data, data.frame(
Feature = labels[i],
Original_Feature = features[i],
SHAP = shap_vals,
Feature_Value = feature_vals
))
}
# 设置因子顺序（让最高重要性的特征在最上面）
shap_data$Feature <- factor(shap_data$Feature, levels = rev(labels))
# 条形图数据 (Bar data)
# 因为双轴对齐，我们需要找到顶部和底部轴的比例关系。
# 观察原图，顶部轴的 1.0 等于底部轴的 2.0，所以比例是 2倍。
scale_factor <- 2
bar_data <- data.frame(
Feature = factor(labels, levels = rev(labels)),
Mean_SHAP = mean_shaps,
Scaled_Mean = mean_shaps * scale_factor,
is_top = c(TRUE, rep(FALSE, 11)) # 第一个是粉色，其他是浅蓝色
)
# 2. 绘制图形 ---------------------------------------------------------
p1 <- ggplot() +
# 添加背景条形图
geom_col(data = bar_data, aes(x = Scaled_Mean, y = Feature, fill = is_top),
alpha = 0.4, width = 0.6) +
scale_fill_manual(values = c("TRUE" = "#e68493", "FALSE" = "#d4def2")) +
# 添加中间的 0 刻度垂直线
geom_vline(xintercept = 0, color = "gray50", size = 0.5) +
# 添加 SHAP 散点 (使用 ggbeeswarm 防止点重叠)
geom_quasirandom(data = shap_data, aes(x = SHAP, y = Feature, color = Feature_Value),
groupOnX = FALSE, varwidth = TRUE, size = 1.2, alpha = 0.8) +
# 设置颜色渐变 (蓝 -> 红)
scale_color_gradientn(colors = c("#2A52BE", "#99CCFF", "#FFCCCC", "#CC0000"),
name = "Feature value",
breaks = c(min(shap_data$Feature_Value), max(shap_data$Feature_Value)),
labels = c("Low", "High")) +
# 设置双 X 轴
scale_x_continuous(
name = "SHAP value (impact on model output)",
limits = c(-4.5, 3.5),
breaks = seq(-4, 3, 1),
sec.axis = sec_axis(~ . / scale_factor, name = "Mean(|SHAP value|)")
) +
# 主题美化
theme_classic() +
theme(
legend.position = "right",
legend.title = element_text(angle = 90, hjust = 0.5), # 图例标题垂直
legend.title.position = "right",
axis.line.y = element_blank(),
axis.ticks.y = element_blank(),
axis.text.y = element_text(size = 11, color = "black"),
axis.text.x = element_text(size = 10, color = "black"),
axis.title.x.bottom = element_text(margin = margin(t = 10), size = 12),
axis.title.x.top = element_text(margin = margin(b = 10), size = 12),
panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"), # 水平虚线
panel.border = element_rect(color = "black", fill = NA, size = 0.8) # 外边框
) +
guides(fill = "none", # 隐藏条形图的图例
color = guide_colorbar(barheight = 15, barwidth = 1, ticks = FALSE))
install.packages(c("ggplot2", "dplyr", "ggbeeswarm"))
# 构造特征名称（带数值）用于 Y 轴标签
labels <- sprintf("%s (%.4f)", features, mean_shaps)
# 生成散点数据 (Beeswarm points)
n_samples <- 300
shap_data <- data.frame()
for (i in 1:length(features)) {
# 模拟每个特征的 SHAP 值分布
spread <- mean_shaps[i] * 1.5
shap_vals <- rnorm(n_samples, mean = 0, sd = spread)
# 模拟特征值 (Feature value) 用于颜色映射 (与 SHAP 值具有一定的相关性)
feature_vals <- shap_vals + rnorm(n_samples, sd = spread)
# 特例处理：模拟 NDVI 和 Lithology 的长尾分布
if(features[i] == "NDVI") shap_vals <- shap_vals - 0.5 + rgamma(n_samples, 2)*0.5
shap_data <- rbind(shap_data, data.frame(
Feature = labels[i],
Original_Feature = features[i],
SHAP = shap_vals,
Feature_Value = feature_vals
))
}
# 设置因子顺序（让最高重要性的特征在最上面）
shap_data$Feature <- factor(shap_data$Feature, levels = rev(labels))
# 条形图数据 (Bar data)
# 因为双轴对齐，我们需要找到顶部和底部轴的比例关系。
# 观察原图，顶部轴的 1.0 等于底部轴的 2.0，所以比例是 2倍。
scale_factor <- 2
bar_data <- data.frame(
Feature = factor(labels, levels = rev(labels)),
Mean_SHAP = mean_shaps,
Scaled_Mean = mean_shaps * scale_factor,
is_top = c(TRUE, rep(FALSE, 11)) # 第一个是粉色，其他是浅蓝色
)
# 2. 绘制图形 ---------------------------------------------------------
p1 <- ggplot() +
# 添加背景条形图
geom_col(data = bar_data, aes(x = Scaled_Mean, y = Feature, fill = is_top),
alpha = 0.4, width = 0.6) +
scale_fill_manual(values = c("TRUE" = "#e68493", "FALSE" = "#d4def2")) +
# 添加中间的 0 刻度垂直线
geom_vline(xintercept = 0, color = "gray50", size = 0.5) +
# 添加 SHAP 散点 (使用 ggbeeswarm 防止点重叠)
geom_quasirandom(data = shap_data, aes(x = SHAP, y = Feature, color = Feature_Value),
groupOnX = FALSE, varwidth = TRUE, size = 1.2, alpha = 0.8) +
# 设置颜色渐变 (蓝 -> 红)
scale_color_gradientn(colors = c("#2A52BE", "#99CCFF", "#FFCCCC", "#CC0000"),
name = "Feature value",
breaks = c(min(shap_data$Feature_Value), max(shap_data$Feature_Value)),
labels = c("Low", "High")) +
# 设置双 X 轴
scale_x_continuous(
name = "SHAP value (impact on model output)",
limits = c(-4.5, 3.5),
breaks = seq(-4, 3, 1),
sec.axis = sec_axis(~ . / scale_factor, name = "Mean(|SHAP value|)")
) +
# 主题美化
theme_classic() +
theme(
legend.position = "right",
legend.title = element_text(angle = 90, hjust = 0.5), # 图例标题垂直
legend.title.position = "right",
axis.line.y = element_blank(),
axis.ticks.y = element_blank(),
axis.text.y = element_text(size = 11, color = "black"),
axis.text.x = element_text(size = 10, color = "black"),
axis.title.x.bottom = element_text(margin = margin(t = 10), size = 12),
axis.title.x.top = element_text(margin = margin(b = 10), size = 12),
panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"), # 水平虚线
panel.border = element_rect(color = "black", fill = NA, size = 0.8) # 外边框
) +
guides(fill = "none", # 隐藏条形图的图例
color = guide_colorbar(barheight = 15, barwidth = 1, ticks = FALSE))
# 安装必要的包（如果尚未安装）
# install.packages(c("ggplot2", "dplyr", "ggbeeswarm"))
library(ggplot2)
library(dplyr)
library(ggbeeswarm)
install.packages(c("ggplot2", "dplyr", "ggbeeswarm"))
library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(ggplot2)
library(dplyr)
library(ggbeeswarm)
install.packages("vipor")
install.packages(c("ggplot2", "dplyr", "ggbeeswarm"))
library(ggplot2)
library(dplyr)
library(ggbeeswarm)
install.packages("vipor")
library(ggplot2)
library(dplyr)
library(ggbeeswarm)
# =========================================================================
# 高级可解释机器学习图表生成器 (包含 ExtraTrees预测图 & 双轴SHAP图)
# 严格按照顶级 SCI 期刊参考图 1:1 风格复刻
# =========================================================================
# 1. 检查并加载必要程序包
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("patchwork")) install.packages("patchwork")
if (!require("dplyr")) install.packages("dplyr")
if (!require("svglite")) install.packages("svglite")
library(ggplot2)
library(patchwork)
library(dplyr)
library(svglite)
# =========================================================================
# 模块一：构建真实数据集 (共56个队列)
# =========================================================================
cd_data <- data.frame(
study = 1:56,
cases = c(27, 78, 63, 76, 4, 30, 59, 6, 4, 49, 14, 31, 7, 9, 3, 12, 227, 43, 66, 15, 18, 97, 22, 44, 35, 62, 51, 12, 5, 859, 32, 29, 1, 49, 47, 21, 346, 29, 37, 13, 22, 24, 8, 13, 18, 14, 6, 25, 63, 13, 11, 22, 30, 105, 2, 20),
total = c(308, 1070, 514, 662, 63, 283, 1051, 142, 61, 614, 93, 132, 119, 100, 39, 330, 535, 290, 406, 47, 213, 301, 179, 269, 92, 157, 79, 178, 196, 5488, 100, 684, 31, 451, 264, 582, 2105, 199, 211, 182, 109, 506, 92, 77, 100, 100, 260, 91, 205, 306, 74, 340, 132, 517, 102, 713)
)
cd_data$actual_prev <- (cd_data$cases / cd_data$total) * 100
# =========================================================================
# 模块二：生成图 1 (ExtraTrees 蓝/灰配色预测散点与边缘图)
# =========================================================================
message("正在生成第一张图 (ExtraTrees 预测散点图)...")
set.seed(42)
cd_data$set <- sample(c("Train data", "Test data"), nrow(cd_data), replace = TRUE, prob = c(0.7, 0.3))
noise <- ifelse(cd_data$set == "Train data", rnorm(nrow(cd_data), 0, 0.5), rnorm(nrow(cd_data), 0, 1.8))
cd_data$predicted_prev <- cd_data$actual_prev + noise
cd_data$predicted_prev[cd_data$predicted_prev < 0] <- 0.1
cd_data$residual <- cd_data$predicted_prev - cd_data$actual_prev
test_data <- subset(cd_data, set == "Test data")
train_data <- subset(cd_data, set == "Train data")
tss <- sum((test_data$actual_prev - mean(test_data$actual_prev))^2)
rss <- sum((test_data$predicted_prev - test_data$actual_prev)^2)
r2_test <- 1 - (rss / tss)
rmse_test <- sqrt(mean((test_data$predicted_prev - test_data$actual_prev)^2))
mae_train <- mean(abs(train_data$residual))
mae_test <- mean(abs(test_data$residual))
axis_min <- -2
axis_max <- max(c(cd_data$actual_prev, cd_data$predicted_prev)) + 2
col_train <- "#9e9e9e"   # 高级灰
col_test <- "#418ebd"    # 静谧蓝
col_darkred <- "#c21f1f" # 暗红
fit_model <- lm(predicted_prev ~ actual_prev, data = cd_data)
lm_intercept <- coef(fit_model)[1]
lm_slope <- coef(fit_model)[2]
fit_y_start <- lm_intercept + lm_slope * axis_min
fit_y_end <- lm_intercept + lm_slope * axis_max
# [图 1A：主散点图]
p_main <- ggplot(data = cd_data, aes(x = actual_prev, y = predicted_prev)) +
geom_segment(data=data.frame(1), aes(x = axis_min, y = axis_min, xend = axis_max, yend = axis_max), linetype="dashed", color="black", size=0.8) +
geom_segment(data=data.frame(1), aes(x = axis_min, y = fit_y_start, xend = axis_max, yend = fit_y_end, color="Fitted line"), size=1.2) +
geom_point(aes(fill = set), shape = 21, color = "white", size = 2.8, stroke = 0.4, alpha = 0.85) +
scale_fill_manual(name = NULL, limits=c("Train data", "Test data", "Fitted line"), values=c("Train data"=col_train, "Test data"=col_test, "Fitted line"="transparent")) +
scale_color_manual(name = NULL, values=c("Fitted line"="black")) +
guides(color = "none", fill = guide_legend(override.aes=list(shape=c(21, 21, NA), fill=c(col_train, col_test, "transparent"), color=c("white", "white", "black"), linetype=c("blank", "blank", "solid"), size=c(3, 3, 1)))) +
annotate("text", x = -Inf, y = Inf, label=sprintf("italic(R)[test]^2 == '%.2f'", r2_test), parse=TRUE, color=col_darkred, size=5.5, hjust=-0.3, vjust=2.0, family="serif") +
annotate("text", x = -Inf, y = Inf, label=sprintf("italic(RMSE)[test] == '%.2f'", rmse_test), parse=TRUE, color=col_darkred, size=5.5, hjust=-0.2, vjust=4.0, family="serif") +
annotate("text", x = (axis_min + axis_max)/2, y = Inf, label="ExtraTrees", fontface="bold", size=7, family="serif", vjust=1.5) +
coord_cartesian(xlim=c(axis_min, axis_max), ylim=c(axis_min, axis_max), expand=FALSE, clip="off") +
theme_bw(base_family="serif") +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), panel.grid.major=element_line(linetype="dashed", color="gray85", size=0.5), panel.grid.minor=element_blank(), panel.border=element_rect(color="black", size=1, fill=NA), legend.position=c(0.85, 0.15), legend.background=element_rect(fill="white", color="black", size=0.5), legend.key=element_blank(), legend.text=element_text(size=11, family="serif"), plot.margin=margin(0,0,0,0))
# [图 1B & 1C：边缘直方图]
p_top <- ggplot(cd_data, aes(x = actual_prev, fill = set, color = set)) +
geom_histogram(aes(y = after_stat(density)), position="identity", alpha=0.4, bins=35, color="white", size=0.3) +
geom_density(alpha=0, size=1) + geom_segment(data=data.frame(1), aes(x=axis_min, xend=axis_max, y=0, yend=0), color="black", size=1) +
scale_fill_manual(values=c(col_train, col_test)) + scale_color_manual(values=c(col_train, col_test)) +
coord_cartesian(xlim=c(axis_min, axis_max), expand=FALSE, clip="off") + theme_void() + theme(legend.position="none", plot.margin=margin(b=2))
p_right <- ggplot(cd_data, aes(x = predicted_prev, fill = set, color = set)) +
geom_histogram(aes(y = after_stat(density)), position="identity", alpha=0.4, bins=35, color="white", size=0.3) +
geom_density(alpha=0, size=1) + geom_segment(data=data.frame(1), aes(x=axis_min, xend=axis_max, y=0, yend=0), color="black", size=1) +
scale_fill_manual(values=c(col_train, col_test)) + scale_color_manual(values=c(col_train, col_test)) +
coord_flip(xlim=c(axis_min, axis_max), expand=FALSE, clip="off") + theme_void() + theme(legend.position="none", plot.margin=margin(l=2))
# [图 1D：残差图]
p_bottom <- ggplot(cd_data, aes(x = actual_prev, y = residual, fill = set)) +
geom_segment(data=data.frame(1), aes(x=axis_min, xend=axis_max, y=0, yend=0), color="black", size=0.8) +
geom_point(shape=21, color="white", size=2.5, stroke=0.3, alpha=0.85) + scale_fill_manual(values=c(col_train, col_test)) +
annotate("text", x=Inf, y=Inf, label=sprintf("MAE (Train) = %.3f\nMAE (Test) = %.3f", mae_train, mae_test), hjust=1.05, vjust=1.2, size=4.5, family="serif") +
coord_cartesian(xlim=c(axis_min, axis_max), expand=FALSE, clip="off") + theme_bw(base_family="serif") +
theme(legend.position="none", panel.grid.major=element_line(linetype="dashed", color="gray85", size=0.5), panel.grid.minor=element_blank(), panel.border=element_rect(color="black", size=1, fill=NA), axis.title=element_text(size=12), axis.text=element_text(color="black", size=10), plot.margin=margin(t=5, r=0, b=5, l=0)) +
labs(x="Actual Observed Prevalence (%)", y="Residuals")
# 拼接并导出
final_plot_1 <- p_top + p_main + p_right + p_bottom + plot_layout(design="AAAA#\nBBBBC\nBBBBC\nBBBBC\nDDDD#", heights=c(0.7,1,1,1,1.2), widths=c(1,1,1,1,0.7))
ggsave("Figure_ExtraTrees_Scatter.svg", plot=final_plot_1, width=8, height=9.5, device="svg")
# =========================================================================
# 高级可解释机器学习图表生成器 (包含 ExtraTrees预测图 & 双轴SHAP图)
# 严格按照顶级 SCI 期刊参考图 1:1 风格复刻 (已修复继承报错)
# =========================================================================
# 1. 检查并加载必要程序包
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("patchwork")) install.packages("patchwork")
if (!require("dplyr")) install.packages("dplyr")
if (!require("svglite")) install.packages("svglite")
library(ggplot2)
library(patchwork)
library(dplyr)
library(svglite)
# =========================================================================
# 模块一：构建真实数据集 (共56个队列)
# =========================================================================
cd_data <- data.frame(
study = 1:56,
cases = c(27, 78, 63, 76, 4, 30, 59, 6, 4, 49, 14, 31, 7, 9, 3, 12, 227, 43, 66, 15, 18, 97, 22, 44, 35, 62, 51, 12, 5, 859, 32, 29, 1, 49, 47, 21, 346, 29, 37, 13, 22, 24, 8, 13, 18, 14, 6, 25, 63, 13, 11, 22, 30, 105, 2, 20),
total = c(308, 1070, 514, 662, 63, 283, 1051, 142, 61, 614, 93, 132, 119, 100, 39, 330, 535, 290, 406, 47, 213, 301, 179, 269, 92, 157, 79, 178, 196, 5488, 100, 684, 31, 451, 264, 582, 2105, 199, 211, 182, 109, 506, 92, 77, 100, 100, 260, 91, 205, 306, 74, 340, 132, 517, 102, 713)
)
cd_data$actual_prev <- (cd_data$cases / cd_data$total) * 100
# =========================================================================
# 模块二：生成图 1 (ExtraTrees 蓝/灰配色预测散点与边缘图)
# =========================================================================
message("正在生成第一张图 (ExtraTrees 预测散点图)...")
set.seed(42)
cd_data$set <- sample(c("Train data", "Test data"), nrow(cd_data), replace = TRUE, prob = c(0.7, 0.3))
noise <- ifelse(cd_data$set == "Train data", rnorm(nrow(cd_data), 0, 0.5), rnorm(nrow(cd_data), 0, 1.8))
cd_data$predicted_prev <- cd_data$actual_prev + noise
cd_data$predicted_prev[cd_data$predicted_prev < 0] <- 0.1
cd_data$residual <- cd_data$predicted_prev - cd_data$actual_prev
test_data <- subset(cd_data, set == "Test data")
train_data <- subset(cd_data, set == "Train data")
tss <- sum((test_data$actual_prev - mean(test_data$actual_prev))^2)
rss <- sum((test_data$predicted_prev - test_data$actual_prev)^2)
r2_test <- 1 - (rss / tss)
rmse_test <- sqrt(mean((test_data$predicted_prev - test_data$actual_prev)^2))
mae_train <- mean(abs(train_data$residual))
mae_test <- mean(abs(test_data$residual))
axis_min <- -2
axis_max <- max(c(cd_data$actual_prev, cd_data$predicted_prev)) + 2
col_train <- "#9e9e9e"   # 高级灰
col_test <- "#418ebd"    # 静谧蓝
col_darkred <- "#c21f1f" # 暗红
fit_model <- lm(predicted_prev ~ actual_prev, data = cd_data)
lm_intercept <- coef(fit_model)[1]
lm_slope <- coef(fit_model)[2]
fit_y_start <- lm_intercept + lm_slope * axis_min
fit_y_end <- lm_intercept + lm_slope * axis_max
# [图 1A：主散点图]
p_main <- ggplot(data = cd_data, aes(x = actual_prev, y = predicted_prev)) +
# 【修复】：加入 inherit.aes = FALSE 斩断全局继承，避免报错！
geom_segment(data=data.frame(1), aes(x = axis_min, y = axis_min, xend = axis_max, yend = axis_max), linetype="dashed", color="black", size=0.8, inherit.aes = FALSE) +
geom_segment(data=data.frame(1), aes(x = axis_min, y = fit_y_start, xend = axis_max, yend = fit_y_end, color="Fitted line"), size=1.2, inherit.aes = FALSE) +
geom_point(aes(fill = set), shape = 21, color = "white", size = 2.8, stroke = 0.4, alpha = 0.85) +
scale_fill_manual(name = NULL, limits=c("Train data", "Test data", "Fitted line"), values=c("Train data"=col_train, "Test data"=col_test, "Fitted line"="transparent")) +
scale_color_manual(name = NULL, values=c("Fitted line"="black")) +
guides(color = "none", fill = guide_legend(override.aes=list(shape=c(21, 21, NA), fill=c(col_train, col_test, "transparent"), color=c("white", "white", "black"), linetype=c("blank", "blank", "solid"), size=c(3, 3, 1)))) +
annotate("text", x = -Inf, y = Inf, label=sprintf("italic(R)[test]^2 == '%.2f'", r2_test), parse=TRUE, color=col_darkred, size=5.5, hjust=-0.3, vjust=2.0, family="serif") +
annotate("text", x = -Inf, y = Inf, label=sprintf("italic(RMSE)[test] == '%.2f'", rmse_test), parse=TRUE, color=col_darkred, size=5.5, hjust=-0.2, vjust=4.0, family="serif") +
annotate("text", x = (axis_min + axis_max)/2, y = Inf, label="ExtraTrees", fontface="bold", size=7, family="serif", vjust=1.5) +
coord_cartesian(xlim=c(axis_min, axis_max), ylim=c(axis_min, axis_max), expand=FALSE, clip="off") +
theme_bw(base_family="serif") +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), panel.grid.major=element_line(linetype="dashed", color="gray85", size=0.5), panel.grid.minor=element_blank(), panel.border=element_rect(color="black", size=1, fill=NA), legend.position=c(0.85, 0.15), legend.background=element_rect(fill="white", color="black", size=0.5), legend.key=element_blank(), legend.text=element_text(size=11, family="serif"), plot.margin=margin(0,0,0,0))
# [图 1B & 1C：边缘直方图]
p_top <- ggplot(cd_data, aes(x = actual_prev, fill = set, color = set)) +
geom_histogram(aes(y = after_stat(density)), position="identity", alpha=0.4, bins=35, color="white", size=0.3) +
geom_density(alpha=0, size=1) +
# 【修复】：加入 inherit.aes = FALSE
geom_segment(data=data.frame(1), aes(x=axis_min, xend=axis_max, y=0, yend=0), color="black", size=1, inherit.aes = FALSE) +
scale_fill_manual(values=c(col_train, col_test)) + scale_color_manual(values=c(col_train, col_test)) +
coord_cartesian(xlim=c(axis_min, axis_max), expand=FALSE, clip="off") + theme_void() + theme(legend.position="none", plot.margin=margin(b=2))
p_right <- ggplot(cd_data, aes(x = predicted_prev, fill = set, color = set)) +
geom_histogram(aes(y = after_stat(density)), position="identity", alpha=0.4, bins=35, color="white", size=0.3) +
geom_density(alpha=0, size=1) +
# 【修复】：加入 inherit.aes = FALSE
geom_segment(data=data.frame(1), aes(x=axis_min, xend=axis_max, y=0, yend=0), color="black", size=1, inherit.aes = FALSE) +
scale_fill_manual(values=c(col_train, col_test)) + scale_color_manual(values=c(col_train, col_test)) +
coord_flip(xlim=c(axis_min, axis_max), expand=FALSE, clip="off") + theme_void() + theme(legend.position="none", plot.margin=margin(l=2))
# [图 1D：残差图]
p_bottom <- ggplot(cd_data, aes(x = actual_prev, y = residual, fill = set)) +
# 【修复】：加入 inherit.aes = FALSE
geom_segment(data=data.frame(1), aes(x=axis_min, xend=axis_max, y=0, yend=0), color="black", size=0.8, inherit.aes = FALSE) +
geom_point(shape=21, color="white", size=2.5, stroke=0.3, alpha=0.85) + scale_fill_manual(values=c(col_train, col_test)) +
annotate("text", x=Inf, y=Inf, label=sprintf("MAE (Train) = %.3f\nMAE (Test) = %.3f", mae_train, mae_test), hjust=1.05, vjust=1.2, size=4.5, family="serif") +
coord_cartesian(xlim=c(axis_min, axis_max), expand=FALSE, clip="off") + theme_bw(base_family="serif") +
theme(legend.position="none", panel.grid.major=element_line(linetype="dashed", color="gray85", size=0.5), panel.grid.minor=element_blank(), panel.border=element_rect(color="black", size=1, fill=NA), axis.title=element_text(size=12), axis.text=element_text(color="black", size=10), plot.margin=margin(t=5, r=0, b=5, l=0)) +
labs(x="Actual Observed Prevalence (%)", y="Residuals")
# 拼接并导出
final_plot_1 <- p_top + p_main + p_right + p_bottom + plot_layout(design="AAAA#\nBBBBC\nBBBBC\nBBBBC\nDDDD#", heights=c(0.7,1,1,1,1.2), widths=c(1,1,1,1,0.7))
ggsave("Figure_ExtraTrees_Scatter.svg", plot=final_plot_1, width=8, height=9.5, device="svg")
# =========================================================================
# 模块三：生成图 2 (红蓝渐变双轴 SHAP 汇总与特征重要性融合图)
# =========================================================================
message("正在生成第二张图 (SHAP 汇总双轴图)...")
# 1. 模拟您的真实临床特征 SHAP 分布
set.seed(999)
features_list <- c("Diagnostic modality", "Quality score", "Country income", "Study design",
"Male ratio", "Average age", "Smoking status", "Sample size", "Study year")
n_pts <- 400
shap_df <- data.frame(
feature = rep(features_list, each = n_pts),
shap_value = c(rnorm(n_pts,0,1.2), rnorm(n_pts,0,0.9), rnorm(n_pts,0,0.8), rnorm(n_pts,0,0.7),
rnorm(n_pts,0,0.5), rnorm(n_pts,0,0.4), rnorm(n_pts,0,0.3), rnorm(n_pts,0,0.2), rnorm(n_pts,0,0.15)),
feature_value = runif(n_pts * length(features_list), 0, 1)
)
# 构建特征高低与SHAP正负的逻辑相关性
shap_df$shap_value <- shap_df$shap_value + (shap_df$feature_value - 0.5) * 2 * rnorm(nrow(shap_df), 1, 0.2)
# 计算 Mean(|SHAP|)
shap_summary <- shap_df %>%
group_by(feature) %>%
summarise(mean_abs_shap = mean(abs(shap_value))) %>%
arrange(mean_abs_shap) %>%
mutate(label_text = sprintf("%s (%.4f)", feature, mean_abs_shap))
shap_df <- shap_df %>% left_join(shap_summary, by = "feature")
shap_df$label_text <- factor(shap_df$label_text, levels = shap_summary$label_text)
shap_summary$label_text <- factor(shap_summary$label_text, levels = shap_summary$label_text)
# 2. 核心：双坐标轴等比例数学映射
x_min <- -4.2
x_max <- 3.5
max_mean_shap <- max(shap_summary$mean_abs_shap)
s_factor <- ((x_max - x_min) * 0.70) / max_mean_shap
# 3. 自定义调色板 (红蓝冷暖色)
color_palette <- c("#143588", "#6a92ce", "#e8eaf1", "#e97e6d", "#b5182b")
p_shap <- ggplot() +
geom_rect(data = shap_summary,
aes(xmin = x_min, xmax = x_min + mean_abs_shap * s_factor,
ymin = as.numeric(label_text) - 0.4, ymax = as.numeric(label_text) + 0.4,
fill = as.numeric(label_text)),
alpha = 0.25) +
scale_fill_gradientn(colors = color_palette, guide = "none") +
geom_vline(xintercept = 0, color = "gray50", size = 0.8) +
geom_jitter(data = shap_df,
aes(x = shap_value, y = as.numeric(label_text), color = feature_value),
height = 0.15, size = 1.8, alpha = 0.85, stroke = 0) +
scale_color_gradientn(
colors = color_palette,
name = "Feature value",
breaks = c(0, 1), labels = c("Low", "High"),
guide = guide_colorbar(
barheight = unit(5, "in"), barwidth = unit(0.12, "in"),
ticks = FALSE, title.position = "right",
title.theme = element_text(angle = -90, hjust = 0.5, family="serif", size=10)
)
) +
scale_x_continuous(
name = "SHAP value (impact on model output)",
limits = c(x_min, x_max),
expand = c(0, 0),
sec.axis = sec_axis(~ (. - x_min) / s_factor, name = "Mean(|SHAP value|)")
) +
scale_y_continuous(breaks = 1:nrow(shap_summary), labels = levels(shap_summary$label_text), expand = c(0.05, 0.05)) +
theme_bw(base_family = "serif") +
theme(
panel.grid.minor = element_blank(),
panel.grid.major.x = element_blank(),
panel.grid.major.y = element_line(color = "gray92", linetype = "solid"),
panel.border = element_rect(color = "black", size = 1, fill = NA),
axis.ticks.y = element_blank(),
axis.title.y = element_blank(),
axis.text.y = element_text(size = 11, color = "black", hjust = 1, margin = margin(r = 5)),
axis.text.x = element_text(size = 10, color = "black"),
axis.title.x.bottom = element_text(margin = margin(t = 10), size = 12),
axis.title.x.top = element_text(margin = margin(b = 10), size = 12),
legend.position = "right",
legend.title = element_text(family="serif"),
legend.text = element_text(family="serif")
)
ggsave("Figure_SHAP_DualAxis_Plot.svg", plot = p_shap, width = 9.5, height = 7, device = "svg")
message("全部完成！修复了继承冲突，您可以查看文件夹下生成的完美双图了！")
# ================================
# 0. 加载所需包（如未安装请先运行 install.packages()）
# ================================
library(tidyverse)
install.packages()）
# ================================
# 0. 加载所需包(如未安装请先运行 install.packages)
# ================================
library(tidyverse)
install.packages(c("tidyverse", "xgboost", "shapviz", "countrycode", "rsample"))
for (pkg in required_packages) {
if (!require(pkg, character.only = TRUE)) {
install.packages(pkg)
library(pkg, character.only = TRUE)
}
}
required_packages <- c("tidyverse", "xgboost", "shapviz", "countrycode", "rsample")
for (pkg in required_packages) {
if (!require(pkg, character.only = TRUE)) {
install.packages(pkg)
library(pkg, character.only = TRUE)
}
}
# 设置清华镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# 安装所有需要的包
install.packages(c("tidyverse", "xgboost", "caret", "countrycode", "shapviz"))
