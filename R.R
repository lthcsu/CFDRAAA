install.packages(c("ggpubr","corrplot",
                   "glmnet","caret","CBCgrps",
                   "tidyverse","rms"))
library(corrplot) 
library(glmnet) 
library(caret) 
library(CBCgrps)
library(nortest)
library(tidyverse)
library(ggpubr)
library(rms)
library(pROC)

#数据整理：主要靠Excel，导入后做轻微设置改动
data==na.omit(data)
head(data)

shapiro.test(data$ANL)
shapiro.test(data$ANAA)
shapiro.test(data$ANAB)
shapiro.test(data$VAAA)
shapiro.test(data$VILT)
shapiro.test(data$'VILT/VAAAA')
shapiro.test(data$Tortuosity)
shapiro.test(data$PWS)
shapiro.test(data$WSS)
shapiro.test(data$MD)
shapiro.test(data$OSI)
shapiro.test(data$RRT)
shapiro.test(data$ECAP)
shapiro.test(data$PWRI)

# 安装并加载所需库
install.packages(c("ggplot2", "dplyr", "MASS", "car", "pROC"))
library(ggplot2)
library(dplyr)
library(MASS)
library(car)
library(pROC)

# ANL和ANAB使用t检验
t_test_ANL <- t.test(data$ANL ~ data$diagnosis)
t_test_ANAB <- t.test(data$ANAB ~ data$diagnosis)

# 其他变量使用Wilcoxon秩和检验
wilcoxon_results <- lapply(data[, c("VAAA", "VILT", "Tortuosity", "PWS", "WSS", "MD", "OSI", "RRT", "ECAP", "PWRI")], function(x) wilcox.test(x ~ data$diagnosis))

# 打印结果
print(t_test_ANL)
print(t_test_ANAB)
print(wilcoxon_results)

# 使用glm进行逻辑回归
logit_model <- glm(diagnosis ~ ANL + ANAB + VAAA + VILT + `VILT/VAAA` + Tortuosity + PWS  + MD  + RRT + ECAP + PWRI, data = data, family = binomial)

# 打印回归结果
summary(logit_model)

# 提取模型系数
coefficients <- coef(summary(logit_model))

# 绘制每个变量的森林图数据框
forest_plot_data <- data.frame(
  Variable = rownames(coefficients),
  Estimate = coefficients[, "Estimate"],
  StdError = coefficients[, "Std. Error"],
  Zvalue = coefficients[, "z value"],
  Pvalue = coefficients[, "Pr(>|z|)"]
)

# 计算95%置信区间
forest_plot_data <- forest_plot_data %>%
  mutate(
    CI_Lower = Estimate - 1.96 * StdError,
    CI_Upper = Estimate + 1.96 * StdError,
    CI = paste(round(CI_Lower, 3), "to", round(CI_Upper, 3)),
    PvalueFormatted = format(Pvalue, scientific = TRUE)
  )

library(dplyr)
# 假设您的forest_plot_data数据框中有三列：Estimate, CI_Lower, CI_Upper
# 添加OR值和95% CI
table_data <- forest_plot_data %>%
  mutate(
    OR = exp(Estimate), # 计算OR值
    `Lower 95% CI` = exp(CI_Lower), # 计算95% CI下界
    `Upper 95% CI` = exp(CI_Upper), # 计算95% CI上界
    `95% CI` = paste0("(", round(`Lower 95% CI`, 2), ", ", round(`Upper 95% CI`, 2), ")") # 创建一个新列来显示95% CI
  ) %>%
  select(Variable, OR, `95% CI`, PvalueFormatted) %>%
  rename(`P Value` = PvalueFormatted)

# 现在，您的table_data包含了您需要的所有列
library(gridExtra)

# 将表格数据转换为grob
table_grob <- tableGrob(table_data)

# 显示表格
grid.draw(table_grob)
library(dplyr)
library(gridExtra)

# 您已经有了table_data数据框
# ...

# 使用tableGrob创建一个表格grob，并调整字体大小和行间距
table_grob <- tableGrob(
  table_data,
  rows = NULL,  # 如果您不想要行名，可以设置为NULL
  theme = ttheme_default(
    core = list(fg_params = list(fontsize = 8)),  # 调整主体内容的字体大小
    colhead = list(fg_params = list(fontsize = 8)),  # 调整列标题的字体大小
    padding = unit(c(2, 4), "mm")  # 调整单元格内的间距，这里是上下2mm，左右4mm
  )
)

# 显示表格
grid.newpage()  # 开始一个新的页面
grid.draw(table_grob)



###lasso
#minmax标准化转换
min_max_scale = function(x){
  (x-min(x))/(max(x)-min(x))
}

data2 = data%>%
  mutate(diagnosis = as.factor(diagnosis))%>%
  mutate_if(.predicate = is.numeric,
            .funs = min_max_scale)%>%
  as.data.frame()

#z转化矩阵
set.seed(123) #random number generator
x <- data.matrix(data2[, -1])
y <- data2[, 1]
y<-as.numeric(unlist(y))

#lasso
lasso <- glmnet(x, y, family = "binomial",nlambda = 1000, alpha = 1)
print(lasso)
plot(lasso, xvar = "lambda", label = TRUE)
#交叉验证
lasso.cv = cv.glmnet(x, y,alpha = 1,nfolds =20,family="binomial")
plot(lasso.cv)
lasso.cv$lambda.min #minimum
lasso.cv$lambda.1se #one standard error away
coef(lasso.cv, s = "lambda.1se")


###logistic
##数据集划分
set.seed(1)
train_id = sample(1:nrow(data),0.7*nrow(data))
train=data[train_id,]
test=data[-train_id,]
# 保存train数据框到csv文件
write.csv(train,file = "train")
write.csv(test,file="test")
###logistics和列线图
mydata<-fold_train
attach(mydata) 
dd<-datadist(mydata) 
options(datadist='dd')
fit0<-lrm(diagnosis ~ ANL+ANAA+ANAB+VAAA+VILT+`VILT/VAAA`+Tortuosity+PWS+WSS+MD+OSI+RRT+ECAP+PWRI, 
          data = mydata, x = T, y = T) 
fit0
nom0 <- nomogram(fit0, fun = plogis,fun.at = c(.001,.01,.05,.5, .95, .99,.999),
                 lp = T, funlabel = "diagnosis rate")  
plot(nom0)
fit1<-lrm(diagnosis ~ ANL+`VILT/VAAA`+Tortuosity+MD, 
          data = mydata, x = T, y = T) 
fit1
summary(fit1)
print(fit1, digits = 10)
predicted_probabilities <- predict(fit1, type = "fitted")

head(predicted_probabilities)
head(predicted_classes)


fit2<-lrm(diagnosis ~ ANL+VAAA+Tortuosity+MD, 
          data = mydata, x = T, y = T) 
fit2
summary(fit2)


##nomogram
nom1 <- nomogram(fit1, fun = plogis,fun.at = c(.005, .01, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95, .99),
                 lp = T, funlabel = "Risk rate")  
plot(nom1)  
print(fit1, digits = 10)
nom2 <- nomogram(fit0, fun = plogis,fun.at = c(.001,.01,.05,.5, .95, .99,.999),
                 lp = T, funlabel = "diagnosis rate")  
plot(nom2) 
print(fit0, digits = 10)
###predict预测并做ROC
gd<-predict(fit1, newdata = train,
            se.fit = FALSE, dispersion = NULL, terms = NULL,
            na.action = na.pass)
gd2<-predict(fit1, newdata = test,
             se.fit = FALSE, dispersion = NULL, terms = NULL,
             na.action = na.pass)


gd1<-predict(fit0, newdata = train,
            se.fit = FALSE, dispersion = NULL, terms = NULL,
            na.action = na.pass)
gd3<-predict(fit0, newdata = test,
             se.fit = FALSE, dispersion = NULL, terms = NULL,
             na.action = na.pass)



#ROC
library(pROC)
library(ggplot2)

# 假设您已经有了训练集的真实结果和预测概率
# train$diagnosis 是真实结果，gd 是预测概率

# 计算ROC曲线
roc.list <- roc(train$diagnosis, gd)
#计算最佳阈值
coords <- coords(roc.list, "best", ret="threshold")
best.threshold <- coords$threshold
# 计算AUC和置信区间
auc_value <- auc(roc.list)
ci <- ci(roc.list, conf.level = 0.95) # 指定置信水平为95%

# 打印整个ci对象以便检查
print(ci)
print(auc_value)
# 从ci对象中提取95%置信区间的值
lower.ci <- ci$conf.int[1]
upper.ci <- ci$conf.int[2]

# 将AUC值和置信区间转换为文本形式，用于在图上显示
# 使用更多的小数位来确保精确度
auc_text <- sprintf("AUC = %.3f [%.3f, %.3f]", auc_value, lower.ci, upper.ci)

# 打印结果
print(auc_text)
auc_text<-"AUC=0.819 [0.731-0.906]"
# 使用ggplot2创建ROC曲线图，并添加AUC值和置信区间
g.list <- ggroc(roc.list, alpha = 1, size = 0.8, legacy.axes = TRUE, color="red") +
  theme_classic() +
  
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed", color = "blue") +
  annotate("text", x = 0.6, y = 0.2, label = auc_text, color = "black", size = 5)  # AUC文本颜色改为黑色

# 打印图形
print(g.list)



roc.list <- roc(train$diagnosis, gd1)
roc.list
g.list <- ggroc(roc.list, alpha = 1 ,size = 0.8,legacy.axes = TRUE,color="red")
g.list+theme_classic2() + ggtitle("train")+annotate(geom = "segment", x = 0, y = 0, xend =1, yend = 1)


#验证集的ROC
roc.list <- roc(test$diagnosis, gd2)
roc.list
g.list <- ggroc(roc.list, alpha = 1 ,size = 0.8,legacy.axes = TRUE,color="skyblue")
g.list <- ggroc(roc.list, alpha = 1 ,size = 0.8,legacy.axes = TRUE,color="red")
g.list+theme_classic2() + ggtitle("test")+annotate(geom = "segment", x = 0, y = 0, xend =1, yend = 1)

# 计算ROC曲线
roc.list <- roc(test$diagnosis, gd2)

# 计算AUC和置信区间
auc_value <- auc(roc.list)
ci <- ci(roc.list, conf.level = 0.95) # 指定置信水平为95%

# 打印整个ci对象以便检查
print(ci)
print(auc_value)
# 从ci对象中提取95%置信区间的值
lower.ci <- ci$conf.int[1]
upper.ci <- ci$conf.int[2]

# 将AUC值和置信区间转换为文本形式，用于在图上显示
# 使用更多的小数位来确保精确度
auc_text <- sprintf("AUC = %.3f [%.3f, %.3f]", auc_value, lower.ci, upper.ci)

# 打印结果
print(auc_text)
auc_text<-"AUC=0.630 [0.436-0.824]"
# 使用ggplot2创建ROC曲线图，并添加AUC值和置信区间
g.list <- ggroc(roc.list, alpha = 1, size = 0.8, legacy.axes = TRUE, color="red") +
  theme_classic() +
  
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed", color = "blue") +
  annotate("text", x = 0.6, y = 0.2, label = auc_text, color = "black", size = 5)  # AUC文本颜色改为黑色

# 打印图形
print(g.list)


#新验证集ROC
#验证集的ROC
roc.list <- roc(remain_data$diagnosis, gd2)
roc.list
g.list <- ggroc(roc.list, alpha = 1 ,size = 0.8,legacy.axes = TRUE,color="skyblue")
g.list <- ggroc(roc.list, alpha = 1 ,size = 0.8,legacy.axes = TRUE,color="red")
g.list+theme_classic2() + ggtitle("test")+annotate(geom = "segment", x = 0, y = 0, xend =1, yend = 1)

# 计算ROC曲线
roc.list <- roc(remain_data$diagnosis, gd2)

# 计算AUC和置信区间
auc_value <- auc(roc.list)
ci <- ci(roc.list, conf.level = 0.95) # 指定置信水平为95%

# 打印整个ci对象以便检查
print(ci)
print(auc_value)
# 从ci对象中提取95%置信区间的值
lower.ci <- ci$conf.int[1]
upper.ci <- ci$conf.int[2]

# 将AUC值和置信区间转换为文本形式，用于在图上显示
# 使用更多的小数位来确保精确度
auc_text <- sprintf("AUC = %.3f [%.3f, %.3f]", auc_value, lower.ci, upper.ci)

# 打印结果
print(auc_text)
auc_text<-"AUC=0.768 [0.681-0.855]"
# 使用ggplot2创建ROC曲线图，并添加AUC值和置信区间
g.list <- ggroc(roc.list, alpha = 1, size = 0.8, legacy.axes = TRUE, color="red") +
  theme_classic() +
  
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed", color = "blue") +
  annotate("text", x = 0.6, y = 0.2, label = auc_text, color = "black", size = 5)  # AUC文本颜色改为黑色

# 打印图形
print(g.list)




roc.list <- roc(test$diagnosis, gd3)
roc.list
g.list <- ggroc(roc.list, alpha = 1 ,size = 0.8,legacy.axes = TRUE,color="skyblue")
g.list <- ggroc(roc.list, alpha = 1 ,size = 0.8,legacy.axes = TRUE,color="red")
g.list+theme_classic2() + ggtitle("test")+annotate(geom = "segment", x = 0, y = 0, xend =1, yend = 1)

###校正曲线
ca11 <- calibrate(fit1, cmethod="hare",method="boot", B=1000,
                  xlab = "Nomogram Predicted Risk", ylab = "Actual Risk")
plot(ca11,xlim=c(0,1.0),ylim=c(0,1.0),
     xlab = "Nomogram Predicted Risk", ylab = "Actual Risk")
#Hosmer–Lemeshow test
library(ResourceSelection)

# 确保预测概率是逻辑回归模型的预测概率
gd <- predict(fit1, newdata = fold_train, type = "fitted")


# 运行Hosmer-Lemeshow检验
hoslem.test(y = fold_train$diagnosis, x = gd, g = 2)


# 添加预测概率到训练集
fold_train$predicted_prob <- gd

# 导出训练集（包括预测概率）到CSV
write.csv(fold_train, file = "train_with_predictions.csv", row.names = FALSE)


