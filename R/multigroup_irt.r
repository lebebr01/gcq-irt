# load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(mirt)
library(readxl)

# Function for DIF Extraction
get.dif.items <- function(out.list,p.val,parms){
  dif.items <- NULL
  non.dif.items <- NULL
  for(i in 1:length(out.list)){ # loop list of items
    chi.sq <- out.list[[i]][2,6] #2 groups, so second row, and 6th column
    df <- out.list[[i]][2,7]
    p <- out.list[[i]][2,8]
    i.name <- names(out.list[i])
    d <- c(i.name,chi.sq,df,p,parms[i,])
    
    if(p < p.val){
      dif.items <- rbind(dif.items,d)
    }else{
      non.dif.items <- rbind(non.dif.items,d)
    }
  }
  if (!is.null(dif.items)) {
    dif.items <- data.frame(dif.items, row.names = NULL)
    colnames(dif.items)[1] <- "item"
    colnames(dif.items)[2] <- "chi.sq"
    colnames(dif.items)[3] <- "df"
    colnames(dif.items)[4] <- "p"
    
  }
  if (!is.null(non.dif.items)) {
    non.dif.items <- data.frame(non.dif.items, row.names = NULL)
    colnames(non.dif.items)[1] <- "item"
    colnames(non.dif.items)[2] <- "chi.sq"
    colnames(non.dif.items)[3] <- "df"
    colnames(non.dif.items)[4] <- "p"
  }
  r.list <- list(dif_items = dif.items, no_dif = non.dif.items)
  return(r.list)
}

# load data
belin <- read_excel('Data/IExcel_IRT_All5A.xlsx')


# Items names (column names of items)
english_items <- dplyr::select(belin, starts_with('Eng_AV_'))
math_items <- dplyr::select(belin, starts_with('Math_AV_'))
read_items <- dplyr::select(belin, starts_with('Read_AV_'))
sci_items <- dplyr::select(belin, starts_with("Sci_AV_"))

# group variable
grade <- as.character(belin$Grade)


# Eigenvalue
english_eigen <- eigen(cor(english_items))$values
english_eigen[1] / english_eigen[2]
english_eigen[1] / sum(english_eigen)

math_eigen <- eigen(cor(math_items))$values
math_eigen[1] / math_eigen[2]
math_eigen[1] / sum(math_eigen)

read_eigen <- eigen(cor(read_items))$values
read_eigen[1] / read_eigen[2]
read_eigen[1] / sum(read_eigen)

sci_eigen <- eigen(cor(sci_items))$values
sci_eigen[1] / sci_eigen[2]
sci_eigen[1] / sum(sci_eigen)


# CTT DIF ----
library(difR)

belin_mh <- cbind(math_items, select(belin, Grade))
names = c('5', '6')
difGMH(belin_mh, group = "Grade", focal.name = names,
       p.adjust.method = "BH")

difGenLogistic(belin_mh, group = "Grade", focal.name = names,
               p.adjust.method = "BH")
difGenLogistic(belin_mh, group = "Grade", focal.name = names,
               p.adjust.method = "BH", type = 'nudif')
difGenLogistic(belin_mh, group = "Grade", focal.name = names,
               p.adjust.method = "BH", type = 'udif')

# logistic regression ---
belin_mh <- belin_mh %>%
  mutate(total_score = rowSums(select(., starts_with("Math_AV_"))))

math_it2_log <- glm(Math_AV_2 ~ 1 + total_score + factor(Grade) + total_score:factor(Grade),
    data = belin_mh, family = binomial)
summary(math_it2_log)

#--------------------------------------------------------
# English -----------------------------------------------
# Fit multigroup IRT model grade is group.
mult_group_eng <- multipleGroup(data = english_items, model = 1, group = grade,
                            invariance = c(colnames(english_items), 
                                           'free_means', 'free_var')) # 132 iterations

constrained_parameters_eng <- coef(mult_group_eng, simplify = TRUE)[[1]][[1]]
dif_add_eng <- DIF(mult_group_eng, c('a1', 'd'), scheme = 'drop',
                seq_state = .05)

# dif_out_eng <- get.dif.items(out.list = dif_add_eng, p.val = .05,
#                          parms = constrained_parameters_eng)
# dif_out_eng

# Freely estimate items with DIF.
model_eng <- '
 I = 1-40
 CONSTRAINB = (1-13, 15-20, 22, 23, 25, 28, 30-36, 39, a1), 
    (1-13, 15-20, 22, 23, 25, 28, 30-36, 39, d)
'
mult_group_part_eng <- multipleGroup(data = english_items, model = model_eng, 
                                     group = grade,
                                 invariance = c('free_means', 'free_var'), 
                                 #constrain = constrain, 
                                 technical = list(NCYCLES = 2000)) # 135 Iterations
parameters_eng <- coef(mult_group_part_eng, simplify = TRUE)

#--------------------------------------------------------
# Mathematics -----------------------------------------------
# Fit multigroup IRT model grade is group.
mult_group_math <- multipleGroup(data = math_items, model = 1, group = grade,
                            invariance = c(colnames(math_items), 
                                           'free_means', 'free_var')) # 186 iterations

constrained_parameters_math <- coef(mult_group_math, simplify = TRUE)[[1]][[1]]
dif_add_math <- DIF(mult_group_math, c('a1', 'd'), scheme = 'drop',
               seq_state = .05)
dif_add_math

# dif_out_math <- get.dif.items(out.list = dif_add_math, p.val = .05,
#                          parms = constrained.parameters)
# dif_out_math

# Freely estimate items with DIF.
model_math <- '
 I = 1-30
 CONSTRAINB = (1, 4, 5, 8, 11, 13, 14, 17-18, 21-23, 25-30, a1), 
    (1, 4, 5, 8, 11, 13, 14, 17-18, 21-23, 25-30, d)
'
mult_group_part_math <- multipleGroup(data = math_items, model = model_math, 
                                      group = grade,
                                 invariance = c('free_means', 'free_var'), 
                                 #constrain = constrain, 
                                 technical = list(NCYCLES = 2000)) # 243 Iterations
parameters_math <- coef(mult_group_part_math, simplify = TRUE)

# probability ----
prob_item2 <- data.frame(rbind(
  cbind(probtrace(extract.item(mult_group_part_math, 2, group = 1), Theta = seq(-6, 6, .1)), theta = seq(-6, 6, .1), group = 'Grade 4'),
                         cbind(probtrace(extract.item(mult_group_part_math, 2, group = 2), Theta = seq(-6, 6, .1)), theta = seq(-6, 6, .1), group = 'Grade 5'),
                         cbind(probtrace(extract.item(mult_group_part_math, 2, group = 3), Theta = seq(-6, 6, .1)), theta = seq(-6, 6, .1), group = 'Grade 6')
                         )
                         )
prob_item2$P.0 <- as.numeric(as.character(prob_item2$P.0))
prob_item2$P.1 <- as.numeric(as.character(prob_item2$P.1))
prob_item2$theta <- as.numeric(as.character(prob_item2$theta))

p <- ggplot(prob_item2, aes(x = theta, y = P.1)) + 
  geom_line(aes(group = group, color = group), size = 2) + 
  scale_x_continuous("Aptitude", breaks = seq(-6, 6, 1)) + 
  ylab("Probability Correct") + 
  scale_color_grey("", start = 0) + 
  theme_bw(base_size = 12)
ggsave(filename = "icc_math2.png", plot = p, dpi = 'retina', height = 4, width = 6,
       path = 'paper/figs')

# iteminformation ----
info_item2 <- data.frame(rbind(
  cbind(iteminfo(extract.item(mult_group_part_math, 2, group = 1), Theta = seq(-6, 6, .1)), theta = seq(-6, 6, .1), group = 'Grade 4'),
  cbind(iteminfo(extract.item(mult_group_part_math, 2, group = 2), Theta = seq(-6, 6, .1)), theta = seq(-6, 6, .1), group = 'Grade 5'),
  cbind(iteminfo(extract.item(mult_group_part_math, 2, group = 3), Theta = seq(-6, 6, .1)), theta = seq(-6, 6, .1), group = 'Grade 6')
)
)
info_item2$V1 <- as.numeric(as.character(info_item2$V1))
info_item2$theta <- as.numeric(as.character(info_item2$theta))

p <- ggplot(info_item2, aes(x = theta, y = V1)) + 
  geom_line(aes(group = group, color = group), size = 2) + 
  scale_x_continuous("Aptitude", breaks = seq(-6, 6, 1)) + 
  ylab("Item Information") + 
  viridis::scale_color_viridis("", discrete = TRUE, option = 'magma') + 
  theme_bw(base_size = 12)
ggsave(filename = "info_math2.png", plot = p, dpi = 'retina', height = 4, width = 6,
       path = 'paper/figs')

# Model Implied test score ----
# test_scores <- data.frame(rbind(
#   cbind(expected.test(extract.group(mult_group_part_math, 1), Theta = seq(-6, 6, .1)), theta = seq(-6, 6, .1), group = 'Grade 4'),
#   cbind(expected.test(extract.group(mult_group_part_math, 2), Theta = seq(-6, 6, .1)), theta = seq(-6, 6, .1), group = 'Grade 5'),
#   cbind(expected.test(extract.group(mult_group_part_math, 3), Theta = seq(-6, 6, .1)), theta = seq(-6, 6, .1), group = 'Grade 6')
# )
# )

tmp <- plot(mult_group_part_math)

# tmp$panel.args
pltdata <- data.frame(lapply(tmp$panel.args, function(x) do.call(cbind, x))[[1]])
names(pltdata) <- c("theta", "tcc", "num")
pltdata$group <- paste0('Grade ', tmp$panel.args.common$groups)

p <- ggplot(pltdata, aes(x = theta, y = tcc)) + 
  geom_line(aes(group = group, color = group), size = 2) + 
  scale_x_continuous("Aptitude", breaks = seq(-6, 6, 1)) + 
  scale_y_continuous("Model Implied Number Correct", breaks = seq(0, 30, 5),
                     limits = c(0, 30)) + 
  scale_color_grey("") + 
  theme_bw(base_size = 12)
ggsave(filename = "test_math.png", plot = p, dpi = 'retina', height = 4, width = 6,
       path = 'paper/figs')

# Using plink to get model implied test scores ----
library(plink)

parameters_math <- coef(mult_group_part_math, simplify = TRUE, 
                        IRTpars = TRUE)

item_info <- bind_rows(
  cbind(group = 'Grade 4', drm(data.frame(parameters_math[[1]]$items[, 1:2]), seq(-6, 6, .1))@prob
        ),
  cbind(group = 'Grade 5', drm(data.frame(parameters_math[[2]]$items[, 1:2]), seq(-6, 6, .1))@prob
  ),
  cbind(group = 'Grade 6', drm(data.frame(parameters_math[[3]]$items[, 1:2]), seq(-6, 6, .1))@prob
  )
) %>%
  mutate(tcc = rowSums(select(., contains('item'))))

p <- ggplot(item_info, aes(x = theta1, y = tcc)) + 
  geom_line(aes(group = group, color = group), size = 2) + 
  scale_x_continuous("Aptitude", breaks = seq(-6, 6, 1)) + 
  scale_y_continuous("Model Implied Number Correct", breaks = seq(0, 30, 5),
                     limits = c(0, 30)) + 
  scale_color_grey("", start = 0) + 
  theme_bw(base_size = 12)
ggsave(filename = "test_math.png", plot = p, dpi = 'retina', height = 4, width = 6,
       path = 'paper/figs')

#--------------------------------------------------------
# Reading -----------------------------------------------
# Fit multigroup IRT model grade is group.
mult_group_read <- multipleGroup(data = read_items, model = 1, group = grade,
                                 invariance = c(colnames(read_items), 
                                                'free_means', 'free_var')) # 144 iterations

constrained_parameters_read <- coef(mult_group_read, simplify = TRUE)[[1]][[1]]
dif_add_read <- DIF(mult_group_read, c('a1', 'd'), scheme = 'drop',
                    seq_state = .05)
dif_add_read

# dif_out_read <- get.dif.items(out.list = dif_add_read, p.val = .05,
#                               parms = constrained.parameters)
# dif_out_read

# Freely estimate items with DIF.
model_read <- '
  I = 1-30
  CONSTRAINB = (2-18, 21, 24-28, 30, a1), 
    (2-18, 21, 24-28, 30, d)
'
mult_group_part_read <- multipleGroup(data = read_items, model = model_read, 
                                      group = grade,
                                      invariance = c('free_means', 'free_var'), 
                                      #constrain = constrain, 
                                      technical = list(NCYCLES = 2000)) # 144 Iterations
parameters_read <- coef(mult_group_part_read, simplify = TRUE)

#--------------------------------------------------------
# Science -----------------------------------------------
# Fit multigroup IRT model grade is group.
mult_group_sci <- multipleGroup(data = sci_items, model = 1, group = grade,
                                 invariance = c(colnames(sci_items), 
                                                'free_means', 'free_var')) # 171 iterations

constrained_parameters_sci <- coef(mult_group_sci, simplify = TRUE)[[1]][[1]]
dif_add_sci <- DIF(mult_group_sci, c('a1', 'd'), scheme = 'drop',
                    seq_state = .05)
dif_add_sci

# dif_out_sci <- get.dif.items(out.list = dif_add_sci, p.val = .05,
#                               parms = constrained.parameters)
# dif_out_sci

# Freely estimate items with DIF.
model_sci <- '
  I = 1-28
  CONSTRAINB = (2-14, 16-21, 23-28, a1), 
    (2-14, 16-21, 23-28, d)
'
mult_group_part_sci <- multipleGroup(data = sci_items, model = model_sci, 
                                      group = grade,
                                      invariance = c('free_means', 'free_var'), 
                                      #constrain = constrain, 
                                      technical = list(NCYCLES = 2000)) # 122 Iterations
parameters_sci <- coef(mult_group_part_sci, simplify = TRUE)


# save models
save(mult_group_part_eng, mult_group_part_math, mult_group_part_read, mult_group_part_sci,
     file = 'Data/multiple_group_models.rda')

M2(mult_group_part_eng)
M2(mult_group_part_math)
M2(mult_group_part_read)
M2(mult_group_part_sci)


# Sex DIF ----
sex <- as.character(belin$ngender)


# English -----------------------------------------------
# Fit multigroup IRT model sex is group.
mult_group_eng <- multipleGroup(data = english_items, model = 1, group = sex,
                                invariance = c(colnames(english_items), 
                                               'free_means', 'free_var')) # 30 iterations

constrained_parameters_eng <- coef(mult_group_eng, simplify = TRUE)[[1]][[1]]
dif_add_eng <- DIF(mult_group_eng, c('a1', 'd'), scheme = 'drop',
                   seq_state = .05)

# dif_out_eng <- get.dif.items(out.list = dif_add_eng, p.val = .05,
#                          parms = constrained_parameters_eng)
# dif_out_eng

# Freely estimate items with DIF.
model_eng <- '
I = 1-40
CONSTRAINB = (1-8, 10, 12, 14-16, 18-30, 32, 34, 36-38, 40, a1), 
(1-8, 10, 12, 14-16, 18-30, 32, 34, 36-38, 40, d)
'
mult_group_part_eng <- multipleGroup(data = english_items, model = model_eng, 
                                     group = sex,
                                     invariance = c('free_means', 'free_var'), 
                                     #constrain = constrain, 
                                     technical = list(NCYCLES = 2000)) # 38 Iterations
parameters_eng <- coef(mult_group_part_eng, simplify = TRUE)

# Mathematics -----------------------------------------------
# Fit multigroup IRT model grade is group.
mult_group_math <- multipleGroup(data = math_items, model = 1, group = sex,
                                 invariance = c(colnames(math_items), 
                                                'free_means', 'free_var')) # 33 iterations

constrained_parameters_math <- coef(mult_group_math, simplify = TRUE)[[1]][[1]]
dif_add_math <- DIF(mult_group_math, c('a1', 'd'), scheme = 'drop',
                    seq_state = .05)
dif_add_math

# dif_out_math <- get.dif.items(out.list = dif_add_math, p.val = .05,
#                          parms = constrained.parameters)
# dif_out_math

# Freely estimate items with DIF.
model_math <- '
I = 1-30
CONSTRAINB = (1-8, 10-12, 16-20, 22, 23, 25-28, 30, a1), 
(1-8, 10-12, 16-20, 22, 23, 25-28, 30, d)
'
mult_group_part_math <- multipleGroup(data = math_items, model = model_math, 
                                      group = sex,
                                      invariance = c('free_means', 'free_var'), 
                                      #constrain = constrain, 
                                      technical = list(NCYCLES = 2000)) # 19 Iterations
parameters_math <- coef(mult_group_part_math, simplify = TRUE)

#--------------------------------------------------------
# Reading -----------------------------------------------
# Fit multigroup IRT model grade is group.
mult_group_read <- multipleGroup(data = read_items, model = 1, group = sex,
                                 invariance = c(colnames(read_items), 
                                                'free_means', 'free_var')) # 50 iterations

constrained_parameters_read <- coef(mult_group_read, simplify = TRUE)[[1]][[1]]
dif_add_read <- DIF(mult_group_read, c('a1', 'd'), scheme = 'drop',
                    seq_state = .05)
dif_add_read

# dif_out_read <- get.dif.items(out.list = dif_add_read, p.val = .05,
#                               parms = constrained.parameters)
# dif_out_read

# Freely estimate items with DIF.
model_read <- '
  I = 1-30
  CONSTRAINB = (1, 4-9, 11-13, 15-20, 22, 23, 25-27, 29, a1), 
    (1, 4-9, 11-13, 15-20, 22, 23, 25-27, 29, d)
'
mult_group_part_read <- multipleGroup(data = read_items, model = model_read, 
                                      group = sex,
                                      invariance = c('free_means', 'free_var'), 
                                      #constrain = constrain, 
                                      technical = list(NCYCLES = 2000)) # 54 Iterations
parameters_read <- coef(mult_group_part_read, simplify = TRUE)

#--------------------------------------------------------
# Science -----------------------------------------------
# Fit multigroup IRT model grade is group.
mult_group_sci <- multipleGroup(data = sci_items, model = 1, group = sex,
                                invariance = c(colnames(sci_items), 
                                               'free_mean', 'free_var')) # 27 iterations

constrained_parameters_sci <- coef(mult_group_sci, simplify = TRUE)[[1]][[1]]
dif_add_sci <- DIF(mult_group_sci, c('a1', 'd'), scheme = 'drop',
                   seq_state = .05)
dif_add_sci

# dif_out_sci <- get.dif.items(out.list = dif_add_sci, p.val = .05,
#                               parms = constrained.parameters)
# dif_out_sci

# Freely estimate items with DIF.
model_sci <- '
  I = 1-28
  CONSTRAINB = (1-3, 5, 6, 9-10, 13, 15-18, 20-21, 23-27, a1), 
    (1-3, 5, 6, 9-10, 13, 15-18, 20-21, 23-27, d)
'
mult_group_part_sci <- multipleGroup(data = sci_items, model = model_sci, 
                                     group = sex,
                                     invariance = c('free_mean', 'free_var'), 
                                     #constrain = constrain, 
                                     technical = list(NCYCLES = 2000)) # 39 Iterations
parameters_sci <- coef(mult_group_part_sci, simplify = TRUE)
