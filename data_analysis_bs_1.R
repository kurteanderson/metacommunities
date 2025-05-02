library(ggplot2)

rm(out.frame)

dir.str <- "G:/My Drive/Stream Metacommunities/Local Models/bs_temp_runs/"

data.list <- list.files(dir.str)

for (i in 1:length(data.list)) {
  
  data.str <- paste(dir.str, data.list[i], sep = "")
  
  if (exists("out.frame") == FALSE) {
    out.frame <- read.table(data.str, sep = ",", header = TRUE)
  } 
  else {
    out.frame <- rbind(out.frame,read.table(data.str, sep = ",", header = TRUE))
  }
  
}

out.frame <- out.frame[out.frame$community.size > 0,]
out.frame <- out.frame[complete.cases(out.frame),]

# ggplot(data = out.frame, aes(x = community.size, y = beta.diversity, color = r.in, shape = as.factor(temp))) +
#   geom_point(size = 2) +
#   theme_classic() +
#   scale_x_continuous(limits = c(100,900)) +
#   scale_y_continuous(limits = c(0.1, 0.4)) +
#   labs(x = "Mean Community Size", y = "Monthly Temporal Beta Diversity", colour = "Recruitment\nProbability") +
#   theme(legend.position = "inside", legend.position.inside = c(.99, .99), legend.justification = c("right", "top"), legend.box.just = "left") +
#   scale_color_continuous(type = "viridis") +
#   scale_shape_manual(NULL, values = c(19,3), labels = c("15 °C", "25 °C"))
# 
# ggplot(data = out.frame, aes(x = community.size, y = ss.slope, color = r.in, shape = as.factor(temp))) +
#   geom_point(size = 2) +
#   theme_classic() +
#  # scale_x_continuous(limits = c(100,900)) +
# #  scale_y_continuous(limits = c(0.1, 0.4)) +
#   labs(x = "Mean Community Size", y = "Median Size Spectra Slope", colour = "Recruitment\nProbability") +
#   theme(legend.position = "inside", legend.position.inside = c(.99, .99), legend.justification = c("right", "top"), legend.box.just = "left") +
#   scale_color_continuous(type = "viridis") +
#   scale_shape_manual(NULL, values = c(19,3), labels = c("15 °C", "25 °C"))
# 
# ggplot(data = out.frame, aes(x = community.size, y = ss.range, color = r.in, shape = as.factor(temp))) +
#   geom_point(size = 2) +
#   theme_classic() +
#   # scale_x_continuous(limits = c(100,900)) +
#   #  scale_y_continuous(limits = c(0.1, 0.4)) +
#   labs(x = "Mean Community Size", y = "Median Size Spectra Slope", colour = "Recruitment\nProbability") +
#   theme(legend.position = "inside", legend.position.inside = c(.99, .99), legend.justification = c("right", "top"), legend.box.just = "left") +
#   scale_color_continuous(type = "viridis") +
#   scale_shape_manual(NULL, values = c(19,3), labels = c("15 °C", "25 °C"))

ggplot(data = out.frame, aes(x = as.factor(r.in), y = ss.slope, color = as.factor(temp), fill = as.factor(temp))) +
  geom_boxplot() +
  theme_classic() +
  labs(title = "c)", x = "Recruitment Probability", y = expression(paste("Median  ", lambda))) +
  theme(legend.position = "inside", legend.position.inside = c(.99, .99), legend.justification = c("right", "top"), legend.box.just = "left", plot.title = element_text(size = 18), axis.text = element_text(size = 16), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  scale_fill_manual(NULL, values = c("#0072B2", "#E69F00"), labels = c("15 °C", "25 °C")) +
  scale_color_manual(NULL, values = c("black", "black"), labels = c("15 °C", "25 °C")) 

ggplot(data = out.frame, aes(x = as.factor(r.in), y = ss.range, color = as.factor(temp), fill = as.factor(temp))) +
  geom_boxplot() +
  theme_classic() +
  labs(title = "d)", x = "Recruitment Probability", y = expression(paste("Range  ", lambda))) +
  theme(legend.position = "inside", legend.position.inside = c(.99, .99), legend.justification = c("right", "top"), legend.box.just = "left", plot.title = element_text(size = 18), axis.text = element_text(size = 16), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  scale_fill_manual(NULL, values = c("#0072B2", "#E69F00"), labels = c("15 °C", "25 °C")) +
  scale_color_manual(NULL, values = c("black", "black"), labels = c("15 °C", "25 °C")) 
  