arg <- commandArgs(TRUE)

file <- read.csv(arg[1], header = FALSE, sep = " " )
reported_ancestry <- read.csv("1KG_Phenos_With_GIH.txt", sep = " ", header = FALSE)
ancestry <- c("ASW", "CEU", "CHB", "CHS", "CLM", "FIN", "GBR", "IBS", "JPT", "LWK", "MXL", "PUR", "TSI", "YRI", "GIH")

dat <- merge(file, reported_ancestry, by = "V2", all = TRUE)

dat[14] <- data.frame(dat[14] %>%
                      map(lfactor, levels = c(3,4,5,6,7,8,10,11,12,13,14,15,16,17,18), labels = ancestry ))


dat$ancestry <- as.character(dat$V3.y)



#interactive map


#pc 1 and 2
library(plotly)


p1and2 <- ggplot(data = dat,aes(x = V3.x, y = V4, color = ancestry)) +
  geom_point () +
  scale_x_continuous(breaks = c(-0.040,-0.035, -0.030,-0.025,-0.021,-0.020,-0.015,-0.010, -0.005, -0.0035, 0.000, 0.005, 0.010, 0.015, 0.020, 0.023,0.025, 0.030, 0.035, 0.040)) +
  scale_y_continuous(breaks = c(-0.040,-0.020,0.000,0.013,0.020,0.040,0.060,0.080,0.100,0.120,0.140,0.160)) +
  labs(x = "PC1", y = "PC2")
fig1and2 <- ggplotly(p1and2)
fig1and2


#pc 1 and 3
p1and2 <- ggplot(data = dat,aes(x = V3.x, y = V5, color = ancestry)) +
  geom_point () +
  scale_x_continuous(breaks = c(-0.040,-0.035, -0.030,-0.025,-0.020,-0.015,-0.010,-0.005,0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040)) +
  scale_y_continuous(breaks = c(-0.020,0.000,0.020,0.040,0.060,0.080,0.100,0.120,0.140,0.160)) +
  labs(x = "PC1", y = "PC3")
fig1and2 <- ggplotly(p1and2)
fig1and2


#pc 2 and 3
p2and3 <- ggplot(dat,aes(x = V4, y = V5, color = ancestry)) +
  geom_point () +
  scale_x_continuous(breaks = c(-0.050, -0.045,-0.040,-0.035, -0.030,-0.025,-0.020,-0.015,-0.010,-0.005,0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040)) +
  scale_y_continuous(breaks = c(-0.020,0.000,0.020,0.040,0.060,0.080,0.100,0.120,0.140,0.160)) +
  labs(x = "PC2", y = "PC3")
fig2and3 <- ggplotly(p2and3)
fig2and3


#pc 3 and 4
p3and4 <- ggplot(dat,aes(x = V5, y = V6, color = ancestry)) +
  geom_point () +
  scale_x_continuous(breaks = c(-0.040,-0.035, -0.030,-0.020,-0.010,0.000, 0.010, 0.020, 0.030, 0.040, 0.050, 0.060)) +
  scale_y_continuous(breaks = c(-0.020,0.000,0.020,0.040,0.060,0.080,0.100,0.120,0.140,0.160)) +
  labs(x = "PC3", y = "PC4")
fig3and4 <- ggplotly(p3and4)
fig3and4
