### Author Qamar Feddah
### Ig and Sex association

set.seed(333)
source("00_setup.R")

colData <- readRDS("outputs/colData.rds")

#IgG
t.test(IgG ~ Sex, data = colData)
summary(lm(IgG ~ Sex, data = colData))

p0 <- ggplot(colData, aes(x = Sex, y = IgG)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_minimal() +
  labs(title = "IgG by Sex")

#IgG1
t.test(IgG1 ~ Sex, data = colData)
summary(lm(IgG1 ~ Sex, data = colData))

p1 <- ggplot(colData, aes(x = Sex, y = IgG1)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_minimal() +
  labs(title = "IgG1 by Sex")

#IgG2
t.test(IgG2 ~ Sex, data = colData)
summary(lm(IgG2 ~ Sex, data = colData))

p2 <- ggplot(colData, aes(x = Sex, y = IgG2)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_minimal() +
  labs(title = "IgG2 by Sex")

#IgG3
t.test(IgG3 ~ Sex, data = colData)
summary(lm(IgG3 ~ Sex, data = colData))

p3 <- ggplot(colData, aes(x = Sex, y = IgG3)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_minimal() +
  labs(title = "IgG3 by Sex")

#IgG4
t.test(IgG4 ~ Sex, data = colData)
summary(lm(IgG4 ~ Sex, data = colData))

p4 <- ggplot(colData, aes(x = Sex, y = IgG4)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_minimal() +
  labs(title = "IgG4 by Sex")

#IgA
t.test(IgA ~ Sex, data = colData)
summary(lm(IgA ~ Sex, data = colData))

p5 <- ggplot(colData, aes(x = Sex, y = IgA)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_minimal() +
  labs(title = "IgA by Sex")


pdf("outputs/plots/Ig_by_Sex.pdf", width = 7, height = 6)

print(p0)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)  

dev.off()

# Save t-test and lm outputs into one text file
dir.create("outputs/T_test")
sink("outputs/T_test/Ig_by_Sex_T_tests.txt")

cat("IgG\n")
print(t.test(IgG ~ Sex, data = colData))
print(summary(lm(IgG ~ Sex, data = colData)))

cat("\nIgG1\n")
print(t.test(IgG1 ~ Sex, data = colData))
print(summary(lm(IgG1 ~ Sex, data = colData)))

cat("\nIgG2\n")
print(t.test(IgG2 ~ Sex, data = colData))
print(summary(lm(IgG2 ~ Sex, data = colData)))

cat("\nIgG3\n")
print(t.test(IgG3 ~ Sex, data = colData))
print(summary(lm(IgG3 ~ Sex, data = colData)))

cat("\nIgG4\n")
print(t.test(IgG4 ~ Sex, data = colData))
print(summary(lm(IgG4 ~ Sex, data = colData)))

cat("\nIgA\n")
print(t.test(IgA ~ Sex, data = colData))
print(summary(lm(IgA ~ Sex, data = colData)))

sink()

# Save full objects as RDS for later reuse
results <- list(
  IgG  = list(t = t.test(IgG ~ Sex, data = colData),
              lm = lm(IgG ~ Sex, data = colData)),
  IgG1 = list(t = t.test(IgG1 ~ Sex, data = colData),
              lm = lm(IgG1 ~ Sex, data = colData)),
  IgG2 = list(t = t.test(IgG2 ~ Sex, data = colData),
              lm = lm(IgG2 ~ Sex, data = colData)),
  IgG3 = list(t = t.test(IgG3 ~ Sex, data = colData),
              lm = lm(IgG3 ~ Sex, data = colData)),
  IgG4 = list(t = t.test(IgG4 ~ Sex, data = colData),
              lm = lm(IgG4 ~ Sex, data = colData)),
  IgA  = list(t = t.test(IgA ~ Sex, data = colData),
              lm = lm(IgA ~ Sex, data = colData))
)

