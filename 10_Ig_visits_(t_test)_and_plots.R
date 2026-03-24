### Author Qamar Feddah
### Ig and Visit association

source("00_setup.R")

colData <- readRDS("outputs/colData.rds")

#IgG
t.test(IgG ~ Visit, data = colData)
summary(lm(IgG ~ Visit, data = colData))

p0 <- ggplot(colData, aes(x = Visit, y = IgG)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_minimal() +
  labs(title = "IgG by Visit")

#IgG1
t.test(IgG1 ~ Visit, data = colData)
summary(lm(IgG1 ~ Visit, data = colData))

p1 <- ggplot(colData, aes(x = Visit, y = IgG1)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_minimal() +
  labs(title = "IgG1 by Visit")

#IgG2
t.test(IgG2 ~ Visit, data = colData)
summary(lm(IgG2 ~ Visit, data = colData))

p2 <- ggplot(colData, aes(x = Visit, y = IgG2)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_minimal() +
  labs(title = "IgG2 by Visit")

#IgG3
t.test(IgG3 ~ Visit, data = colData)
summary(lm(IgG3 ~ Visit, data = colData))

p3 <- ggplot(colData, aes(x = Visit, y = IgG3)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_minimal() +
  labs(title = "IgG3 by Visit")

#IgG4
t.test(IgG4 ~ Visit, data = colData)
summary(lm(IgG4 ~ Visit, data = colData))

p4 <- ggplot(colData, aes(x = Visit, y = IgG4)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_minimal() +
  labs(title = "IgG4 by Visit")

#IgA
t.test(IgA ~ Visit, data = colData)
summary(lm(IgA ~ Visit, data = colData))

p5 <- ggplot(colData, aes(x = Visit, y = IgA)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_minimal() +
  labs(title = "IgA by Visit")


pdf("outputs/plots/Ig_by_Visit.pdf", width = 7, height = 6)

print(p0)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)  

dev.off()

# Save t-test and lm outputs into one text file
sink("outputs/T_test/Ig_by_Visit_tests.txt")

cat("IgG\n")
print(t.test(IgG ~ Visit, data = colData))
print(summary(lm(IgG ~ Visit, data = colData)))

cat("\nIgG1\n")
print(t.test(IgG1 ~ Visit, data = colData))
print(summary(lm(IgG1 ~ Visit, data = colData)))

cat("\nIgG2\n")
print(t.test(IgG2 ~ Visit, data = colData))
print(summary(lm(IgG2 ~ Visit, data = colData)))

cat("\nIgG3\n")
print(t.test(IgG3 ~ Visit, data = colData))
print(summary(lm(IgG3 ~ Visit, data = colData)))

cat("\nIgG4\n")
print(t.test(IgG4 ~ Visit, data = colData))
print(summary(lm(IgG4 ~ Visit, data = colData)))

cat("\nIgA\n")
print(t.test(IgA ~ Visit, data = colData))
print(summary(lm(IgA ~ Visit, data = colData)))

sink()

# Save full objects as RDS for later reuse
results <- list(
  IgG  = list(t = t.test(IgG ~ Visit, data = colData),
              lm = lm(IgG ~ Visit, data = colData)),
  IgG1 = list(t = t.test(IgG1 ~ Visit, data = colData),
              lm = lm(IgG1 ~ Visit, data = colData)),
  IgG2 = list(t = t.test(IgG2 ~ Visit, data = colData),
              lm = lm(IgG2 ~ Visit, data = colData)),
  IgG3 = list(t = t.test(IgG3 ~ Visit, data = colData),
              lm = lm(IgG3 ~ Visit, data = colData)),
  IgG4 = list(t = t.test(IgG4 ~ Visit, data = colData),
              lm = lm(IgG4 ~ Visit, data = colData)),
  IgA  = list(t = t.test(IgA ~ Visit, data = colData),
              lm = lm(IgA ~ Visit, data = colData))
)

