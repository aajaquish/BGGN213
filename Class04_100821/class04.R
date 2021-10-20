x <- 1:50
plot(x, sin(x), type="l", col="blue", lwd="3", xlab="My nice text")

# install.packages("tinytex")
# tinytex::install_tinytex()
x <- rnorm(1000)
mean(x)
sd(x)

summary(x)
boxplot(x)
hist(x)
rug(x)
