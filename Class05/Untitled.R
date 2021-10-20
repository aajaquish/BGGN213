#Class 05: Data Visualization 10-13-21
#Today we are going to use ggplot2 package

#first you need to load the package/library before you can use it
library(ggplot2)
ggplot(cars)
cars

#We will use this inbuilt "cars" dataset first
head(cars)

#All ggplots have at least 3 layers
# data + aes (aesthetics) + geoms (geometries)
#labs = different labels for features
#geom_line, geom_smooth (method="lm")to make it linear, grey is STerror
ggplot(data=cars) + 
  aes(x=speed, y=dist) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(title="Stopping Distance of Old Cars", x= "Speed (MPH)",
       y="Stopping Distance (feet)")

#Side-note: ggplot is not the only graphics system
# a very popular one is good old "base" R graphics
plot(cars)


#Lab Project Week4
#url <- "https://bioboot.github.io/bimm143_S20/class-material/up_down_expression.txt"
genes <- read.delim(url)
head(genes)

#Q. How many genes in the dataset?
nrow(genes)

#What are the column names and how many there are?
colnames(genes)
ncol(genes)

#How many upregulated genes in the State column?
table(genes$State)

#What fraction of total genes is upregulated?
round(table(genes$State)/nrow(genes)*100, 2)

#Use the genes dataset in ggplot function
#set the x and y aesthetic mappings to the Condition1 and Condition2
p <- ggplot(genes) +
  aes(x=Condition1, y=Condition2, col=State) +
  geom_point()

#Adding new colors
p + scale_colour_manual( values=c("blue", "gray", "red"))


p + geom_point(col="blue")
p + aes(col=State) + geom_point(col="blue")
p + geom_point(col="blue") + aes(col=State)

#Labeling x and yand title
p + scale_colour_manual( values=c("blue", "gray", "red")) + 
  labs(title= "Gene Expression Changes Upon Drug Treatment", x="Control (no drug)", y="Drug Treatment")

#6. Optional--Let's explore gapminder (a dataset)
#install.packages("gapminder")
library(gapminder)
# File location online
#url <- "https://raw.githubusercontent.com/jennybc/gapminder/master/inst/extdata/gapminder.tsv"

#gapminder <- read.delim(url)
library(gapminder)
head(gapminder)

#Lets make a new plot of year vs lifeexp
#alpha=transparency
ggplot(gapminder) +
  aes(x=year, y=lifeExp, col=continent) +
  geom_jitter(width=0.3,alpha=0.4) +
  geom_boxplot(aes(group=year), alpha=0.2)

ggplot(gapminder) +
  aes(x=year, y=lifeExp, col=continent) +
  geom_jitter(width=0.3,alpha=0.4) +
  geom_violin(aes(group=year), alpha=0.2)

ggplot(gapminder) +
  aes(x=year, y=lifeExp, col=continent) +
  geom_jitter(width=0.3,alpha=0.4) 
  #geom_violin(aes(group=year), alpha=0.2,
              #draw_quantiles=0.5

#Install the plotly
#install.packages("plotly")

ggplotly()

#Plot Animation
#install.packages("gifski")
#install.packages("gganimate")
