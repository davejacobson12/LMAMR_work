#R intro tutorial (using ggplot2 for boxplots)

#Working with an alpha diversity file, with metadata included, from Qiime (output of add_alpha_to_mapping.py). We will be working with this file in R (I prefer to use R studio - any version of R should be ok for our purposes)

#Using nano, or another text editor, remove # from #SampleID. R ignores all lines that begin with # and we don't want this to happen, so we must remove the # from #SampleID. 

#Additionally, it is best if the sample names do not start with a numeric character. R will change the sample names if they start with a number (usually will prepend an x. to the start of the sample name, I think) but this can cause problems when merging files. It is easier to fix this problem before loading the file into R. 

#Loading the file into R. This code is operating under the assumption that your sample names are in the first column and the first row has all of your metadata categories. 

#Open R studio

#read in the file 

myFile <- read.table("path to file", header =T, row.names =1, sep = "\t")

#read.table is the command to read your file and will make your data frame (myFile)
#always provide the full path to your file
#header =T is saying that you want your first row to be the column headers. I typically use header = T but you may not want to in some circumstances
#row.names =1 is saying that the first column are you sample names. You don't need to include this, but I prefer to. It makes things easier for merging files later on (in my opinion). If your sample names are in row 2 or 3 or 4, you would put row.names =2, or row.names =3, etc. 
#sep = "\t" tells R that your file is tab separated. If it is comma separated or something else, you would use "," or however your file is separated

#Look at your file

View(myFile)

#Make a basic boxplot in R to get a first glance at the data

boxplot(myFile$observed_otus ~ myFile$metadata, ylab = "Y axis Title", xlab = "X axis Title", main = "Plot Title", las =2)

#when specifying variables, give the name of the data frame followed by $ and name of variable column
#the continuous variable goes first, whether its observed OTUs or any other variable
# ~ is necessary to make a formula
#the metadata category comes second
#ylab allows you to label the y axis, otherwise it will be blank
#xlab allows you to label the x axis, otherwise it will be blank
#main allows you to label the title, otherwise it will be blank
# las =2 makes the x axis tick mark labels to be horizontal. This is useful when the metadata has long variable names

#There are many other options, but this is just to look at the data. We'll make a nicer boxplot in ggplot2
#You can make this boxplot for whichever metadata categories you're interested in. If there are a couple that look interesting, you can perform a statistical test

#Kruskal-Wallis test is a non-parametric ANOVA test. So you can test a continuous variable against a metadata variable with more than 2 categories. You can also do a one-way ANOVA

kruskal.test(myFile$observed_otus ~ myFile$metadata)

#if you have a significant p-value, you can perform a Dunn post hoc test to determine which pairwise comparisons are significantly different. 

#install and load FSA
install.packages("FSA")
library(FSA)

dunnTest(myFile$observed_otus ~ myFile$metadata, method = "bh")

#The formula is the same as used for the kruskal.test
#I prefer to use the Benjamini-Hochberg False Discover Rate test when using adjusted p-values for multiple comparisons.

#install and load the GGplot2

install.packages("ggplot2")
library(ggplot2)

#Assing ggplot object (this can be used to make a variety of different graphs)

myGG <- ggplot(data = myFile, mapping = aes(x = myFile$age_cat y = myFile$observed_otus, fill = myFile$other_metadata)))

myGG + geom_boxplot()

#This is the basic code. Notice that you use x= and y= for your categorical and continuous variables, respectively
#fill is used to add coloring based on metadata categories. You can use the same variable as your x =, and this will result in a boxplot that has a different color for each x tick mark.

#You can add many to the plot from this based code. 

#to manually select the colors

myGG + geom_boxplot() + scale_fill_manual(values = c("red", "blue", "whatever colors you want")
                                          
                                          #I usually go to colorbrewer2.org and pick my colors based on their scheme
                                          
                                          #To add labels, change their size, make them bold 
                                          
                                          myGG + geom_boxplot() + xlab("Metadata") + ylab("OTUs") + theme(axis.title.x = element_text(face = "bold", size = 15)
                                                                                                          
                                                                                                          #To center the title
                                                                                                          
                                                                                                          myGG + geom_boxplot() + ggtitle("Title") + theme(plot.title = element_text(hjust = 0.5)
                                                                                                                                                           
                                                                                                                                                           #To rename the legend title
                                                                                                                                                           myGG + geom_boxplot() + labs(fill = "legend title")
                                                                                                                                                           
                                                                                                                                                           
                                                                                                                                                           #there are limitless options to altering the plot, but this is a good place to start
                                                                                                                                                           
                                                                                                                                                           
                                                                                                                                                           