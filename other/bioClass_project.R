###Jacobson R code. This document has my R code and comments, the tables and plots are found in Jacobson_assignment2.pdf

### I am using a linear mixed effects model to approach the problem in assignment 2. This is because the data contains both fixed and random effects and you need to account for pseudoreplicaiton in this problem. Cohort is a random effect because it is only expected to influence the variance, as mice were randomly grouped into different cohorts. The remainder of the factors are fixed effects. We are interested in testing whether age has an effect on memory and learning. I will test the impact of age on memory by testing whether certain age group(s) of mice are able to remember the location of the reward longer (i.e. after more testing runs) than other age groups. My null hypothesis is that age group has no effect on mouse memory. Additionally, I will test age on learning by testing whether certain age groups of mice learn better (i.e. more mice go to the reward in testing phase at x number of trial runs). My null hypothesis is that age has no impact on mice and that the same number mice from each age group will go to the reward at x number of trial runs. 

>memory <- read.table(“memory.txt”, sep=“\t”, header = T) ### read the dataset into R

#### Must build the model before testing the assumptions

library(nlme) ### load the library

attach(memory)

>memory.lme <- lme(COUNT ~ TREATMENT * STRAIN * AGE * NUM_TRAIN * RUN, data = memory, random = ~ 1 | COHORT) ### Cohort is the only random effect, test all of the main effects and the interactions of the main effects

>summary(memory.lme) ### test the assumption that within group errors are centered around 0


## Random effects are independent in different groups and within group errors are independent of random effects is a good assumption because different the cohorts are all independent of each other and grouping of cohorts does not impact how the fixed effects impact the model.

> getVarCov(memory.lme, type = "random.effects") ## Random effects covariance matrix


> qqnorm(memory.lme, ~ranef(.)) ### random effects are normally distributed with a mean of 0. 



>anova(memory.lme) ### because my memory.lme accounts for all main effects and interactions, this anova is sufficient for both the memory and learning hypotheses. 


### Here we see that TREATMENT, TREATMENT:RUN, TREATMEN:NUM_TRAIN, and TREATMENT:RUN:NUM_TRAIN all have significant effects (p<0.0001). I will first address these results because they are not directly what we are testing but provide interesting information about the population. It is important to note that we are only interested in the results with treatment, because this is what tells us about whether the mice are ending up in the reward or control section. 

> tapply(COUNT, list(RUN,TREATMENT), mean) #### see the direction of effect of run on treatment. I am using the mean value to see what the average number of mice in control and reward for each comparison.

### output ### first column is number of runs, second column is mice in control section, third column is mice in reward section



### here we see that mice are more likely to find the reward after the first run and then find control and reward at equal rates after 15 runs. This shows the depreciation in memory over time in mice in general. Again this isn’t what we are testing but it orients us with how the mice perform in general.

> tapply(COUNT, list(NUM_TRAIN,TREATMENT), mean) ### direction of effect of number of trainings on treatment (learning)

###output### first column is the number of trainings, second column is mice in control section, third column is mice in reward section


### We see that mice with just 1 training choose reward and control equally, while mice that have more training runs are more likely to go to the reward section than the control section of the maze.

> ftable(tapply(COUNT, list(NUM_TRAIN,RUN,TREATMENT), mean))

###output### First column is number of trainings, second column is number of runs, third column is mice in control section, fourth column is mice in reward section



##### Without diving too far into the tables, we can see that all mice are more likely to go to the reward section in the first couple runs, yet the effect is more pronounced in those mice with more training. Those with more training also appear to have better memory, that is, it takes more trials for them to return to a state of being no more likely to choose the reward or control side. 


##### Now to test the effect that age has on memory and learning. The original anova showed us that TREATMENT:AGE, TREATMENT:AGE:RUN, TREATMENT:AGE:NUM_TRAIN, and TREATMENT:AGE:NUM_TRAIN:RUN all have a significant p-value, which indicates that age does have an effect on age and memory/learning. 

#### I will use tapply and graphical displays to determine direction of the effect of age on memory and learning

> tapply(COUNT, list(AGE, TREATMENT), mean) ### what is the general impact of age on whether mice end up in the control or reward


### These results indicate that older mice, especially 10 week old mice, are more likely to go to the reward section during the testing phase. We next want to see if this is due to learning/memory or both

> ftable(tapply(COUNT, list(NUM_TRAIN,AGE, TREATMENT), mean)) ### what is the direction of the effect of age on learning?

###output#### First column is number of trainings, second column is age, third is mice in control section, fourth is mice in reward section


##### in mice with 1 or 2 trainings, age appears to have no impact on mice learning (mice are no more likely to navigate to the reward vs. the control region). However, the impact of age is pronounced in mice that have 3, 4, or 5 trainings. In each training number, all ages groups mice are more likely to go to the reward section than the control, yet age 10 mice are the most likely to go to the reward for each training section. The results indicate that mice learn better with increased age and this effect peaks at age 10, with a slight decrease in learning ability from 10 to 15 weeks. 

### use ggplot2 to show this information graphically

> aggregate(COUNT ~ AGE + TREATMENT + NUM_TRAIN, FUN = mean) -> num_train_agg ### use aggregate to make a table that has the same information as tapply

library(ggplot2)

> p <- ggplot(num_train_agg, aes(AGE, COUNT))

> p + geom_point(aes(colour = factor(TREATMENT))) + geom_point(aes(shape = factor(NUM_TRAIN))) + scale_color_brewer(palette = "Set1") + scale_shape(solid=F) + ggtitle("Age impacts learning ability") + scale_x_continuous(breaks = c(1,5,10,15)) ### create the plot

##### In the plot we see that age 10 has higher mean counts of mice than all other ages in the reward section after 3, 4, and 5 trainings. Across all age groups, mice with 1 or 2 trainings are no more likely to remember the reward than the control. Additionally, we can see that age 15 mice with 3, 4, or 5 trainings are more likely to find the reward section, but the effect is not as large as what we see with age 10 mice. Overall, this plot reinforces the data from tapply. Age 10 mice learn the best, while age 1 mice show little difference in finding the reward section when they have 1, 2, 3, 4 or 5 trainings. 

###### investigate  the impact of age on memory 
> ftable(tapply(COUNT, list(RUN, AGE, TREATMENT), mean)) #### Run is a proxy for memory because we want to see how many runs are needed before the mice forget the location of the reward. 

### output#### First column is number of runs, second column is age, third is mice in control section, fourth is mice in reward section
con

#### This table show the impact of age on memory. In the first two runs, there is really no impact of age on memory; all aged mice find the location of the reward in a similar proportion. At run 4, age 1 mice are no more likely to remember the reward than the control, while ages 5, 10, and 15 are more likely to remember the reward, with age 10 being the most likely to remember. At run 5, both ages 10 and 15 are more likely to remember the reward but age 10 is only the age to remember the reward at a higher rate than the control in runs 6 and 7. This finding is similar to the learning example. Age 1 mice are the least likely to learn or remember, followed by age 5, then age 15, while age 10 appears to be the age at which the mouse has the best neural function. 

> aggregate(COUNT ~ AGE + TREATMENT + RUN, FUN = mean) -> run_agg
> q <- ggplot(run_agg, aes(AGE, COUNT))
>  q + geom_point(aes(colour = factor(TREATMENT))) + geom_point(aes(shape = factor(RUN))) + scale_color_brewer(palette = "Set1") + scale_shape(solid=F) + scale_shape_manual(values = c(0:15)) + scale_x_continuous(breaks = c(1,5,10,15)) + ggtitle("Age impacts memory")

##### This plot shows the same trends as the tapply table. We see that at all ages mice are more likely to go to the reward section in the first two runs. We see this trend continue until the 7th run for age 10 mice but only to run three for age one and five mice, and run 5 for age 15 mice. This indicates that it takes longer for age 10 mice to forget the location of the reward, when compared to other aged mice. Aged 15 mice perform better than ages 1 and 5 but worse than age 10 mice, therefore we see a peak of memory performance at age 10.  

> ftable(tapply(COUNT, list(RUN, NUM_TRAIN, AGE, TREATMENT), mean))

###output### First column is run, second column is number of trainings, 3rd column is age, fourth is control and fifth is reward. 

###This table is very long but it tells the same story as the previous two tests. Mice are most likely to find the reward section in the first run when they have 5 trainings, and there is no difference in age in the first two runs. However, in the 3rd-8th run we see that age 10 mice perform better than their 1, 5, and 15 week old mice when they have 3 or more trainings.  
