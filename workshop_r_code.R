##############################################################################################
##This code represents instructional material, created by Octavia Dancu for                 ##
##the MiCM workshop. Please do not distribute this without express instructor permission.   ##
##############################################################################################
#Download your data!
uvm_counts<-read.table("/Users/octaviadancu/Downloads/UVM_count_data.txt", stringsAsFactors=FALSE, header=TRUE, sep="\t")
uvm_clin<-read.table("/Users/octaviadancu/Downloads/uvm_clinical_data.txt", stringsAsFactors=FALSE, header=TRUE, sep="\t")
#for windows, it will look something like this!
uvm_counts<-read.table("C:\\Users\\ariana\\Downloads\\UVM_count_data.txt", stringsAsFactors=FALSE, header=TRUE, sep="\t")
uvm_clin<-read.table("C:\\Users\\ariana\\Downloads\\uvm_clinical_data.txt", stringsAsFactors=FALSE, header=TRUE, sep="\t")
#uvm_clin<-uvm_clin[-1,]
#let's take a quick first look at our datasets to get a sense of what we're working with!
dim(uvm_counts)
dim(uvm_clin)
str(uvm_counts)
str(uvm_clin)

#here's a quick example on how to subset your dataframe
#to get the first 10 rows and first two columns of the clinical data:

a<-uvm_clin[1:10,1:2]
uvm_clin[1:4,1:3]

a<-uvm_clin$gender
a
a<-uvm_clin$bcr_patient_barcode
install.packages("ggplot2")
library("ggplot2")

#let's have a quick example of how to generate a scatter plot in R with ggplot2
ggplot(uvm_counts, aes(x=uvm_counts$DNMT3A, y=uvm_counts$DNMT3B)) + geom_point()

#here's a barplot representing eye color across the uvm cohort
ggplot(uvm_clin, aes(uvm_clin$eye_color))+geom_bar()
#histogram
ggplot(uvm_clin, aes(uvm_clin$age_at_initial_pathologic_diagnosis)) + geom_histogram()
#boxplot
ggplot(uvm_clin, aes(x=uvm_clin$eye_color, y=uvm_clin$age_at_initial_pathologic_diagnosis))+geom_boxplot()
#violin plot
ggplot(uvm_clin, aes(x=uvm_clin$eye_color, y=uvm_clin$age_at_initial_pathologic_diagnosis))+geom_violin()
#density plot
ggplot(uvm_clin, aes(uvm_clin$age_at_initial_pathologic_diagnosis)) + 
  geom_density()

####################################################################

install.packages("ggpubr")
library(ggpubr)
#global p-value with anova:
p<-ggplot(iris, aes(x=Species, y=Sepal.Length))+geom_boxplot()
p
p+ stat_compare_means(method = "anova")

a<-c("setosa", "virginica")
a
#what if we want to make pairwise comparisons?
my_comparisons <- list( c("setosa", "virginica"), c("setosa", "versicolor"), c("virginica", "versicolor") )
my_comparisons
#let's add the global pvalue on the plot as well!
p+stat_compare_means(comparisons = my_comparisons)
p+stat_compare_means(comparisons = my_comparisons)+ stat_compare_means(method = "anova")
#is there anything off with this?
#what do we want to change with this to make the labels clearer?
p+stat_compare_means(comparisons = my_comparisons)+ stat_compare_means(method = "anova", label.y = 9.5)
#let's store the plot we made as an object named fig
fig<-p+stat_compare_means(comparisons = my_comparisons)+ stat_compare_means(method = "anova", label.y = 9.5)

#what if we wanted to say that setosa is the reference, and we want to see if there is a significant
#difference between the sepal length of the other two groups vs the reference?
p+stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "setosa")

#I want my significance asterisks to be in red!
p+stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "setosa", color="red")
#https://www.r-bloggers.com/2017/06/add-p-values-and-significance-levels-to-ggplots/

####################################################################

#error bars and standard deviation represented:
#mean standard deviation
ggerrorplot(data = iris,
            desc_stat = "mean_sd",
            y = "Sepal.Length",
            x = "Species")

#but what if we want to see the underlying data?
ggerrorplot(data = iris,
            desc_stat = "mean_sd",
            y = "Sepal.Length",
            x = "Species",
            add = "jitter")

#here's how to add error bars!
ggerrorplot(data = iris,
            y = "Sepal.Length",
            x = "Species",
            add = "jitter", error.plot = "errorbar")
#very nice! but what about some color to differentiate the groups visually a bit better?
ggerrorplot(data = iris,
            y = "Sepal.Length",
            x = "Species",
            color="Species",
            add = "jitter", error.plot = "errorbar")

#to change the color scheme used:

ggerrorplot(data = iris,
            y = "Sepal.Length",
            x = "Species",
            color="Species",
            add = "jitter", error.plot = "errorbar", palette=c("black", "pink", "purple"))

#now let's see how to paste together plots in R to make a multi-panel figure!
#let's first save our new figure as fig2
fig2<-ggerrorplot(data = iris,
               y = "Sepal.Length",
               x = "Species",
               color="Species",
               add = "jitter", error.plot = "errorbar")



####################################################################

ggarrange(fig,fig2)
#what if we want them one on top of the other?
ggarrange(fig,fig2, ncol=1, nrow=2)
#let's add some labels!
ggarrange(fig,fig2, labels = c("A", "B"))

#https://www.r-bloggers.com/2017/07/ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/

####################################################################

#let's save this plot in pdf
pdf(file="eye_color.pdf")
ggplot(uvm_clin, aes(x=uvm_clin$eye_color, y=uvm_clin$age_at_initial_pathologic_diagnosis))+geom_boxplot()
dev.off()
####################################################################

ggplot(uvm_counts, aes(x=uvm_counts$NSD1, y=uvm_counts$TSPAN6))+geom_point()

#let's make a dendrogram!

uvm_counts_t<-t(uvm_counts)
#keep only rows of the expression data that represent genes that are expressed 
#in at least one sample (mean expression > 0)
sub_UVM <- uvm_counts_t[rowMeans(uvm_counts_t) > 0,]
#get the variance for each gene
vars <- apply(sub_UVM,1,var)
#get the 100 genes with the most variance
select <- order(vars,decreasing = T)[1:100]
sub_UVM <- sub_UVM[select,]
#to cluster:
clust_UVM<-hclust(dist(t(sub_UVM)))
plot(as.dendrogram(clust_UVM))
####################################################################

#cladogram

install.packages("ape") 
library("ape")

plot(as.phylo(clust_UVM), type = "cladogram" , cex = 0.6)
dev.off()
plot(as.phylo(clust_UVM), type = "fan")
par(mar = c(2, 2, 2, 2)) 
plot(as.phylo(clust_UVM), type = "fan")
####################################################################
####################################################################
#fit <- survfit(Surv(time, status) ~ sex, data = lung)

require("survival")
install.packages("survminer")
library(survminer)

#let's take a look at the myeloid dataset!
head(myeloid)
#let's see if sex made a difference in terms of outcome for these patients

#now fit survival info to a curve, stratifying info based on the patients' sex
fit <- survfit(Surv(futime, death) ~ sex, data = myeloid)
#let's plot the survival curves
ggsurvplot(fit, data = myeloid, pval = TRUE)



#now, let's see if the different treatment type made a difference!


fit <- survfit(Surv(futime, death) ~ trt, data = myeloid)
ggsurvplot(fit, data = myeloid, pval = TRUE)

####################################################################

uvm_clin$gender<-factor(uvm_clin$gender, levels = c("FEMALE", "MALE"), labels = c("F", "M"))
uvm_clin$NSD1_status<-"LOW"
for(i in 1:nrow(uvm_counts)){
  if(uvm_counts$NSD1[i]>=3000){
    uvm_clin$NSD1_status[i]<-"MEDIUM"
  }
  if(uvm_counts$NSD1[i]>=3158.5){
    uvm_clin$NSD1_status[i]<-"HIGH"
  }
}
uvm_clin$NSD1_status<-factor(uvm_clin$NSD1_status, levels = c("HIGH", "MEDIUM", "LOW"), labels = c("H", "M", "L"))

#if there was an futime column, it would replace it with
#days_to_death. In this case, given that futime doesn't yet 
#exist in uvm_clin, it is being generated and added as the
#last column of uvm_clin
uvm_clin$futime<-uvm_clin$days_to_death
#making a fustat column filled with 1s. 1 means that the
#patient is not censored
uvm_clin$fustat<-c(1)
for(i in 1:nrow(uvm_clin)){
  if(uvm_clin$futime[i]=="[Not Applicable]"){
    #this patient is still alive last time we checked
    uvm_clin$futime[i]<-uvm_clin$days_to_last_followup[i]
    #censor this patient after the last follow-up, we don't know when they will die/if they died
    #after the last check up
    uvm_clin$fustat[i]<-0
  }
}
surv_object <- Surv(time = as.numeric(uvm_clin$futime), event = uvm_clin$fustat)


fit1<- survfit(surv_object ~ gender, data = uvm_clin)
fit2<- survfit(surv_object ~ NSD1_status, data = uvm_clin)
#now let's use ggsurvplot function to plot the survival curves!
ggsurvplot(fit1, data = uvm_clin, pval = TRUE)
ggsurvplot(fit2, data = uvm_clin, pval = TRUE)

####################################################################

#skipping this exercise, replacing with Manhattan plots

# from the chord diagram demo here: https://jokergoo.github.io/circlize_book/book/the-chorddiagram-function.html
install.packages("circlize")
library(circlize)
set.seed(999)
mat = matrix(sample(18, 18), 3, 6) 
rownames(mat) = paste0("S", 1:3)
colnames(mat) = paste0("E", 1:6)
mat
#let's see the directionality of the data a bit more clearly

df = data.frame(from = rep(rownames(mat), times = ncol(mat)),
                to = rep(colnames(mat), each = nrow(mat)),
                value = as.vector(mat),
                stringsAsFactors = FALSE)
df

#you can plot the chord diagram using whichever dataformat is
#easier for you!
chordDiagram(mat)
chordDiagram(df)
#########################################################

#manhattan plots:
install.packages("qqman")
library(qqman)
head(gwasResults)
manhattan(gwasResults, chr="CHR", bp="BP", snp="SNP", p="P" )
manhattan(gwasResults, annotatePval = 0.01)


####################################################################
####################################################################

#considering color:
#here is an example of uniform color scheme:
ggplot(uvm_clin, aes(uvm_clin$eye_color))+geom_bar(fill="red")
ggplot(uvm_clin, aes(uvm_clin$eye_color))+geom_bar(color="red")
#what about for scatter plots? What should we do to get a scatterplot with red dots?
#fill or color?
ggplot(uvm_counts, aes(x=uvm_counts$DNMT3A, y=uvm_counts$DNMT3B)) + geom_point(fill="red")
ggplot(uvm_counts, aes(x=uvm_counts$DNMT3A, y=uvm_counts$DNMT3B)) + geom_point(color="red")

#let's revisit the eye colour plot and try to make it a bit more user friendly using the color scheme!
#here is the plain one, that doesn't have any changes made to the color scheme
ggplot(uvm_clin, aes(uvm_clin$eye_color))+geom_bar()
#is it easy to see at a glance what each of these means and compare?
#let's change the colors to reflect the eye color of the bar!
ggplot(uvm_clin, aes(uvm_clin$eye_color))+geom_bar(fill=c("gray", "black", "blue", "brown", "green"))
#okay, great! But the colours are not quite what I want them to be (brown looks a bit too red to me!)
#let's use the table of ggplot colours to improve this!
ggplot(uvm_clin, aes(uvm_clin$eye_color))+geom_bar(fill=c("gray", "black", "cornflowerblue", "sienna4", "forestgreen"))
#much better!
#Bonus: how do we only represent the known eye colours (omitting the NA and Unk)?
which(uvm_clin$eye_color=="[Not Available]")
which(uvm_clin$eye_color=="[Unknown]")
a<-c(which(uvm_clin$eye_color=="[Not Available]"), which(uvm_clin$eye_color=="[Unknown]"))
ggplot(uvm_clin[-a,], aes(uvm_clin[-a,]$eye_color))+geom_bar(fill=c("cornflowerblue", "sienna4", "forestgreen"))

#
#ggplot(uvm_clin, aes(as.numeric(uvm_clin$age_at_initial_pathologic_diagnosis)))+geom_histogram()


#boxplot example
ggplot(uvm_clin, aes(x=uvm_clin$eye_color, y=as.numeric(uvm_clin$age_at_initial_pathologic_diagnosis)))+geom_boxplot()
ggplot(uvm_clin, aes(x=uvm_clin$eye_color, y=as.numeric(uvm_clin$age_at_initial_pathologic_diagnosis), fill=uvm_clin$eye_color))+geom_boxplot()+scale_fill_brewer()
#ggplot(uvm_clin, aes(x=uvm_clin$eye_color, y=as.numeric(uvm_clin$age_at_initial_pathologic_diagnosis)))+geom_boxplot()


#density plot

ggplot(uvm_clin, aes(as.numeric(uvm_clin$age_at_initial_pathologic_diagnosis)))+geom_density()


#bar plot, changing fill
#ggplot(uvm_clin, aes(uvm_clin$eye_color))+geom_bar(fill="red")

#example of brewer use
ggplot(uvm_clin, aes(x=uvm_clin$eye_color, y=as.numeric(uvm_clin$age_at_initial_pathologic_diagnosis), fill=uvm_clin$eye_color))+geom_boxplot()+scale_fill_brewer("Blues")
ggplot(iris, aes(iris$Sepal.Length, iris$Sepal.Width, color=iris$Species))+geom_point()+scale_color_brewer(palette="Paired")
#revisiting the pathologic stage
ggplot(uvm_clin, aes(pathologic_stage, fill=uvm_clin$pathologic_stage)) + geom_bar()+scale_fill_brewer(palette="Reds")
#ggplot(uvm_clin, aes(age_at_initial_pathologic_diagnosis, fill=uvm_clin$gender)) + geom_density(alpha=0.4)

#viridis
install.packages("viridis")
library(viridis)
ggplot(iris, aes(iris$Sepal.Length, iris$Sepal.Width, color=iris$Species))+geom_point()+scale_color_viridis(discrete=TRUE)

ggplot(iris, aes(iris$Sepal.Length, iris$Sepal.Width, color=iris$Sepal.Length))+geom_point()+scale_color_viridis()

#wesanderson 
install.packages("devtools")
library(devtools)
devtools::install_github("karthik/wesanderson")
library(wesanderson)
ggplot(iris, aes(iris$Sepal.Length, iris$Sepal.Width, color=iris$Species))+geom_point()+scale_color_manual(values = wes_palette("GrandBudapest1", n = 3, "discrete"))

#showing how metbrewer works
install.packages("MetBrewer")
library(MetBrewer)
ggplot(iris, aes(x=iris$Species, y=iris$Sepal.Length, fill=iris$Species))+geom_boxplot()+ scale_fill_manual(values=met.brewer("Ingres", 3))
#let's change the direction of the colors from the palette!
ggplot(iris, aes(x=iris$Species, y=iris$Sepal.Length, fill=iris$Species))+geom_boxplot()+ scale_fill_manual(values=met.brewer("Ingres", 3, direction=-1))



###################################################################################
devtools::install_github("AndreaCirilloAC/paletter")
library(paletter)
#here, you need to fill in the path to the photo you want to get the colors from!
image_path <- "/Users/octaviadancu/Downloads/marie_antoinette.jpeg"
colours_vector <- create_palette(image_path = image_path,
                                 number_of_colors =3,
                                 type_of_variable = "categorical")

ggplot(iris, aes(x=Species, y=Sepal.Length))+geom_boxplot(fill=colours_vector)

###################################################################################
#Exercise 3 Answer:
ggplot(uvm_clin, aes(age_at_initial_pathologic_diagnosis, fill=uvm_clin$gender)) + geom_density()
ggplot(uvm_clin, aes(age_at_initial_pathologic_diagnosis, fill=uvm_clin$gender)) + geom_density(alpha=0.4)

ggplot(uvm_clin, aes(age_at_initial_pathologic_diagnosis)) + geom_density()
###################################################################################
#ghibli color schemes

install.packages('ghibli')
library(ghibli)

fig <- ggplot(iris, aes(Sepal.Length, fill = Species)) + geom_density(alpha = 0.7)+scale_fill_ghibli_d("PonyoMedium")
fig

###################################################################################
#colorblindr

#install the different packages we'll need
remotes::install_github("wilkelab/cowplot")
library(cowplot)
#install.packages("colorspace", repos = "http://R-Forge.R-project.org")
#library(colorspace)
remotes::install_github("clauswilke/colorblindr")
library(colorblindr)
library(ggplot2)
#let's make a ggplot!
fig <- ggerrorplot(data = iris,
                   y = "Sepal.Length",
                   x = "Species",
                   color="Species",
                   add = "jitter", error.plot = "errorbar")
fig

#let's now see how it looks for someone colorblind - grid with multiple different types of colorblindness
cvd_grid(fig)
#interactive display showing plot as seen by someone with colorblindness - can change up the severity
view_cvd(fig) 

# from https://github.com/clauswilke/colorblindr

#pretty good, but can we make it better?
p<-ggerrorplot(data = iris,
               y = "Sepal.Length",
               x = "Species",
               color="Species",
               palette="Paired",
               add = "jitter", error.plot = "errorbar")
cvd_grid(p)
#much better!

###################################################################################


#sizes
ggplot(iris, aes(iris$Sepal.Length, iris$Sepal.Width, color=iris$Sepal.Length))+geom_point(size=2)+scale_color_viridis()
ggplot(iris, aes(iris$Sepal.Length, iris$Sepal.Width, color=iris$Sepal.Length))+geom_point(size=4)+scale_color_viridis()


#boxplot, editing titles and axes
ggplot(uvm_clin, aes(x=uvm_clin$eye_color, y=as.numeric(uvm_clin$age_at_initial_pathologic_diagnosis)))+geom_boxplot()
ggplot(uvm_clin, aes(x=uvm_clin$eye_color, y=as.numeric(uvm_clin$age_at_initial_pathologic_diagnosis)))+geom_boxplot(fill=c("gray", "black", "cornflowerblue", "sienna4", "forestgreen"))
#adding title
ggplot(uvm_clin, aes(x=uvm_clin$eye_color, y=as.numeric(uvm_clin$age_at_initial_pathologic_diagnosis)))+geom_boxplot(fill=c("gray", "black", "cornflowerblue", "sienna4", "forestgreen"))+ggtitle("Age at Initial Diagnosis, stratified by Eye Colour")

###################################################################################

#cleaning up the example plot
ggplot(iris, aes(iris$Sepal.Length, iris$Sepal.Width, color=iris$Species))+geom_point()+scale_color_brewer(palette="Paired")

#background
ggplot(iris, aes(iris$Sepal.Length, iris$Sepal.Width, color=iris$Species))+geom_point()+scale_color_brewer(palette="Paired")+theme_minimal()

ggplot(iris, aes(iris$Sepal.Length, iris$Sepal.Width, color=iris$Species))+geom_point()+scale_color_brewer(palette="Paired")+theme_light()

#getting rid of the legend
ggplot(iris, aes(iris$Sepal.Length, iris$Sepal.Width, color=iris$Species))+geom_point()+scale_color_brewer(palette="Paired")+theme_light()+theme(legend.position = "none")

ggplot(iris, aes(iris$Sepal.Length, iris$Sepal.Width, color=iris$Species))+geom_point()+scale_color_brewer(palette="Paired")+theme_light()+theme(legend.position = "none")+ggtitle("Iris")+theme(plot.title=element_text(hjust = 0.5))

###################################################################################
#upgrading the boxplot
#making up some dataset
a<-c(1,2,1,3,2,1,1,2)
b<-c(10,13,12,10,14,13,15)
ab<-c(rep("A", 8), rep("B", 7))
abc<-c(a,b)
d<-cbind.data.frame(abc, ab)
ggplot(d, aes(x=ab,y=d$abc))+geom_boxplot()
ggplot(d, aes(x=ab,y=d$abc))+geom_violin()
ggplot(d, aes(x=ab,y=d$abc))+geom_boxplot()+geom_jitter(alpha=0.5)
ggplot(d, aes(x=ab,y=d$abc))+geom_jitter(position=position_jitter(0.1))
ggplot(d, aes(x=ab,y=d$abc))+geom_violin()+geom_jitter(position=position_jitter(0.1))

###################################################################################
###################################################################################

#facet grid demo:

a<-ggplot(uvm_clin, aes(x=uvm_clin$new_tumor_event_after_initial_treatment))+geom_bar()
a+facet_grid(cols = vars(eye_color))
#let's add some colors to make this a bit easier!
a<-ggplot(uvm_clin, aes(x=uvm_clin$new_tumor_event_after_initial_treatment, fill=uvm_clin$eye_color))+geom_bar()
a+facet_grid(cols = vars(eye_color))+ scale_fill_manual(values=c("Black", "grey", "cornflowerblue", "sienna4", "forestgreen"))

###################################################################################

#linear regression
a<-ggplot(uvm_counts, aes(x=uvm_counts$DNMT3A, y=uvm_counts$DNMT3B)) + geom_point()+ stat_smooth(method = lm)
a
model <- lm(uvm_counts$DNMT3B ~ uvm_counts$DNMT3A, data = uvm_counts)
summary(model)
#correlation
cor(uvm_counts$DNMT3B, uvm_counts$DNMT3A)

###################################################################################

#let's make a heatmap!

sub_UVM <- uvm_counts[rowMeans(uvm_counts) > 0,]
vars <- apply(sub_UVM,1,var)
select <- order(vars,decreasing = T)[1:80]
sub_UVM <- sub_UVM[select,1:100]

clust_UVM<-hclust(dist(sub_UVM))
install.packages('pheatmap')
library("pheatmap")
pheatmap(sub_UVM)
###################################################################################
install.packages("heatmaply")
library(heatmaply)
heatmaply(sub_UVM)
###################################################################################
install.packages("plotly")
library(plotly)
a<-ggplot(uvm_counts, aes(x=uvm_counts$DNMT3A, y=uvm_counts$DNMT3B, color=rownames(uvm_counts))) + geom_point()+theme(legend.position = "none") + stat_smooth(method = lm)
ggplotly(a)

devtools::install_github("octaviamd/interactmapper")
library(interactmapper)

#example of use of interactmapper function, interact_multi
interact_multi(iris[,1:4], iris$Species, iris[,1:2], "UMAP", "viridis", "Species", c("Sepal Length", "Sepal Width"))
interact_multi(uvm_counts, as.factor(uvm_clin$new_tumor_event_after_initial_treatment), uvm_clin$eye_color, "UMAP", "viridis", "New Tumor Event", "Eye Color")
