setwd("C:/Users/anave/Desktop/Bioinfromatics/Semester 2/DNA RNA Dynamics/Report")

library(minfi)
list.files("C:/Users/anave/Desktop/Bioinfromatics/Semester 2/DNA RNA Dynamics/Report/Input_data_report")

###STEP 1: Load raw data with minfi and create an object called RGset storing a RGChannelSet object 

#Loading the sample sheet into a data frame
SampleSheet <- read.csv("C:/Users/anave/Desktop/Bioinfromatics/Semester 2/DNA RNA Dynamics/Report/Input_data_report/Samplesheet_report_2022.csv",header=T, stringsAsFactors = T)
SampleSheet

#Setting a directory from which my raw data will be imported
baseDir <- ("C:/Users/anave/Desktop/Bioinfromatics/Semester 2/DNA RNA Dynamics/Report/Input_data_report")

#Loading the sample sheet into a data frame.
#Using the read.metharray.sheet() function. 
?read.metharray.sheet
#It takes as input a directory where it starts searching for an object containing pheno-data information for the samples in an Illumina methylation experiment.
targets <- read.metharray.sheet(baseDir)
targets

#Reading the methylation array experiment raw data from a data frame containing information about our samples and corresponding arrays.
#Using the read.metharray.exp() function.
?read.metharray.exp
#Providing the data frame containing the data about the experiment.
#The function reads it into an object of class RGChannelSet. 
#It stores data for both red and green channel. Data is be organized on the probe level, i.e. according to addresses of probes, not according to IDs of CpG sites. From this object, by applying appropraite accessor functions, we can obtain information about fluorescence emitted by each probe.
RGset <- read.metharray.exp(targets = targets)
#Saving it
save(RGset,file="RGset.RData")

#Exploring my object
RGset
str(RGset)

###STEP 2: Create data frames Red and Green to store the red and green fluorescence respectively 

#Using the getRed function of the minfi library, which takes as input an RGChannelSet object, and returns the information from the red channel of the object. Saving the output as a data frame
Red <- data.frame(getRed(RGset))
#Checking what the data frame looks like
head(Red)
#Columns correspond to arrays, rows represent addresses of probes for which signal intensity is given

#Repeat the same for the green channel with the getGreen function
Green <- data.frame(getGreen(RGset))
head(Green)

###STEP 3: Get fluorescence for address 46801437

#Checking if my address corresponds to a type I or type II probe. I do this by looking at the Illumina 450k Manifest
load("C:/Users/anave/Desktop/Bioinfromatics/Semester 2/DNA RNA Dynamics/R workingdirectory/Illumina450Manifest_clean.RData")
#First checking if it has an address A, if it does, it can be either Type I or type II, if it doesn't correspond to and Address A, then it's surely a type I since only those probes have address B
myAddress <-Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID=="46801437",c("Infinium_Design_Type","Color_Channel")]
myAddress
#I see it's type II

#Getting the red fluorescence of my address by sub-setting the data frame containing fluorescence data from the red channel
addressRed <- Red[rownames(Red)=="46801437",]
addressRed

#Getting the green fluorescence of my address by sub-setting the data frame containing fluorescence data from the green channel
addressGreen <- Green[rownames(Green)=="46801437",]
addressGreen

###STEP 4: Create the object MSet.raw 

#A MethylSet object is obtained by applying the preprocessRaw function on the RGChannelSet.
?preprocessRaw
#MethylSet object stores prepossessed data about the methylation signal (data without normalization) organized according to CpG sites. 
#Information about which CpG was found methylated and which unmethylated can be derived from it by applying appropriate accessor functions.
MSet.raw <- preprocessRaw(RGset)
MSet.raw
save(MSet.raw,file="MSet_raw.RData")

#In next steps I'm extracting data about CpGs found methylated/unmethylated from the MethylSet object to see how it is organized.
?getMeth
#getMeth() is an accessor function of a MethylSet object. It is used to the obtain intensity of fluorescence emitted by methylated signals.
Meth <- getMeth(MSet.raw)
#In the output object rows represent indivudual CpG sites, and columns represent samples.
head(Meth)

#getUnmeth() is homologous to getMeth(), but it is used to obtain the intensity of fluorescence emitted by unmethylated signals. 
Unmethylated <- getUnmeth(MSet.raw)
#The output is organized in the same way as output of getMeth()
head(Unmethylated)



###STEP 5: Quality control - QC plot, negative control probes, detection p-values

#To get a QC plot I need the median intensities of mehtylated and unmethylated CpGs of each sample. 
#I can obtain that from MSet.raw by using the function getQC(). I store them in the QC data frame.
?getQC
QC <- getQC(MSet.raw)
QC
#plotQC() is a minfi function that outputs the actual QC plot. 
#I expect all my samples to be above the threshold (dotted line). If some samples are located bellow the threshold, this indicated their quality isn't good. 
#pdf() will save the plots I generate as a pdf file in my current working directory
pdf("QC-plot.pdf",width=7,height=7)
plotQC(QC)
dev.off()
#In the plot I see all my samples are of good quality. 

#Checking the intensities of negative probes is done to make sure the background is not too strong. 
#All negative probes are expected to have signal intensity between 100 and 1000 units.
#I will check this by making a Strip Plot with the controlStripPlot function.
#First I look at what this function needs as arguments.
?controlStripPlot
#First argument is an RGChannel object, and I have to extract the exact name of controls I want to plot from that object.
#getProbeInfo() function takes as argument a RGChannel object and returns information on the probes specified in argument 'type' 
allControlProbes <- data.frame(getProbeInfo(RGset, type = "Control"))
#I check the what are the columns of the dataframe I created
str(allControlProbes)
#I see that type of control probe is saved in column 'Types'
#I want to see all categorical variables in column 'Types', so I can get the exact label of negative probes
table(allControlProbes$Type)
#Now I have all information needed to use controlStripPlot()
pdf("Negative-controls.pdf",height=7,width=7)
controlStripPlot(RGset, controls="NEGATIVE")
dev.off()
#I see all for all samples negative probes in both color channels have log2 intensity below 10, which indicated that non of the samples have problematic background

#Detection p-values give the probability that the fluorescence you are measuring is different from the fluorescence of the background.
#All probes with a detection p-value higher than the threshold cannot be considered as giving reliable information.
#Samples shouldn't have a large portion of probes with detection p-values above the thresholds.If they do, we will discard such samples. 
#Now I'm checking for each sample the number of probes with an unsatisfactory detection p-value.
#Function detectionP() identifies can help us identify failed positions from an RGChannel object.In its output CpG site IDs are rows, samples are columns, and each entry is a detection p-value of for a CpG site in the corresponding sample.
detectP <- detectionP(RGset)
save(detectP,file="detP.RData")
load("detP.RData")
#Inspecting how the object looks
head(detectP)

#Select the entries where detection p-values are above the threshold of 0,01. 
#In the output object TRUE represents the failed positions, and FALSE represents not failed positions.
failed <- detectP>0.01
head(failed)
#with the summary function we can check how many positions have failed and how many haven't in each sample. For each sample it will return a count of both cases.
summary(failed)

#Calculating if some samples have >5% of failed positions.
#colMeans function will be used on the object containing Boolean values illustrating which position failed and which is good.
#It will calculate means of all columns (TRUE/(TRUE+FALSE)), and since each column represents a sample, we can check if each sample has a proportion of failes positions >0.05
failed_percentage <- colMeans(failed) 
head(failed_percentage)
#Checking the output of colMeans shows all samples are of sufficient quality, so all will be retained for further analysis.

#Now we will create an object containing all probes with a detection p-value >0.01. 
#They will be removed from the data set after normalization, since the normalization procedure is disturbed by missing data.
#rowMeans function will calculate the proportion of failed positions for each probe. 
#In the output object each value will represent TRUE/(TRUE+FALSE) for each probe across all 8 samples.
means_of_rows <- rowMeans(failed)
head(means_of_rows)

#Creating an object containing a TRUE or FALSE value for each probe, depending on if it's proportion of failed positions in all samples is <0.01
retain_probes <- means_of_rows<0.01
head(retain_probes)
#Checking how many probes have failed
table(retain_probes)
#730 probes can be filtered out after normalization, since the functions used in the normalization procedure don't look kindly on missing data. 
#We save the CpG IDs corresponding to these probes in a vector. And later on these probes can be removed from the normalized data set.
#Achieving this by sub-setting the object containing FALSE/TRUE info for each position
remove_probes <- retain_probes[retain_probes==FALSE]
#Checking if the vector is of the same length as is the number of FALSE from table(retain_probes)
length(remove_probes)



###STEP 6: Get raw beta and M values. Plot the densities of mean methylation values for WT and MUT samples separately
#Using the getBeta(MethylSet object) function to get a matrix of beta values. Beta values will range from 0 to 1
beta <- getBeta(MSet.raw)
head(beta)

#Using the getM(MethylSet object) function to get M values.M values will range from -inf to +inf
M <- getM(MSet.raw)
head(M)

#I checked in the Sample Sheet which samples are MUT and which are WT. 
#Sub-setting the beta and M matrices to obtain matrices containing M and beta values only of MUT and only of WT samples
MUT_M <- M[,c(2,5,6,8)]
head(MUT_M)
MUT_beta <- beta[,c(2,5,6,8)]
head(MUT_beta)
WT_M <- M[,c(1,3,4,7)]
head(WT_M)
WT_beta <- beta[,c(1,3,4,7)]
head(WT_beta)

#Calculating means for beta and M values for both MUT and WT samples.
#Using the apply(matrix,MARGIN,mean,na.rm) function to calculate mean of all rows in the input matrix. 
#MARGIN parameter will be set to 1 to indicate that means should be calculated over rows. Also, na.rm parameter will be set to true so it removes missing data when it is encountered. 
#Removing missing values is needed because missing values will influence the density function we will use later
#Output of this function is a vector containing mean values of each row in input matrix.

#This returns a vector containing mean beta values for each position in WT samples.
beta_mean_WT <- apply(WT_beta,1,mean,na.rm=T)
head(beta_mean_WT)
#Removing NA values if they remained in the output vector
beta_mean_WT <- beta_mean_WT[!is.na(beta_mean_WT)]

#This returns a vector containing mean M values for each position in WT samples.
M_mean_WT <- apply(WT_M,1,mean,na.rm=T)
head(M_mean_WT)
#Removing NA values if they remained in the output vector
M_mean_WT <- M_mean_WT[!is.na(M_mean_WT)]

#This returns a vector containing mean beta values for each position in MUT samples.
beta_mean_MUT <- apply(MUT_beta,1,mean,na.rm=TRUE)
head(beta_mean_MUT)
#Removing NA values if they remained in the output vector
beta_mean_MUT <- beta_mean_MUT[!is.na(beta_mean_MUT)]

#This returns a vector containing mean M values for each position in MUT samples.
M_mean_MUT <- apply(MUT_M,1,mean,na.rm=T)
head(M_mean_MUT)
#Removing NA values if they remained in the output vector
M_mean_MUT <- M_mean_MUT[!is.na(M_mean_MUT)]

#Calculating the density distribution for each vector which was just created. We use the denisty() function
d_mean_beta_WT <- density(beta_mean_WT)
d_mean_beta_MUT <- density(beta_mean_MUT)
d_mean_M_WT <- density(M_mean_WT)
d_mean_M_MUT <- density(M_mean_MUT)

pdf("Step-6.pdf",width=10,height=5)
#par(mfrow=) arguemnt enables having more then 1 plot in the same image
par(mfrow=c(1,2))
#plot is a generic R function which produces a plot of R objects
plot(d_mean_beta_WT,main="Density of Beta Values",col="orange")
#lines() is a generic function taking coordinates from the input, and joining the points with  lines. It enables us to produce plots containing multiple lines
lines(d_mean_beta_MUT,col="green")
plot(d_mean_M_WT,main="Density of M Values",col="orange")
lines(d_mean_M_MUT,col="green")
dev.off()
#The density plot of beta values shows that MUT and WT samples have similar methylation profiles. In the both extremes of the distribution (around beta = 0 and beta = 0.9) the WT samples show a higher density compared to MUT samples. In middle ranges of beta values (0.2<beta<0.7) MUT samples have a higher density.
#The density plot of M values corroborates what the beta values says in all parts except the peak at lower values. In the range of lower M values MUT samples have higher density. 

##STEP 7: Normalize the data according to SWAN normalization. Produce 6 plots.

#First we sub-sample the Illumina450Manifest to start separating type I and type II probes in our data
dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="I",]
dfI <- droplevels(dfI)
dim(dfI)
dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="II",]
dfII <- droplevels(dfII)
dim(dfII)

#Normalization with the preprocessSWAN function contained in the minfi package
#The preprocessSWAN function takes as input either an object of class RGChannelSet and optionally an object of class MethylSet which is the result of the same RGChannelSet given as argument.
#preprocessSWAN returns a MethylSet which contains normalized data
preprocessSWAN_results <- preprocessSWAN(RGset)
class(preprocessSWAN_results)
preprocessSWAN_results
#Getting the beta values of normalized data with the getBeta accessor function of the MethylSet class object
beta_preprocessSWAN <- getBeta(preprocessSWAN_results)
head(beta_preprocessSWAN)
#Checking if normalized beta values of all samples are okay (checking if they are all in 0-1 range)
summary(beta_preprocessSWAN)

#Now I need to subset the matrix containing the normalized beta values to separate the type I and type II probes
#I subset according to dfI
beta_preprocessSWAN_I <- beta_preprocessSWAN[rownames(beta_preprocessSWAN) %in% dfI$IlmnID,]
#Checking if the number of rows in the type I subset matrix is the same as in dfI
dim(beta_preprocessSWAN_I)
#I subset according to dfII
beta_preprocessSWAN_II <- beta_preprocessSWAN[rownames(beta_preprocessSWAN) %in% dfII$IlmnID,]
#Checking if the number of rows in the type I subset matrix is the same as in dfI
dim(beta_preprocessSWAN_II)

#Calculating the mean normalized beta values of each probe in both type I and type II subsets. 
#Goal is to plot the density of mean normalized beta values of type I and type II probes
#Using again the apply function. Output is a vector
mean_preprocessSWAN_beta_I <- apply(beta_preprocessSWAN_I,1,mean)
mean_preprocessSWAN_beta_II <- apply(beta_preprocessSWAN_II,1,mean)
#Calculating the densities
d_mean_of_preprocessSWAN_beta_I <- density(mean_preprocessSWAN_beta_I,na.rm=T)
d_mean_of_preprocessSWAN_beta_II <- density(mean_preprocessSWAN_beta_II,na.rm=T)

#Calculating the standard deviations of normalized beta values for subsets of type I and type II probes
#Using again the apply function. Output is a vector
sd_preprocessSWAN_beta_I <- apply(beta_preprocessSWAN_I,1,sd,na.rm=T)
sd_preprocessSWAN_beta_II <- apply(beta_preprocessSWAN_II,1,sd,na.rm=T)

#Now getting density distribution of standard deviation for normalized beta values of type I and type II probes
d_sd_preprocessSWAN_beta_I <- density(sd_preprocessSWAN_beta_I,na.rm=T)
d_sd_preprocessSWAN_beta_II <- density(sd_preprocessSWAN_beta_II,na.rm=T)

#From the matrix containing raw beta values of all samples we take just the rows with names that are in dfI. 
#This way we sub-sample the raw beta values of type I probes
beta_I <- beta[rownames(beta) %in% dfI$IlmnID,]
#Checking if the number of rows is the same as in dfI
dim(beta_I)

#Repeating the same procedure to sub-sample raw beta values of type II probes
beta_II <- beta[rownames(beta) %in% dfII$IlmnID,]
dim(beta_II)

#Calculating the mean raw beta values of each probe in both type I and type II subsets. 
#Goal is to plot the density of mean beta values of type I and type II probes
#Using again the apply function. Output is a vector
mean_beta_I <- apply(beta_I,1,mean)
mean_beta_II <- apply(beta_II,1,mean)
#Calculating the densities
d_mean_of_beta_I <- density(mean_beta_I,na.rm=T)
d_mean_of_beta_II <- density(mean_beta_II,na.rm=T)

#Calculating the standard deviations of raw beta values for subsets of type I and type II probes
#Using again the apply function. Output is a vector
sd_beta_I <- apply(beta_I,1,sd,na.rm=T)
sd_beta_II <- apply(beta_II,1,sd,na.rm=T)

#Now getting density distribution of standard deviation for raw beta values of type I and type II probes
d_sd_beta_I <- density(sd_beta_I,na.rm=T)
d_sd_beta_II <- density(sd_beta_II,na.rm=T)

#Now I plot everything
pdf("Step-7.pdf",width=12,height=5)
par(mfrow=c(2,3))
plot(d_mean_of_beta_I,col="blue",main="Density plot of raw beta mean values ")
lines(d_mean_of_beta_II,col="red")
plot(d_sd_beta_I,col="blue",main="Density plot of raw beta SD values")
lines(d_sd_beta_II,col="red")
boxplot(beta)
plot(d_mean_of_preprocessSWAN_beta_I,col="blue",main="Density plot of preprocessSWAN beta mean values")
lines(d_mean_of_preprocessSWAN_beta_II,col="red")
plot(d_sd_preprocessSWAN_beta_I,col="blue",main="Density plot of preprocessSWAN beta SD values")
lines(d_sd_preprocessSWAN_beta_II,col="red")
boxplot(beta_preprocessSWAN)
dev.off()
#Looking at the plots I obtained, I wouldn't say that SWAN is the optimal normalization, since it doesn't seem to have changed the results very much. The densities of mean beta values of the two chemistries are slightly more overlapping after normalization, but only in the beta values in range 0-0.15, after that normalization didn't result in overlapping profiles for type I and type II probes.
#Density of SD of type I and type II probes remained almost the same after normalization .
#Since SWAN is useful in cases when we cannot presume our samples have similar distributions of methylation, I think this is not the case with these samples since the box plots of raw beta values to show a similar distribution in all samples.Perheps the Quantile normalization would perform better.



###STEP 8: identify differentially methylated probes between group WT and group MUT 

#Since I want to identify the differentially methylated sites between two sets of samples, I have to draw information about this from the column 'Group' of the SampleSheet data frame I created in the first step
#Buliding a custom function, so I can use it as input for the apply() function, in order to apply it an all the rows (=all sites) of my matrix with normalized beta values. This function returns p-value for the site it tests
#Inside this function I use the wilcox.tect() function which takes as input a vector of data (which is why it must be run row-by-row) and performs a Wilcoxon test on two samples (which is in fast a Mann-Whitney test). Output of this function is a list and p-values of sites are extracted from it
My_mannwhitney_function <- function(x) {
  wilcox <- wilcox.test(x~ SampleSheet$Group)
  return(wilcox$p.value)
} 

pValues_wilcox <- apply(beta_preprocessSWAN,1, My_mannwhitney_function)
length(pValues_wilcox)
#Merging p-values and beta values into a data frame
final_wilcox <- data.frame(beta_preprocessSWAN, pValues_wilcox)
head(final_wilcox)
#Now sub-setting the data frame I created only for significantly differentially methylated sites, i.e., for sites with p-values <=0.05
significant <- final_wilcox[final_wilcox$pValues_wilcox<=0.05,]
#Checking if I filtered correctly. I expect the p-values I see to be <=0.05
head(significant)
#Checking how many probes are differentially methylated
dim(significant)
#53380 probes are differentially methylated with a p-value <=0.05



###STEP 9: Do multiple testing correction
#Using the function p.adjust() which takes as input a set of p-values and returns p-values adjusted with the correction technique you specify in the method parameter
corrected_Bonferroni <- p.adjust(final_wilcox$pValues_wilcox,"bonferroni")
corrected_BH <- p.adjust(final_wilcox$pValues_wilcox,"BH")
#Creating a data frame with all beta values, unadjusted p-values and both corrected p-values
final_wilcox_corrected <- data.frame(final_wilcox,corrected_BH,corrected_Bonferroni)
head(final_wilcox_corrected)

#Now I will check how many probes are still found significant with a p-value threshold of 0.05 for each type of p-value I calculated
#First checking the unadjusted p-values. The survivning probes are given by the first dimension of the dim() function output
dim(final_wilcox_corrected[final_wilcox_corrected$pValues_wilcox<=0.05,])
#53380 sites are found as differentially methylated when using the unadjusted p-value
#Next I check the p-values adjusted with Benjamin-Hockberg correction
dim(final_wilcox_corrected[final_wilcox_corrected$corrected_BH<=0.05,])
#1 site is found significant after this correction
#Finally, checking the p-values adjusted with the Bonferroni correction
dim(final_wilcox_corrected[final_wilcox_corrected$corrected_Bonferroni<=0.05,])
#Again only 1 probe remains significant after correction



###STEP 10: Producing a volcano plot and a Manhattan plot
library(qqman)
#Volcano plots have on y-axis the -log10(p-value) and on x-axis the difference between the averages of sample groups. So, I need to get this difference
#First I extract the normalized beta values for WT and MUT samples by subseting the beta_preprocessSWAN matrix
beta_norm_MUT <- beta_preprocessSWAN[,SampleSheet$Group=="MUT"]
head(beta_norm_MUT)
beta_norm_WT <-beta_preprocessSWAN[,SampleSheet$Group=="WT"]
head(beta_norm_WT)
#Now I calculate the mean beta values for each probe of both groups.
#I will use the apply() function combined with mean() again
beta_norm_mean_MUT <- apply(beta_norm_MUT,1,mean)
beta_norm_mean_WT <- apply(beta_norm_WT,1,mean)
#Now I can calculae the difference wit a simple - operation
difference_WT_MUT <- beta_norm_mean_WT-beta_norm_mean_MUT
head(difference_WT_MUT)
#Creating a data frame which will be the input for the plot() function. It will contain all values from difference_WT_MUT and -log10(p-value)
toVolcPlot <- data.frame(difference_WT_MUT, -log10(final_wilcox_corrected$pValues_wilcox))
head(toVolcPlot)

#Now I can finally create the plot
pdf("Volcano-plot.pdf",width=10,height=5)
plot(toVolcPlot[,1], toVolcPlot[,2])
dev.off()

#Now I need to create a Manhattan plot.
#Since the Manhattan plot has chromosomes on x-axis I need to annotate the data frame containing the p-values with genomic information about the probes
#I get this information from the IlluminaManifest
#I will use the merge() function. It is a generic function which takes as input two data frames and merges them by common columns or row names.
#First I have to process my data objects, so they are of acceptable format for this
#Since we want to merge based on columns, and probe IDs are stored in row names, this needs to be fixed.
#The data frame I'm creating will be a modification of the data frame containing beta and unadjusted p-values. It will have an additional column which will contain now the name of each row
final_wilcox_manhattan <- data.frame(rownames(final_wilcox),final_wilcox)
head(final_wilcox_manhattan)
#Now I want to rename the column I just added, so it overlaps with the name of the column in the Illumina Manifest where probe IDs are stored.
#I get back the names of columns in Illumina Manifest to see what I can use
colnames(Illumina450Manifest_clean)
#Making shure 'IlmnID' indeed stores probe CpG site IDs
head(Illumina450Manifest_clean)
#Now I rename the column I added in line 402
colnames(final_wilcox_manhattan)[1] <- "IlmnID"
#Checking names of columns in the data frame
colnames(final_wilcox_manhattan)

#Now I can use the merge() function
final_wilcox_annotated <- merge(final_wilcox_manhattan,Illumina450Manifest_clean,by="IlmnID")
#Checking if all rows and columns are present. Number of columns should be the sum of columns in the data frame with p-values and columns in the Illumina manifest - 1 (because one is in common)
dim(Illumina450Manifest_clean)
dim(final_wilcox_manhattan)
dim(final_wilcox_annotated)

#From this I create input object for the function which creates a Manhattan plot. It should contain columns with 4 types of information: probe ID, chromosome where probe is located, p-value, position of the site the probe is targeting on the chromosome
#Checking to see the correct names of columns in the annotated data frame
colnames(final_wilcox_annotated)
#Subsetting
input_Manhattan <- final_wilcox_annotated[colnames(final_wilcox_annotated) %in% c("IlmnID","CHR","MAPINFO","pValues_wilcox")]
head(input_Manhattan)
#Chromosomes should be in order in the Manhattan plot, so I'll check the order of chromosomes in the input data frame
levels(input_Manhattan$CHR)
#I have to reorder chromosomes how I like. I will use a vector as help
order_chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
#Now I use the factor() function which will turn my vector values into factor values and assign them to the CHR column
input_Manhattan$CHR <- factor(input_Manhattan$CHR,levels=order_chr)
#Checking if they are ordered how I wanted
levels(input_Manhattan$CHR)
#Now I can produce a Manhattan plot with the data frame I created using the manhattan() function. 
#This is a function from the qqman library.
?manhattan
#This function takes as input a data frame which contains chromosome, position, and p-value and outputs a Manhattan plot. chr parameter values need to be numeric
#So, to use it I must convert my chromosome numbers into numeric values with the as.numeric() function which turns input values into numeric values. 
input_Manhattan$CHR <- as.numeric(input_Manhattan$CHR)
#The function is returning an error message. Could be related to NA values
summary(input_Manhattan)
# I have to remove NA values.
input_Manhattan <- input_Manhattan[complete.cases(input_Manhattan$pValues_wilcox),]

#Finally the plot. The manhattan() function 
pdf("Manhattan.pdf",height=10,width=15)
manhattan(input_Manhattan, snp="IlmnID",chr="CHR", bp="MAPINFO", p="pValues_wilcox",col=rainbow(24))
dev.off()
#Both the Manhattan plot and the Volcano plot have untypical looks to them. 
#This is because a non-parametric test was used to test for significance. These tests tend to give discrete p-values in data sets such as this one, so a typical Manhattan plot or Volcano plot is impossible to obtain from these values since there isn't a continuous distribution of p-values.


###STEP 11: Create a heat map of top 100 differentially mehtylated probes 
#A heat map can be created using the heatmap.2() function from the gplots library.
library(gplots)
?heatmap.2
#This function requires a matrix containing beta values of our samples as input. 
#Since I need to produce a heat map of top 100 differentially methylated probes, I will sort the data frame containing probes with significant differential methylation and their beta values. After that, I will extract the top 100 rows into a new matrix.
sorted_significants <- significant[order(significant$pValues_wilcox, decreasing = TRUE),]
input_heatmap <- as.matrix(sorted_significants[1:100,1:8])
head(input_heatmap)
#I need to create color coding for the heat map. The order of colors in the vector will correspond to the order of WT and MUT samples.
#MUT samples will be labeled as orange, and WT as green
colorbar <- c("green","orange","green","green","orange","orange","green","orange")
#Parameters of the heatmap.2() function are set to default values. This means distance will be calculated as euclidean distance, and the hierarchical clustering is complete.
pdf("Heat-map.pdf",height=10,width=7)
heatmap.2(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F,main="Complete linkage")
dev.off()
#The dendrogram which I obtained has 2 clusters, corresponding to MUT and WT samples. So, the clustering seems to be correct.  
#The basis of the clustering can perhaps be in two places - bottom area of the heat map and the are just above the middle. 
#It's evident that CpG sites targeted by probes clustered in the bottom of the probe dendrogram are less methylated in the MUT cluster, than in the WT cluster.
#Furthermore, just above the middle part is noticeable that the WT cluster has mostly very high beta values of corresponding probes, while the MUT cluster has middle-range beta values.

