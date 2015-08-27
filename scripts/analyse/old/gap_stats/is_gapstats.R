# This script has been written with RStudio.
# It computes a graph of the number of reads matching an intron size for each repetition (3).
# It uses 3 files.
# To use it: only file1, file2 and file3 should be changed.

# Files
file1="/home/drahull/Master/stage_M1/Data_test/J1onJ_star/is/is.gap_stats.txt";
file2="/home/drahull/Master/stage_M1/Data_test/J2onJ_star/is/is.gap_stats.txt";
file3="/home/drahull/Master/stage_M1/Data_test/J3onJ_star/is/is.gap_stats.txt";

# Files content
valuesF1="nbReadsJ1";
valuesF2="nbReadsJ2";
valuesF3="nbReadsJ3";


#Mapping tool used for the have these files informations
#toolName="STAR"; # Name of the tool
maxIntronSizeOfTool=1000000; # The max intron size set in the tool (by default: tophat2 is 500Kb, STAR 1Mb and CRAC 300Kb)
#toolComment="";

# Work
setwd("/home/drahull/Master/stage_M1/Scripts_R/gap_stats"); # Set the work folder


# Reading files as tables, skipping the first line and the second as header
J1table=read.table(file1, header=TRUE, sep="\t", skip=1);
J2table=read.table(file2, header=TRUE, sep="\t", skip=1);
J3table=read.table(file3, header=TRUE, sep="\t", skip=1);

#Regroup value on log10()=0.1 scapes to smooth graph. log10 => 0 to 7.
nbStepLog=71;

#Jvalues is an array of 2-Dim with maxStepLog values each (default value : 0). The first is the higher size of introns the second the number of reads of this introns
Jvalues=array(0, dim=c(nbStepLog,4), dimnames=list(NULL,c("sizes",valuesF1,valuesF2,valuesF3)));

#Initialization of some works values
i=0; #iterator of the log10
stepLog=0.1; #the step value of the iterator
for(x in 1:nbStepLog){Jvalues[x,"sizes"]=i; i=i+stepLog;} #Initialize log10 sizes

j=1; #Iterator of the table Jvalues;

s1=1; #iterator of the table J1
s2=1; #iterator of the table J2
s3=1; #iterator of the table J3

length_s1=length(J1table$size);
length_s2=length(J2table$size);
length_s3=length(J3table$size);

# If log10val of the size of the intron < at the actual, we add the number of reads for each Jx.
z=1;
while(z <= nbStepLog){
  if(z < nbStepLog){
    #J1
    while(s1 <= length_s1 && log10(J1table$size[s1]) <= Jvalues[z,"sizes"]){
      Jvalues[z,valuesF1] = Jvalues[z,valuesF1] + J1table$reads[s1];
      s1=s1+1;
    }
    #J2
    while(s2 <= length_s2 && log10(J2table$size[s2]) <= Jvalues[z,"sizes"]){
      Jvalues[z,valuesF2] = Jvalues[z,valuesF2] + J2table$reads[s2];
      s2=s2+1;
    }
    #J3
    while(s3 <= length_s3 && log10(J3table$size[s3]) <= Jvalues[z,"sizes"]){
      Jvalues[z,valuesF3] = Jvalues[z,valuesF3] + J3table$reads[s3];
      s3=s3+1;
    }
  }
  #Add the rest of s1 s2 s3 for the last size
  else{
    #J1
    while(s1 <= length_s1){
      Jvalues[z,valuesF1] = Jvalues[z,valuesF1] + J1table$reads[s1];
      s1=s1+1;
    }
    #J2
    while(s2 <= length_s2){
      Jvalues[z,valuesF2] = Jvalues[z,valuesF2] + J2table$reads[s2];
      s2=s2+1;
    }
    #J3
    while(s3 <= length_s3){
      Jvalues[z,valuesF3] = Jvalues[z,valuesF3] + J3table$reads[s3];
      s3=s3+1;
    }
  }
  z=z+1;
}

#Get some basics informations
nbReadsMax=Jvalues[which.max(Jvalues)];
nbReadsMean=mean(c(Jvalues[,valuesF1],Jvalues[,valuesF2],Jvalues[,valuesF3]))

#Create the plot of J1
textSize=0.5 # cex values
plot.default(Jvalues[,"sizes"], Jvalues[,valuesF1], type='l', col="blue", log="y", cex=textSize);
  #Add J2 and J3
lines(Jvalues[,"sizes"], Jvalues[,valuesF2], col="orange");
lines(Jvalues[,"sizes"], Jvalues[,valuesF3], col="green");
  #Add the vertical line of the limite of the tool introns
maxIntronToolLog=log10(maxIntronSizeOfTool);
abline(v=maxIntronToolLog, col="red");
  # Add legend
  #lines
legend(x="topleft",inset=c(0,0), c(valuesF1, valuesF2, valuesF3),lty=c(1,1,1), lwd=c(2,2,2), col=c("blue","orange","green"), cex=textSize);
  #vertical line label
labelPosY=10^(log10(nbReadsMax)/1.1);
text(maxIntronToolLog, labelPosY," Introns size\n limit of the\n tool", cex=textSize, col="red", adj=0);

