#!/bin/bash
#$ -S /bin/bash
#$ -j y
#$ -N bam2cram
#$ -l h_rt=100:00:00
#$ -l h_vmem=20G
#$ -cwd
#$ -l mem_free=20G

#This script takes a bam file or a list of bam files (filename must end ".list") and creates cram files with GATK PrintReads 
#    InpFil - (required) - Path to Bam file. Alternatively a file with a list of bams can be provided and the script run as an array job. List file name must end ".list"
#    RefFil - (required) - shell file containing variables with locations of reference files, jar files, and resource directories; see list below for required variables
#    LogFil - (optional) - File for logging progress
#    Flag - B - BadET - prevent GATK from phoning home
#    Flag - F - Fix mis-encoded base quality scores - see GATK manual. GATK will subtract 31 from all quality scores; used to fix encoding in some datasets (especially older Illumina ones) which starts at Q64 (see https://en.wikipedia.org/wiki/FASTQ_format#Encoding)
#    Help - H - (flag) - get usage information

#list of required vairables in reference file:
# $REF - reference genome in fasta format - must have been indexed using 'bwa index ref.fa'
# $GATK - GATK jar file 
# $ETKEY - GATK key file for switching off the phone home feature, only needed if using the B flag

#list of required tools:
# java <http://www.oracle.com/technetwork/java/javase/overview/index.html>
# GATK <https://www.broadinstitute.org/gatk/> <https://www.broadinstitute.org/gatk/download>

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
(-t <X>-<Y> [if providing a list]) CramCreate.sh -i <BAMFiles> -r <referencefile> -l <logfile> -FBH

     -i (required) - Path to tumor Bam file for tumor variant calling or \".list\" file containing a multiple paths
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -l (optional) - Log file
     -B (flag) - Prevent GATK from phoning home
     -F (flag) - Fix mis-encoded base quality scores - see GATK manual
     -H (flag) - echo this message and exit
"

BadET="false"
FixMisencoded="false"

while getopts i:r:a:l:BFH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";; 
        a) ArrNum="$OPTARG";; 
        l) LogFil="$OPTARG";;
        B) BadET="true";;
        F) FixMisencoded="true";;
        H) echo "$usage"; exit;;
  esac
done

#check all required paramaters present
if [[ ! -e "$InpFil" ]] || [[ ! -e "$RefFil" ]]; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi

if [[ -z "${ArrNum}" ]]
then
    ArrNum=$SGE_TASK_ID
fi

#Call the RefFil to load variables
RefFil=`readlink -f $RefFil`
source $RefFil

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func" #library functions begin "func"

#Set local Variables
funcGetTargetFile

#Sample inputs and names
InpFil=`readlink -f $InpFil` #resolve absolute path to bam
BamFil=$(tail -n+$ArrNum $InpFil | head -n 1) 
BamNam=`basename $BamFil | sed s/.bam//`
BamNam=${BamNam/.bam/} # a name for the output files

#CRAM and log files named after sample
if [[ -z $LogFil ]]; then LogFil=$BamNam.cram.log; fi # a name for the log file
CramFil=$BamNam.cram #Output File
GatkLog=$BamNam.cram.gatklog #a log for GATK to output to, this is then trimmed and added to the script log
TmpLog=$BamNam.cram.temp.log #temporary log file
TmpDir=$BamNam.cram.tempdir; mkdir -p $TmpDir #temporary directory

echo "Reference Genome File is $REF"

#Start Log File
ProcessName="CRAM from BAM generatation with GATK" # Description of the script - used in log
funcWriteStartLog

##Run CRAM generation
StepName="CRAM generation with GATK"
StepCmd="java -Xmx16G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=$TmpDir -jar $GATKJAR PrintReads 
 -R $REF
 -I $BamFil
 -O $CramFil" #command to be run
funcGatkAddArguments # Adds additional parameters to the GATK command depending on flags (e.g. -B or -F)
echo $StepCmd
funcRunStep

#End Log
funcWriteEndLog

#End
