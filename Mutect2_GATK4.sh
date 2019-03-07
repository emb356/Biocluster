#!/bin/bash
#$ -S /bin/bash
#$ -j y
#$ -N Mutect2
#$ -l h_rt=100:00:00
#$ -l h_vmem=30G
#$ -cwd
#$ -l mem_free=30G

#This script takes a bam file or a list of bam files (filename must end ".list") and runs variant calling using GATK4.0's Mutect2
#Can be used on both exome and genome bams (exome needs TgtBed parameter)
#    InpFil - (required) - Path to Tumor Bam file to be aligned. Alternatively a file with a list of bams can be provided and the script run as an array job. List file name must end ".list"
#    NormFil - (required) - Path to Normal Bam file to be aligned. Alternatively a file with a list of bams can be provided (must be in same order as matched tumor sample in InpFil) and the script run as an array job. List file name must end ".list"
#    RefFil - (required) - shell file containing variables with locations of reference files, jar files, and resource directories; see list below for required variables
#    TgtBed - (optional) - Exome capture kit targets bed file (must end .bed for GATK compatability) ; may be specified using a code corresponding to a variable in the RefFil giving the path to the target file- only required if calling pipeline
#    LogFil - (optional) - File for logging progress
#    Flag - B - BadET - prevent GATK from phoning home
#    Flag - F - Fix mis-encoded base quality scores - see GATK manual. GATK will subtract 31 from all quality scores; used to fix encoding in some datasets (especially older Illumina ones) which starts at Q64 (see https://en.wikipedia.org/wiki/FASTQ_format#Encoding)
#    Help - H - (flag) - get usage information

#list of required vairables in reference file:
# $REF - reference genome in fasta format - must have been indexed using 'bwa index ref.fa'
# $MUTECTREF - population allele frequency vcf from GATK
# $EXOMPPLN - directory containing exome analysis pipeline scripts
# $GATK - GATK jar file 
# $ETKEY - GATK key file for switching off the phone home feature, only needed if using the B flag

#list of required tools:
# java <http://www.oracle.com/technetwork/java/javase/overview/index.html>
# GATK <https://www.broadinstitute.org/gatk/> <https://www.broadinstitute.org/gatk/download>
# bgzip 
# tabix

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
(-t <X>-<Y> [if providing a list]) Mutect2_GATK4.sh -i <TumorFile> -n <NormFile> -r <referencefile> -t <targetfile> -l <logfile> -FBH

     -i (required) - Path to tumor Bam file for tumor variant calling or \".list\" file containing a multiple paths
     -n (required) - Path to normal Bam file for normal variant calling or \".list\" file containing a multiple paths
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -t (optional) - Exome capture kit targets or other genomic intervals bed file (must end .bed for GATK compatability)
     -l (optional) - Log file
     -B (flag) - Prevent GATK from phoning home
     -F (flag) - Fix mis-encoded base quality scores - see GATK manual
     -H (flag) - echo this message and exit
"

BadET="false"
FixMisencoded="false"

while getopts i:r:a:n:l:t:BFH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        n) NormFil="$OPTARG";;
        r) RefFil="$OPTARG";; 
        a) ArrNum="$OPTARG";; 
        l) LogFil="$OPTARG";;
        t) TgtBed="$OPTARG";;
        B) BadET="true";;
        F) FixMisencoded="true";;
        H) echo "$usage"; exit;;
  esac
done

#check all required paramaters present
if [[ ! -e "$InpFil" ]] || [[ ! -e "$NormFil" ]] || [[ ! -e "$RefFil" ]]; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi

# add Target Bed param if present
if [[ ! -e "$TgtBed" ]]; then 
    TgtParam=""
else
    TgtParam="-L $TgtBed
              --interval-padding 100"
fi

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

#Tumor sample inputs and names
InpFil=`readlink -f $InpFil` #resolve absolute path to bam
BamFil=$(tail -n+$ArrNum $InpFil | head -n 1) 
BamNam=`basename $BamFil | sed s/.bam//`
BamNam=${BamNam/.bam/} # a name for the output files

#Normal sample inputs and names
NormFil=`readlink -f $NormFil` #resolve absolute path to bam
NBamFil=$(tail -n+$ArrNum $NormFil | head -n 1) 
NBamNam=`basename $NBamFil | sed s/.bam//`
NBamNam=${NBamNam/.bam/} # a name for the output files

#Vcf and log files named after tumor sample
if [[ -z $LogFil ]]; then LogFil=$BamNam.M2VCF.log; fi # a name for the log file
VcfFil=$BamNam.Mutect2.g.vcf #Output File
GatkLog=$BamNam.Mutect2.gatklog #a log for GATK to output to, this is then trimmed and added to the script log
TmpLog=$BamNam.Mutect2.temp.log #temporary log file
TmpDir=$BamNam.Mutect2.tempdir; mkdir -p $TmpDir #temporary directory
Bamout=$BamNam.Mutect2.bam
#infofields="-A AlleleFraction -A BaseQualityRankSumTest -A Coverage -A MappingQualityRankSumTest -A MappingQualityZero -A QualByDepth -A RMSMappingQuality -A FisherStrand -A InbreedingCoeff -A ClippingRankSumTest -A DepthPerSampleHC" #Annotation fields to output into vcf files
infofields="-A Coverage -A DepthPerAlleleBySample -A MappingQuality -A BaseQuality -A StrandArtifact"

echo "Reference Genome File is $REF"

#Start Log File
ProcessName="Somatic VCF generatation with GATK 4.0 Mutect2" # Description of the script - used in log
funcWriteStartLog

##Run VCF generation
StepName="VCF generation with GATK 4.0 Mutect2"
StepCmd="java -Xmx16G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=$TmpDir -jar $GATKJAR Mutect2 
 -R $REF
 $TgtParam
 -I $BamFil
 -tumor $BamNam 
 -I $NBamFil
 -normal $NBamNam
 -O $VcfFil
 --germline-resource $MUTECTREF
 --af-of-alleles-not-in-resource 0.0000025
 $infofields" #command to be run
funcGatkAddArguments # Adds additional parameters to the GATK command depending on flags (e.g. -B or -F)
echo $StepCmd
funcRunStep

##gzip and index the VCF
StepName="gzip and index the VCF"
StepCmd="bgzip -f $VcfFil; tabix -f -p vcf $VcfFil.gz"
funcRunStep
rm $VcfFil.idx

#End Log
funcWriteEndLog

#End
