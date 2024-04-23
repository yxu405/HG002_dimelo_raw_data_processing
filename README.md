# HG002_dimelo_raw_data_processing
This script presents steps to process raw pod5 files after the sequencing step from nanopore

~ = /private/groups/migalab/
#This script is dedicated to process dimelo-seq or any methylation based analysis
#basecalling
$ dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v4.2.0 pod5s/ --modified-bases 5mC 6mA > <bam>

#example command on phoenix: 
~/migalab/dan/dorado_v0.5.0/dorado-0.5.0-linux-x64/bin/dorado basecaller \ 
/private/home/jmmenend/software/dorado-0.3.4-linux-x64/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 \  
~/migalab/dan/01_09_24_R1041_DiMeLoAdaptive_CENPC/01_09_24_R1041_DiMeLoAdaptive_CENPC/20240109_1200_6B_PAS52674_0adbae11/pod5_pass \
--modified-bases 5mC 6mA \ 
> ~/migalab/dan/01_09_24_R1041_DiMeLoAdaptive_CENPC/20240109_1200_6B_PAS52674_0adbae11/pod5_pass/01_09_24_R1041_DiMeLoAdaptive_CENPC_5mC_6mA.bam

#Alignment 
$ sbatch alignment_script.slurm

#Filter out mods that are under a certain threshold 
$ python3 ModBamClean_ChrParallel_v0.22.py --bam <bam> --out <name_withouot_.bam> --modqA 250 --modqC 260 
--threads 64

#example command on pheonix 
python3 /private/groups/migalab/dan/script/ModBamClean_ChrParallel_v0.22.py --bam \  
11_28_23_R1041_DiMeLoAdaptive_CENPB_5mC_6mA_winnowmap_sorted_sorted_clean_MD_filtered.bam --out \ 
11_28_23_R1041_DiMeLoAdaptive_CENPB_5mC_6mA_winnowmap_sorted_sorted_clean_MD_filtered_mC --modqA 250 --modqC 260 --threads 64


#clean out secondary reads and sort
samtools view -b -@ 96 -F 2308 <filtered.bam> > <filtered.clean.bam>

samtools sort -@ 10 <filtered.clean.bam> > <filtered.sorted.bam>

samtools index -@ 10 <filtered.sorted.bam>

# At this point, the file is ready for any downstream analysis.
# However, if the dimelo seq package is intended to be run on the exsisting bam file, more processing steps need 
to be performed

# decapitalize MM and ML tags 
#1. convert bam file into sam file 
samtools view -h <bam> -@ 24 > <sam>

#2. change Mm  
sed 's/MM/Mm/g' <sam> > <Mm.sam> 

#3. change Ml
sed 's/ML/Ml/g'<Mm.sam> > <Mm_Ml.sam>

#4. convert file back to bam 
samtools view -@ 128 -S -b <Mm_Ml.sam> > <decapped_bam> 

# find reads that do not contain tags. Using python code here  
python3  
import sys
import pysam

filename = <decapped_bam>
samfile = pysam.AlignmentFile(filename, "rb")
no_Ml=[]
for read in samfile:
    try:
        ml_tag = read.get_tag("Ml")
    except KeyError:
        no_Ml.append(read.query_name)

samfile.close()


with open("reads_to_filter.txt", "w") as file:
    for item in no_Ml:
        # Write each item on a new line
        file.write(f"{item}\n")
        
bash
samtools view -h <bam> | grep -vFf reads_to_filter.txt | samtools view -b -o <final_bam>
