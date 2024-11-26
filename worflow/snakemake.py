#################################
#### Short reads co-assembly ####
#################################


#################################
#### Import Python libraries ####
#################################
import glob


#################################
#### Define Global Variables ####
#################################
SAMPLE,RF=glob_wildcards("../resources/Data/{sample}_{rf}.fastq.gz")


#########################
#### General Outputs ####
#########################
rule all:
    input: 
        expand("../resources/Outputs/fastqc_rr/{sample}_{rf}.html", sample=SAMPLE, rf=RF), 
        expand("../resources/Outputs/fastqc_rr/{sample}_{rf}_fastqc.zip", sample=SAMPLE, rf=RF),
        expand("../resources/Outputs/trimmed_reads/{sample}_R{pe}.fastq.gz", sample=SAMPLE, pe=["1","2"]), 
        expand("../results/fastqc_tr/{sample}_{rf}.html", sample=SAMPLE, rf=RF),
        expand("../resources/Outputs/kraken2_rr/reports/{sample}.txt", sample=SAMPLE),
        expand("../results/bracken/{sample}.txt", sample=SAMPLE),
        expand( "../resources/Outputs/mergepe/{sample}.fastq.gz", sample=SAMPLE)

#################################
#### Quality Check raw reads ####
#################################
rule fastqc_rr:
    input:
        reads="../resources/Data/{sample}_{rf}.fastq.gz"
    output: 
        html="../resources/Outputs/fastqc_rr/{sample}_{rf}.html",
        zip="../resources/Outputs/fastqc_rr/{sample}_{rf}_fastqc.zip"
    params:
        extra="--quiet"
    threads:
        4
    resources:
        mem_mb=4000
    log:
        "../resources/Logs/fastqc_rr/{sample}_{rf}.log"
    conda:
        "../envs/fastqc_env.yaml"
    wrapper:
        "v4.5.0/bio/fastqc"


###################################
#### Quality Control raw reads ####
###################################
rule trimmomatic_rr:
    input:
        r1="../resources/Data/{sample}_R1.fastq.gz",
        r2="../resources/Data/{sample}_R2.fastq.gz"
    output:
        r1="../resources/Outputs/trimmed_reads/{sample}_R1.fastq.gz", 
        r2="../resources/Outputs/trimmed_reads/{sample}_R2.fastq.gz", 
        r1_unpaired="../resources/Outputs/trimmomatic/{sample}_R1.unpaired.fastq.gz", 
        r2_unpaired="../resources/Outputs/trimmomatic/{sample}_R2.unpaired.fastq.gz"
    log:
        "../resources/Logs/trimmomatic/{sample}.log"
    conda:
        "../envs/trimmomatic_env.yaml"
    params:
        trimmer=["TRAILING:15", 
                 "LEADING:15",
                 #"SLIDINGWINDOW:4:20",
                 "ILLUMINACLIP:../resources/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads",
                 "MINLEN:50"], #Need to add more arguments
        compression_level="-5"
    threads: 
        4
    resources:
        mem_mb=3024
    wrapper:
        "v4.5.0/bio/trimmomatic/pe"

#####################################
#### Quality check trimmed reads ####
#####################################
rule fastqc_tr:
    input:
        reads="../resources/Outputs/trimmed_reads/{sample}_{rf}.fastq.gz"
    output: 
        html="../results/fastqc_tr/{sample}_{rf}.html",
        zip="../results/fastqc_tr/{sample}_{rf}_fastqc.zip"
    params:
        extra="--quiet"
    threads:
        4
    resources:
        mem_mb=4000
    log:
        "../resources/Logs/fastqc_tr/{sample}_{rf}.log"
    conda:
        "../envs/fastqc_env.yaml"
    wrapper:
        "v4.5.0/bio/fastqc"
