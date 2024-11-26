#################################
#### Short reads co-assembly ####
#################################


#################################
#### Import Python libraries ####
#################################
import glob
import yaml

#################################
#### Define Global Variables ####
#################################
# Define groups and their respective samples
GROUPS = {
    "PrimaryLagoon": ["AP"],
    "MacrophyteRoots": ["AR1", "AR2", "AR3"],
    "Algae": ["M1", "M4", "M9"], 
    "Lemna" : ["M3", "M5", "M8"],
    "FTW" : ["M6", "M7"],
    "SecondaryLagoon": ["AS", "JWC1", "JWC2", "JWF1", "JWF2"], 
    "FTWRoots": ["JR1", "JR2"], 
    "Sediments" : ["JS1", "JS2"]
}

# Extract all samples
SAMPLES = [sample for group in GROUPS.values() for sample in group]
RF = ["R1", "R2"]

#########################
#### General Outputs ####
#########################
rule all:
    input: 
        # FastQC for raw reads
        expand("../resources/Outputs/fastqc_rr/{sample}_{rf}.html", sample=SAMPLES, rf=RF), 
        expand("../resources/Outputs/fastqc_rr/{sample}_{rf}_fastqc.zip", sample=SAMPLES, rf=RF),
        # Trimmed reads
        expand("../resources/Outputs/trimmed_reads/{sample}_R{pe}.fastq.gz", sample=SAMPLES, pe=["1","2"]), 
        # FastQC for trimmed reads
        expand("../results/fastqc_tr/{sample}_R{pe}.html", sample=SAMPLES, pe=["1","2"]),
        expand("../results/fastqc_tr/{sample}_R{pe}_fastqc.zip", sample=SAMPLES, pe=["1","2"]),
        # Co-assembly outputs
        expand("../results/coassembly/{group}.fasta", group=GROUPS.keys())

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
        r1_trimmed="../resources/Outputs/trimmed_reads/{sample}_R1.fastq.gz", 
        r2_trimmed="../resources/Outputs/trimmed_reads/{sample}_R2.fastq.gz", 
        r1_unpaired="../resources/Outputs/trimmomatic/{sample}_R1.unpaired.fastq.gz", 
        r2_unpaired="../resources/Outputs/trimmomatic/{sample}_R2.unpaired.fastq.gz"
    log:
        "../resources/Logs/trimmomatic/{sample}.log"
    conda:
        "../envs/trimmomatic_env.yaml"
    params:
        trimmer=["TRAILING:15", 
                 "LEADING:15",
                 "ILLUMINACLIP:../resources/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads",
                 "MINLEN:50"],
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
        reads="../resources/Outputs/trimmed_reads/{sample}_R{pe}.fastq.gz"
    output: 
        html="../results/fastqc_tr/{sample}_R{pe}.html",
        zip="../results/fastqc_tr/{sample}_R{pe}_fastqc.zip"
    params:
        extra="--quiet"
    threads:
        4
    resources:
        mem_mb=4000
    log:
        "../resources/Logs/fastqc_tr/{sample}_R{pe}.log"
    conda:
        "../envs/fastqc_env.yaml"
    wrapper:
        "v4.5.0/bio/fastqc"

#####################
#### Co-Assembly ####
#####################
rule coassemble:
    input:
        r1=lambda wildcards: ",".join([f"../resources/Outputs/trimmed_reads/{sample}_R1.fastq.gz" for sample in GROUPS[wildcards.group]]),
        r2=lambda wildcards: ",".join([f"../resources/Outputs/trimmed_reads/{sample}_R2.fastq.gz" for sample in GROUPS[wildcards.group]])
    output:
        fasta="../results/coassembly/{group}.fasta"
    log:
        "../resources/Logs/coassembly/{group}.log"
    conda:
        "../envs/spades_env.yaml"
    shell:
        """
        mkdir -p $(dirname {output.fasta}) $(dirname {log})
        spades.py --meta -o ../results/coassembly/{wildcards.group} \
                  -1 {input.r1} -2 {input.r2} \
                  > {log} 2>&1
        """

