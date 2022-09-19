version 1.0

#Mohammed Alhusayan
#10/20/2021


##Normal tasks done in almost all pipelines

## Quality check of fastq files before proceeding with the analysis.
task QC_Trimming {
    input {
        Array[File] fastqs
        String outdir
        Int threads}

    command <<<

    FASTQC=~{outdir}/FASTQC

        ##Check if the directory exists, otherwise create the directory.
        [ -d $FASTQC ] && echo $(date) "===================> FASTQC Directory Exists" | tee -a out.log \
        || echo $(date) "===================> FASTQC directory does not exist, Creating FASTQC Directory" | tee -a out.log \
        && mkdir -p  $FASTQC

        fastqc ~{sep =' ' fastqs}  -t ~{threads} -o $FASTQC  > $FASTQC/FASTQC.log 2>&1
        trim_galore --paired --illumina --fastqc \
        -o $FASTQC  ~{sep =' ' fastqs}
    >>>

    parameter_meta {
        fastqs: "Fastq files to be quality checked."
        outdir: "Directory to store the output files from FASTQC."
        threads: "Number of threads to use for this analysis."
    }

    meta {
        description: "This tasks performs a quality check on the fastq files."
        author: "Mohammed Alhusayan"
        email: "mohammed.alhusayan@pennmedicine.upenn.edu"
    }


    output {
        Array[File]+ cleaned_fastqs = glob("${outdir}/FASTQC/*.gz")
    }

}

task align_reads {
    input {
        Array[File]+ fastqs
        String? genomeDir
        Int threads
        String outdir
        String? genomeName
    }

    command <<<

        source activate ngsmo
        name=$(basename ~{fastqs[0]} | cut -d_ -f1-3)
        BOWTIE=~{outdir}/BOWTIE


        ##Create Bowtie directory if it does not exist.
        [ -d $BOWTIE ] && echo $(date) "===================> BOWTIE Directory in ${name} Exists" | tee -a out.log \
        || echo $(date) "===================> BOWTIE directory in ${name} does not exist, Creating BOWTIE Directory" | tee -a out.log \
        && mkdir -p  $BOWTIE


        ## Align reads and remove duplicates, skip SAM file and directly pipe the output of BOWTIE2 into samtools
        ## To get the bam file
        bowtie2 -x ~{genomeDir} -U ~{sep =',' fastqs}  -p ~{threads} 2>$BOWTIE/BOWTIE.log \
        | samtools view -b -F 3340  -q 30 -o ${BOWTIE}/${name}_~{genomeName}_filtered.bam


        ##  Sort alignments by coordinates on the genome.
        samtools sort -@ ~{threads} -o ${BOWTIE}/${name}_~{genomeName}_filtered_sorted.bam ${BOWTIE}/${name}_~{genomeName}_filtered.bam


        ## Index the bam file for fast random access for further analysis.
        samtools index -@ ~{threads} ${BOWTIE}/${name}_~{genomeName}_filtered_sorted.bam



        conda deactivate
    >>>


    output {
        Array[File] Bam = glob("${outdir}/BOWTIE/*filtered_sorted.bam")
    }
    runtime {
        cpu: threads
    }
    ##Add docker image after installation of docker.

}

task make_bigwig {
    input {
        File bam
    }

    command <<<



    >>>
}


##Variant Calling ##
#############################################


############Data_PreProcessing#############
task align_reads_vc {
    input {
        String fastq_dir
        Int threads
        String outdir
        String refgenome
        String name
    }

    command <<<
        source activate ngsmo
        out_dir=~{outdir}/bams
        [ -d $out_dir ] || mkdir $out_dir

        bwa mem -t ~{threads} ~{refgenome} ~{fastq_dir}/*.gz |  samtools view -b -q20 | samtools sort > ${out_dir}/~{name}.bam
        conda deactivate
    >>>

    output {
        File bam = outdir + "/bams/" + name + ".bam"

}

}


task markDup {
    input {
        File bam
        String outdir
        String name
    }

    command <<<

        [ -d ~{outdir}/bams/dup_stats ] || mkdir ~{outdir}/bams/dup_stats
        gatk MarkDuplicates -I ~{bam} -O  ~{outdir}/bams/~{name}_dedup.bam -M ~{outdir}/bams/dup_stats/~{name}_dup_stats.txt

    >>>

    output {
        File dedup_bam = outdir + "/bams/" + name + "_dedup.bam"
    }
}


task add_rg {
    input {
        String name
        String outdir
        File bam
        String LB
        String PU
        Int ID

    }

    command <<<
        gatk AddOrReplaceReadGroups -I ~{bam} -O ~{outdir}/bams/~{name}_dedup_ann.bam \
        -PL ILLUMINA -LB ~{LB} -SO coordinate -SM ~{name} -PU ~{PU}  --CREATE_INDEX True -ID ~{ID}

    >>>


    output {
        File ann_bam = outdir + "/bams/" + name + "_dedup_ann.bam"
    }
}


##BQSR##
##Base quality score recalibration (BQSR) is a process in which we apply machine learning to model these errors empirically and adjust
##the quality scores accordingly. For example we can identify that, for a given run, whenever we called two A nucleotides in a row, the next base
##we called had a 1% higher rate of error. So any base call that comes after AA in a read should have its quality score reduced by 1%.
##We do that over several different covariates (mainly sequence context and position in read, or cycle) in a way that is additive.
##So the same base may have its quality score increased for one reason and decreased for another.
task BQSR {

    input {
        String refgenome
        File indels
        File snps
        String outdir
        String name
        File bam

        ## Inputs below are not supplied to the actual command, it is supplied as inputs to the workingdir
        File snpsIndex
        File indelIndex
    }

    command <<<

    
    gatk BaseRecalibrator -R ~{refgenome} --known-sites ~{indels} --known-sites ~{snps}  -I ~{bam} -O ~{outdir}/bams/~{name}_recal.table



        >>>


    output {

        File recalibration_table =  outdir + "/bams/" + name + "_recal.table"

    }


}

task Apply_BQSR {
    input {
        File recalibration_table
        File bam
        String outdir
        String name
    }

    command <<<
    
        gatk ApplyBQSR -I ~{bam}  -bqsr ~{recalibration_table} -O ~{outdir}/bams/~{name}_recal.bam
    
    >>>


    output {
        File recal_bam = outdir +"/bams/" + name + "_recal.bam"
    }
}

task Mutect2 {
    input {
        File refgenome
        Array[File] bams
        String outdir
    }
}