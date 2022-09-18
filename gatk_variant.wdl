version 1.0

import "../tasks.wdl" as tsk

workflow variant_calling {
    input {
        Array[String] fastq_dirs
        Int threads
        String outdir
        String refgenome
        File indels
        File snps
        File snpsIndex
        File indelIndex

        Array[Pair[Int,String]] zipped = zip(range(length(fastq_dirs)),fastq_dirs)
    }


    scatter (pair in zipped ) {

        ##The way of indexing pairs, using "right" or "left"
        String dir = pair.right
        Int i = pair.left
        String name = basename(dir)
        String LB = "LIB" + i
        String PU = "BRC" + i
        Int ID = i + 1
        
        call tsk.align_reads_vc as align_reads {
        input: fastq_dir=dir,threads=threads, outdir=outdir, refgenome=refgenome, name=name
        }

        call tsk.markDup as mark_duplicates {
            input: bam=align_reads.bam, outdir=outdir,name=name
        }


        call tsk.add_rg as add_rg {
            
            input: name=name, outdir=outdir, bam=mark_duplicates.dedup_bam, LB = LB, PU = PU, ID=ID
        }

        call tsk.BQSR as BQSR {
            input: refgenome=refgenome,bam=add_rg.ann_bam, indels=indels,snps=snps,outdir=outdir,name=name, snpsIndex=snpsIndex,indelIndex=indelIndex
        }

        call tsk.Apply_BQSR as Apply_BQSR {
            input: recalibration_table=BQSR.recalibration_table, bam = add_rg.ann_bam,outdir=outdir,name=name
        }

    }
}


