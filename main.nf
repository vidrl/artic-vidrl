#!/usr/bin/env nextflow


/* 
vim: syntax=groovy
-*- mode: groovy;-*-

ARTIC-VIDRL workflow

*/

nextflow.enable.dsl=2


include { validate_primer_scheme } from './modules/utils'
include { get_fastq_gather } from './modules/utils'
include { get_fastq_dirs } from './modules/utils'
include { get_samples } from './modules/utils'

include { ArticNanoq } from './modules/artic'
include { ArticMinion } from './modules/artic'
include { ArticCovtobed } from './modules/artic'
include { ArticReport } from './modules/artic'
include { ArticParams } from './modules/artic'

workflow {

    started = String.format('%tF %<tH:%<tM', java.time.LocalDateTime.now())

    println("""
    Pipeline parameters:

    outdir:           $params.outdir
    version           $params.version
    
    sample_sheet:     $params.sample_sheet
    fastq_gather:     $params.fastq_gather
    fastq_id:         $params.fastq_id
    fastq_dir:        $params.fastq_dir
    fastq_ext:        $params.fastq_ext
    barcodes:         $params.barcodes

    scheme_dir:       $params.scheme_dir
    medaka_model:     $params.medaka_model
    min_length:       $params.min_length
    max_length:       $params.max_length
    min_quality:      $params.min_quality
    normalise:        $params.normalise
    report_title:     $params.report_title
    """)

    if (!params.medaka_model){
        println("Please provide a Medaka model (--medaka_model)")
        System.exit(1)
    }

    if (!params.scheme_dir){
        println("Please provide a primer scheme directory (--scheme_dir)")
        System.exit(1)
    }
    primer_scheme = validate_primer_scheme(params.scheme_dir)
    println("Primer scheme directory: ${primer_scheme[0]} (scheme: ${primer_scheme[1]})")
    
    primer_scheme_minion = tuple(primer_scheme[0], primer_scheme[1])
    primer_scheme_bed = primer_scheme[2]

    if (params.fastq_gather) {
        println("Collecting read files: $params.fastq_gather")
        if (!params.fastq_id){
            println("Please provide an output identifier for the sample (--fastq_id)")
            System.exit(1)
        }
        fastq_files = get_fastq_gather(params.fastq_gather)
    } else if (params.fastq_dir) {
        println("Collecting read files with extension '$params.fastq_ext' in subdirectories of: $params.fastq_dir")
         if (!params.fastq_ext){
            println("Please provide a file extension to gather read files in subdirectories (--fastq_ext)")
            System.exit(1)
        }
        if (params.sample_sheet){
            println("Collecting read files from sample sheet: $params.sample_sheet")
            // Subset barcodes by sample sheet
            fastq_files = get_samples(params.fastq_dir, params.fastq_ext, params.sample_sheet)
        } else {
            // Collect all barcodes
            fastq_files = get_fastq_dirs(params.fastq_dir, params.fastq_ext, params.barcodes)
        }

    } else {
        println("Please provide read file inputs (--fastq_gather | --fastq_dir)")
        System.exit(1)
    }

    fastq_files | view

    artic_nanoq = ArticNanoq(
        fastq_files
    )
    artic_medaka = ArticMinion(
        artic_nanoq[0], 
        primer_scheme_minion
    )
    artic_coverage = ArticCovtobed(
        artic_medaka[0]
    )

    artic_params = ArticParams(started)

    artic_report = ArticReport(
        artic_coverage | collect, // coverage
        artic_medaka[1] | collect, // masks
        artic_nanoq[1] | collect, // nanoq file tuples
        primer_scheme_bed,
        artic_params
    )

    

}