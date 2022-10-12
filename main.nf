#!/usr/bin/env nextflow


/* 
vim: syntax=groovy
-*- mode: groovy;-*-

ARTIC-VIDRL workflow

*/

nextflow.enable.dsl=2


include { validate_primer_scheme } from './modules/utils'
include { get_fastq_files } from './modules/utils'

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
    (primer_scheme, primer_bed) = validate_primer_scheme(params.scheme_dir)
    println("Primer scheme directory: ${primer_scheme[0]} (scheme: ${primer_scheme[1]})")
    

    fastq_files = get_fastq_files(
        params.fastq_gather, 
        params.fastq_id, 
        params.fastq_dir, 
        params.fastq_ext, 
        params.barcodes, 
        params.sample_sheet
    )

    artic_nanoq = ArticNanoq(
        fastq_files
    )
    artic_medaka = ArticMinion(
        artic_nanoq[0], 
        primer_scheme
    )
    artic_coverage = ArticCovtobed(
        artic_medaka[0]
    )

    artic_params = ArticParams(started)

    artic_report = ArticReport(
        artic_coverage | collect,
        artic_medaka[1] | collect,
        artic_nanoq[1] | collect,
        primer_bed,
        artic_params
    )

}