#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { SOMASV_CALLER } from './workflows/somasvcaller.nf'

//
// WORKFLOW: Run main ensemblesomasvcaller pipeline
//
workflow {
    SOMASV_CALLER ()
}