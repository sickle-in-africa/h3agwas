#!/usr/bin/env nextflow
/*
 *  REMOVE SAMPLES WITH OUTLYING ANCESTRY
 *  =====================================
 *
 *******************************************************************/
nextflow.enable.dsl = 2

include {
    checkInputParams;
    getInputChannels;
    convertCohortDataToEigensoftFormat;
    selectPrincipalComponents;
} from "${projectDir}/modules/outlyingAncestry.nf"

workflow {

    checkInputParams()

    cohortData = getInputChannels()

    eigensoftCohortData\
        = convertCohortDataToEigensoftFormat(
            cohortData)

    selectPrincipalComponents(eigensoftCohortData)
}
