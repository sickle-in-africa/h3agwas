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
    testPrincipalComponentToPhenotypeAssociations;
} from "${projectDir}/modules/outlyingAncestry.nf"

workflow {

    checkInputParams()

    cohortData = getInputChannels()

    eigensoftCohortData\
        = convertCohortDataToEigensoftFormat(
            cohortData)

    principalComponents = selectPrincipalComponents(eigensoftCohortData)

    testPrincipalComponentToPhenotypeAssociations(
        cohortData,
        principalComponents)

}
