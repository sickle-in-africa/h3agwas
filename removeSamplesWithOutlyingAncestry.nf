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
    concatenateAcrossPrincipalComponents;
    calculateVarianceAddedByNewPrincipalComponent;
    drawScreePlotForPrincipalComponents;
    getNumberOfSignificantlyAssociatedPrincipalComponents;
} from "${projectDir}/modules/outlyingAncestry.nf"

workflow {

    checkInputParams()

    (cohortData,
     principalComponentIds) \
        = getInputChannels()

    eigensoftCohortData \
        = convertCohortDataToEigensoftFormat(
            cohortData)

    principalComponents \
        = selectPrincipalComponents(eigensoftCohortData)

    pvalueAndVarianceTuples \
        = testPrincipalComponentToPhenotypeAssociations(
            principalComponentIds
                .combine(principalComponents
                    .combine(cohortData)))

    pvalueAndVarianceTable \
        = concatenateAcrossPrincipalComponents(
            pvalueAndVarianceTuples.collect())

    pvalueAndVarianceAndVarianceAddedTable \
        = calculateVarianceAddedByNewPrincipalComponent(
            pvalueAndVarianceTable)

    drawScreePlotForPrincipalComponents(
        pvalueAndVarianceAndVarianceAddedTable) | view()

    getNumberOfSignificantlyAssociatedPrincipalComponents(
        pvalueAndVarianceAndVarianceAddedTable)


}
