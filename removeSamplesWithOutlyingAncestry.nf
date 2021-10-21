#!/usr/bin/env nextflow
/*
 *  REMOVE SAMPLES WITH OUTLYING ANCESTRY
 *  =====================================
 *
 *******************************************************************/
nextflow.enable.dsl = 2

include {
    printWorkflowExitMessage;
    collectPlotsTogetherAndZip;
    sendWorkflowExitEmailWithPlots;
} from "${projectDir}/modules/base.nf"

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
    removeOutlyingSamples;
    drawPrincipalComponentPlotForSamples;
    drawPrincipalComponentPlotForOutliers;
    extractOutliers;
    rebuildCohortData;
} from "${projectDir}/modules/outlyingAncestry.nf"

workflow {

    checkInputParams()

    (cohortData,
     principalComponentIds,
     covariatesReport) \
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

    screePlot \
        = drawScreePlotForPrincipalComponents(
            pvalueAndVarianceAndVarianceAddedTable)

    numberOfSignificantPrincipalComponents \
        = getNumberOfSignificantlyAssociatedPrincipalComponents(
            pvalueAndVarianceAndVarianceAddedTable)

    principalComponentsWithoutOutliers \
        = removeOutlyingSamples(
            eigensoftCohortData,
            numberOfSignificantPrincipalComponents)

    outlyingSamples \
        = extractOutliers(principalComponentsWithoutOutliers)

    principalComponentPlotForSamples \
        = drawPrincipalComponentPlotForSamples(
            principalComponents)

    principalComponentPlotForOutliers \
        = drawPrincipalComponentPlotForOutliers(
            principalComponentsWithoutOutliers)

    filteredCohortData \
        = rebuildCohortData(
            cohortData.combine(outlyingSamples))

    filteredCovariatesReport \
        = rebuildCovariatesReport(
            'ancestry',
            covariatesReport,
            filteredCohortData)

    plots = channel
        .empty().mix(
            principalComponentPlotForSamples,
            principalComponentPlotForOutliers)
        .collect()

    collectPlotsTogetherAndZip(
        "ancestry",
        plots)
}

workflow.onComplete {
    printWorkflowExitMessage()
    sendWorkflowExitEmailWithPlots("ancestry")
}