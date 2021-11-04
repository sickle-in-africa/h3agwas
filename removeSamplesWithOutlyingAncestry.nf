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
    rebuildCovariatesReport;
} from "${projectDir}/modules/base.nf"

include {
    checkInputParams;
    getInputChannels;
    convertCohortDataToEigensoftFormat;
    mapToPrincipalComponents;
    getNumberOfAvailableComponents;
    selectNumbersSmallerThanLimit;
    testPrincipalComponentToPhenotypeAssociations;
    concatenateAcrossPrincipalComponents;
    calculateVarianceAddedByNewPrincipalComponent;
    drawScreePlotForPrincipalComponents;
    getNumberOfSignificantlyAssociatedPrincipalComponents;
    addSignificantPrincipleComponentsToCovariates;
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
        = mapToPrincipalComponents(eigensoftCohortData)

    numberOfAvailableComponents \
        = getNumberOfAvailableComponents(principalComponents)

    availableComponentIds \
        = selectNumbersSmallerThanLimit(
            principalComponentIds,
            numberOfAvailableComponents)

    pvalueAndVarianceTuples \
        = testPrincipalComponentToPhenotypeAssociations(
            availableComponentIds
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

    extendedCovariatesReport \
        = addSignificantPrincipleComponentsToCovariates(
            numberOfSignificantPrincipalComponents,
            principalComponents,
            covariatesReport)

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
            extendedCovariatesReport,
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