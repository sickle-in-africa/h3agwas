nextflow.enable.dsl=2

include {
    printWorkflowExitMessage;
    sendWorkflowExitEmail;
    rebuildCovariatesReport;
} from "${projectDir}/modules/base.nf"

include {
    checkInputParams;
    getInputChannels;
    selectAutosomes;
    indexByChromosome;
    decompressGeneticMapsArchive;
    selectAutosomalGenotypeSet;
    indexReferencePanel;
    concatenateFiles;
    selectBiallelicSnvsWithBcftools;
    alignWithConformGt;
    selectOppositeStrandSnvs;
    rebuildCohortData;
    flipGenotypesOnOppositeStrand;
    indexWithTabix;
    concatenateWithBcftools;
} from "${projectDir}/modules/strandErrors.nf"


workflow {

    checkInputParams()

    (inputCohortData,
     referencePanels,
     geneticMapsArchive,
     covariatesReport) \
        = getInputChannels()

    referencePanelsWithIndexes \
        = indexReferencePanel(
            referencePanels)

    geneticMaps \
        = decompressGeneticMapsArchive(geneticMapsArchive)
            .flatten()

    autosomalGeneticMaps \
        = indexByChromosome(selectAutosomes(geneticMaps))

    filteredReferencePanels \
        = selectBiallelicSnvsWithBcftools(
        referencePanelsWithIndexes)

    genotypeSet \
        = selectAutosomalGenotypeSet(
        inputCohortData)

    (alignedGenotypeSubsets,
     alignmentLogFiles) \
        = alignWithConformGt(
        genotypeSet.combine(filteredReferencePanels))

    alignmentLogFiles | view()

    listsOfSnvsOnTheOppositeStrand \
        = selectOppositeStrandSnvs(
            alignmentLogFiles)

    snvsOnTheOppositeStrand \
        = concatenateFiles(
            listsOfSnvsOnTheOppositeStrand.collect())

    indexedGenotypeSubsets \
        = indexWithTabix(
            alignedGenotypeSubsets)

    alignedGenotypeSet \
        = concatenateWithBcftools(
            indexedGenotypeSubsets.collect())

    alignedCohortData \
        = rebuildCohortData(
            alignedGenotypeSet,
            inputCohortData)

    flippedCohortData \
        = flipGenotypesOnOppositeStrand(
            alignedCohortData.combine(snvsOnTheOppositeStrand))

    filteredCovariatesReport \
        = rebuildCovariatesReport(
            'strandFlipped',
            covariatesReport,
            flippedCohortData)
}

workflow.onComplete {
    printWorkflowExitMessage()
    sendWorkflowExitEmail()
}