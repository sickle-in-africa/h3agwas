include {
    userEmailAddressIsProvided;
    checkOutputDir;
    checkCohortName;
    checkEmailAdressProvided;
    getBasicEmailSubject;
    getBasicEmailMessage;
    getCohortData;
    getCovariatesReport;
    checkInputCohortData;
    checkAssociationInput;
    checkCovariatesReport;
} from "${projectDir}/modules/base.nf"

def checkInputParams() {
    checkCohortName()
    checkOutputDir()
    checkEmailAdressProvided()
    checkAssociationInput()
    checkInputCohortData(params.associationInput)
    checkCovariatesReport()
}

def getInputChannels() {
    return [
        getCohortData(params.associationInput),
        getCovariatesReport(params.associationInput)]
}

process getAssociationReport {
    label 'plink'

    tag "${params.associationInput}CohortData"

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
    output:
        publishDir "${params.outputDir}/${params.associationInput}/reports", mode: 'copy'
        path "${params.cohortName}.assoc"
    script:
        """
        plink \
        --keep-allele-order \
        --bfile ${cohortBed.getBaseName()} \
	    --assoc \
	    --maf 0.01 \
	    --out ${params.cohortName}
	"""
}

process fitGenotypicAssociationModel {
    label 'plink2'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam), path(cohortCov)
    output:
        publishDir "${params.outputDir}/${params.associationInput}/reports", mode: 'copy'
        path "${cohortBed.getBaseName()}.genotypic.PHENO1.glm.logistic"
    script:
        """
        plink2 \
            --bfile ${cohortBed.getBaseName()} \
            --logistic genotypic \
            --covar ${cohortCov} \
            --out ${cohortBed.getBaseName()}.genotypic
        """
}

process fitHethomAssociationModel {
    label 'plink2'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam), path(cohortCov)
    output:
        publishDir "${params.outputDir}/${params.associationInput}/reports", mode: 'copy'
        path "${cohortBed.getBaseName()}.hethom.PHENO1.glm.logistic"
    script:
        """
        plink2 \
            --bfile ${cohortBed.getBaseName()} \
            --logistic hethom \
            --covar ${cohortCov} \
            --out ${cohortBed.getBaseName()}.hethom
        """
}

process fitDominantAssociationModel {
    label 'plink2'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam), path(cohortCov)
    output:
        publishDir "${params.outputDir}/${params.associationInput}/reports", mode: 'copy'
        path "${cohortBed.getBaseName()}.dominant.PHENO1.glm.logistic"
    script:
        """
        plink2 \
            --bfile ${cohortBed.getBaseName()} \
            --logistic dominant \
            --covar ${cohortCov} \
            --out ${cohortBed.getBaseName()}.dominant
        """
}

process fitRecessiveAssociationModel {
    label 'plink2'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam), path(cohortCov)
    output:
        publishDir "${params.outputDir}/${params.associationInput}/reports", mode: 'copy'
        path "${cohortBed.getBaseName()}.recessive.PHENO1.glm.logistic"
    script:
        """
        plink2 \
            --bfile ${cohortBed.getBaseName()} \
            --logistic recessive \
            --covar ${cohortCov} \
            --out ${cohortBed.getBaseName()}.recessive
        """
}

process fitHetonlyAssociationModel {
    label 'plink2'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam), path(cohortCov)
    output:
        publishDir "${params.outputDir}/${params.associationInput}/reports", mode: 'copy'
        path "${cohortBed.getBaseName()}.hetonly.PHENO1.glm.logistic"
    script:
        """
        plink2 \
            --bfile ${cohortBed.getBaseName()} \
            --logistic hetonly \
            --covar ${cohortCov} \
            --out ${cohortBed.getBaseName()}.hetonly
        """
}

process drawManhattanPlot {
    label 'qqman'

    tag "associationReport"

    input:
        path associationReport
    output:
        publishDir "${params.outputDir}/${params.associationInput}/plots", mode: 'copy'
        path manhattanPlot
    script:
        manhattanPlot = "${params.cohortName}.manhattan.png"
        """
        #!/usr/bin/env Rscript --vanilla
        library(qqman)
        assoc <- read.table("${associationReport}", header=TRUE)
        png("${manhattanPlot}", width = 1000, height = 1000, units = "px")
        manhattan(
	    assoc,
	    chr="CHR",
	    bp="BP",
	    snp="SNP",
	    p="P",
	    logp=TRUE,
	    ylim = c(0, 8))
        dev.off()
        """
}

process drawQqPlot {
    label 'qqman'
	
    tag "associationReport"

    input:
        path associationReport
    output:
        publishDir "${params.outputDir}/${params.associationInput}/plots", mode: 'copy'
        path qqplot
    script:
        qqplot = "${params.cohortName}.qqplot.png"
        """
        #!/usr/bin/env Rscript --vanilla
        library(qqman)
        assoc <- read.table("${associationReport}", header=TRUE)
        png("${qqplot}", width = 1000, height = 1000, units = "px")
        qq(assoc\$P)
        dev.off()
        """
}
