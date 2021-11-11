include {
    checkCohortName;
    checkOutputDir;
    checkIlluminaGenotypeReports;
    checkIlluminaSampleReport;
    checkIlluminaLocusReport;
    checkClinicalPhenotypeFam;
    checkCovariatesReport;
    userEmailAddressIsProvided;
    checkEmailAdressProvided;
    getBasicEmailSubject;
    getBasicEmailMessage;
} from "${projectDir}/modules/base.nf"

def checkInputParams() {
    checkCohortName()
    checkOutputDir()
    checkEmailAdressProvided()
    checkIlluminaGenotypeReports()
    checkIlluminaSampleReport()
    checkIlluminaLocusReport()
    checkClinicalPhenotypeFam()
    checkCovariatesReport()
}

def getInputChannels() {
    return [
        getIlluminaGenotypeReports(),
        getIlluminaSampleReport(),
        getIlluminaLocusReport(),
        getClinicalPhenotypeReport(),
        getCovariatesReport()]
}

def getIlluminaGenotypeReports() {
    return channel
        .fromPath(params.input.genotypeReports)
}

def getIlluminaSampleReport() {
    return channel
        .fromPath(params.input.sampleReport)
}

def getIlluminaLocusReport() {
    return channel
        .fromPath(params.input.locusReport)
}

def getClinicalPhenotypeReport() {
    return channel
        .fromPath(params.input.clinicalPhenotypeFam)
}

def getCovariatesReport() {
    return channel 
        .fromPath(params.input.covariatesReport)
}

process mergeSampleReports {

    input:
        path sampleReports
    output:
        path "${params.cohortName}-sample-report.csv"
    script:
        """
        awk '(NR == 1) || (FNR > 1)' ${sampleReports} \
            > ${params.cohortName}-sample-report.csv
        """
}

process mergeLocusReports {

    input:
        path locusReports
    output:
        path "${params.cohortName}-locus-report.csv"
    script:
        """
        awk '(NR == 1) || (FNR > 1)' ${locusReports} \
            > ${params.cohortName}-locus-report.csv
        """
}

def splitTextFiles(inputFiles) {
    return inputFiles
        .splitText(
            by: 5000000,
            keepHeader: false,
            file: true,
            compress: false )
}

process convertGenotypeReportToLongFormat {
    label 'smallMemory'
    label 'perl'

    tag "${genotypeReport.getBaseName()}"

    cache 'lenient'

    input:
        path genotypeReport
    output:
        path "*.lgen"
    script:
        template 'convertGenotypeReportToLongFormat.pl'
}

def concatenateLgenFiles(inputLgenFiles) {
    return inputLgenFiles
        .collectFile(
            name: "concatenated.lgen",
            sort: true)
}

process convertLocusReportToMap {
    label 'smallMemory'

    tag "${locusReport.getBaseName()}"

    input:
        path locusReport
    output:
        path "${outputFile}"
    script:
        outputFile = "${locusReport.getBaseName()}.map"
        """
        cut \
            -f2-4 \
            -d',' ${locusReport} \
            | sed 's/,/ /g' \
                | awk \
                    '{print \$2,\$1,"0",\$3}' \
                    | awk '\$1!="0"' \
                        | sed '1d' \
                            | sort \
                                | uniq -u \
            > ${outputFile}
        """
}

process removeSamplesWithFailedGenotypes {
    label 'smallMemory'

    tag "${sampleReport.getBaseName()}"

    input:
        path sampleReport
    output:
        path "${outputFile}"
    script:
        outputFile = "${sampleReport.getBaseName()}.filtered.${sampleReport.getExtension()}"
        """
        grep \
            -v "Failed Sample" ${sampleReport} \
            > ${outputFile}
        """
}

process convertSampleReportToFam {
    label 'smallMemory'

    tag "${sampleReport.getBaseName()}"

    input:
        path sampleReport
    output:
        path "${outputFile}"
    script:
        outputFile = "${sampleReport.getBaseName()}.fam"
        """
        cat ${sampleReport} \
            | cut \
                -f2,13-14 \
                -d',' \
                | awk \
                    'FS="," \
                    {print \$1,\$1,"0","0","-9","-9"}' \
                    | sed '1d' \
                        | sort \
                            | uniq -u
            > ${outputFile}
        """
}

process intersectFamFilesBySampleId {

    tag "${fam1}, ${fam2}"

    input:
        path fam1
        path fam2
    output:
        path "intersection.fam"
    script:
        """
        awk 'NR==FNR{a[\$2];next} \$2 in a' ${fam2} ${fam1} > intersection.fam
        """
}

process buildCohortData {
    label 'mediumMemory'
    label 'plink'

    tag "lgen, map, fam"

    cache 'lenient'

    input:
        path cohortLgen
        path cohortMap
        path illuminaFam
    output:
        path "nophenotypes.{bed,bim,fam}"
    script:
        """
        plink \
            --keep-allele-order \
            --lgen ${cohortLgen} \
            --map ${cohortMap} \
            --fam ${illuminaFam} \
            --no-parents \
            --no-sex \
            --no-pheno \
            --threads $task.cpus \
            --make-bed \
            --out nophenotypes
        """
}

process rebuildCohortDataWithPhenotypes {
    label 'plink'
    label 'mediumMemory'

    tag "cohortData, phenotypesFam"

    input:
        tuple path(cohortBed), path(cohortBim), path(unphenotypedFam)
        path phenotypesFam
    output:
        publishDir "${params.outputDir}/input/cohortData", mode: 'copy'
        path "${params.cohortName}.input.{bed,bim,fam}"
    script:
        """
        plink \
            --keep-allele-order \
            --bfile ${cohortBed.getBaseName()} \
            --keep ${phenotypesFam} \
            --pheno ${phenotypesFam} \
            --mpheno 4 \
            --update-sex ${phenotypesFam} 3 \
            --threads $task.cpus \
            --make-bed \
            --out ${params.cohortName}.input
        """
}

process rebuildCovariatesReport {
    label 'plink2'
    label 'mediumMemory'

    input:
        path covariatesReport
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
    output:
        publishDir "${params.outputDir}/input/cohortData", mode: 'copy'
        path "${params.cohortName}.input.cov"
    script:
        """
        plink2 \
            --bfile ${cohortBed.getBaseName()} \
            --covar ${covariatesReport} \
            --threads $task.cpus \
            --write-covar cols=fid \
            --out ${params.cohortName}.input
        """
}
