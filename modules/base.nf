def checkOutputDir() {
    if (stringIsNull(params.outputDir)) {
        exit 1, 'params.outputDir not set -> please provide an output directory with enough space to save the results'
    }
}
def checkSelectedSnv() {
    if (stringIsNull(params.selectedSnv)) {
        exit 1, 'please provide a snv id!'
    }
}
def checkAssociationInput() {
    posibleInputs = ["input", "basicFiltered", "sampleFiltered", "snvFiltered", "phased", "ancestry"]
    if (stringIsNull(params.associationInput)) {
        exit 1, "params.associationInput not set -> please set it to one of the following items: ${posibleInputs}"
    }
    if (!posibleInputs.contains(params.associationInput)) {
        exit 1, "params.associationInput value not recognised -> please set it to one of the following items: ${posibleInputs}"
    }
}
def checkIlluminaGenotypeReports() {
    checkParamHoldingPathOrArrayOfPaths(
        params.input.genotypeReports,
        'params.input.genotypeReports not set -> please provide a genotype report path glob pattern')
}
def checkIlluminaSampleReport() {
    checkParamHoldingPathOrArrayOfPaths(
        params.input.sampleReport,
        'params.input.sampleReport not set -> please provide a sample report file path')
}
def checkIlluminaLocusReport() {
    checkParamHoldingPathOrArrayOfPaths(
        params.input.locusReport,
        'params.input.locusReport not set -> please provide a locus report file path')
}
def checkClinicalPhenotypeFam() {
    checkParamHoldingPathOrArrayOfPaths(
        params.input.clinicalPhenotypeFam,
        'params.input.clinicalPhenotypeFam not set -> please provide a clinical phenotype fam file path')
}
def checkCovariatesReport() {
    checkParamHoldingPathOrArrayOfPaths(
        params.input.covariatesReport,
        'params.input.covariatesReport not set -> please provide a covariates cov file path')
}
def checkEmailAdressProvided() {
    if (!userEmailAddressIsProvided()) {
        println 'You have not specified an email address; ' \
        	+ 'we will not email you the results of this workflow.'
    }
}
def checkCohortName () {
	if (stringIsNull(params.cohortName)) {
		exit 1, 'params.cohortName not set -> please provide a short cohort name to label your output files'
	}
}
def checkReferencePanelsDir() {
    if (stringIsNull(params.phase.referencePanelsDir)) {
        exit 1, 'params.phase.referencePanelsDir not set -> please provide a directory of reference panels e.g. from 1000 Genomes'
    }
    checkDirPath(params.phase.referencePanelsDir)
}
def checkGeneticMapsDir() {
    if (stringIsNull(params.phase.geneticMapsDir)) {
        exit 1, 'params.phase.geneticMapsDir not set -> please provide a directory of genetic maps'
    }
    checkDirPath(params.phase.geneticMapsDir)
}
def checkReferenceSequence() {
    checkParamHoldingPathOrArrayOfPaths(
        params.baseQC.referenceSequence,
        'params.baseQC.referenceSequence not set -> please provide a reference sequence fasta file path')
}
def checkInputCohortData(inputDataTag) {

    bedFile = file(params.outputDir + "${inputDataTag}/cohortData/" + params.cohortName + ".${inputDataTag}.bed")
    bimFile = file(params.outputDir + "${inputDataTag}/cohortData/" + params.cohortName + ".${inputDataTag}.bim")
    famFile = file(params.outputDir + "${inputDataTag}/cohortData/" + params.cohortName + ".${inputDataTag}.fam")

    if (!(bedFile.exists())) {
        exit 1, "could not find input bed file at ${bedFile} please check your output directory and try again"
    }
    if (!(bimFile.exists())) {
        exit 1, "could not find input bim file at ${bimFile} please check your output directory and try again"
    }
    if (!(famFile.exists())) {
        exit 1, "could not find input fam file at ${famFile} please check your output directory and try again"
    }
}
def checkParamHoldingPathOrArrayOfPaths(param, errorMessage) {
    if(param instanceof String) {
        if (stringIsNull(param)) {
            exit 1, errorMessage
        }
        checkFilePath(param)
    } else if(param instanceof Collection) {
        param.each {
            checkFilePath(it)
        }
    }
}
def checkFilePath(inputPath) {
    if (!stringIsGlobPattern(inputPath)) {
        checkFileOrDirExists(file(inputPath))
        checkFileNonempty(file(inputPath))
        checkFileIsFile(file(inputPath))
    }
}
def checkDirPath(inputPath) {
    checkDirPathStringFormat(inputPath)
    checkFileOrDirExists(file(inputPath))
    checkDirIsDir(file(inputPath))
}
def checkFileOrDirExists(inputFile) {
    if (!inputFile.exists()) {
        exit 1, "the file or directory:\n  ${inputFile}\ndoes not exist. Please check the path you entered is correct, and try again"
    }
}
def checkFileNonempty(inputFile) {
    if (inputFile.isEmpty()) {
        exit 1, "the file:\n  ${inputFile}\nis empty. Please check the file you specified, and try again"
    }
}
def checkFileIsFile(inputFile) {
    if (!inputFile.isFile()) {
        exit 1, "the file:\n  ${inputFile}\nis not actually a file. Maybe you entered a directory by mistake? Please check the file you specified, and try again"
    }
}
def checkDirPathStringFormat(inputPath) {
    if(inputPath[-1] != '/') {
        exit 1, "the directory path:\n  ${inputPath}\n is not formatted correctly. Please end all directory paths with a \'/\'"
    }
}
def checkDirIsDir(inputFile) {
    if (!inputFile.isDirectory()) {
        exit 1, "the directory:\n  ${inputFile}\n is not actually a directory. Perhaps you forgot to unzip an archive?"
    }
}
def userEmailAddressIsProvided() {
	return !(stringIsNull(params.email))
}
def stringIsNull(string) {
	return ( string =~ /NULL/ )
}
def stringIsGlobPattern(string) {
    return( string =~ /\*/ || string =~ /\?/ || string =~ /\[/ || string =~ /\{/ )
}
def getBasicEmailSubject() {
    return "[nextflow|h3agwaws] run ${workflow.runName} has finished"
}
def getBasicEmailMessage() {
    return """\
        Hi there, 

        Your nextflow job ${workflow.scriptName}: ${workflow.runName} has finished.
        Please check the attachments to this email,
        and the execution summary below. 

        All the best,
        H 3 A G W A S



        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()
}

def getCohortData(inputDataTag) {

    bed = channel.fromPath(
        params.outputDir + "${inputDataTag}/cohortData/" + params.cohortName + ".${inputDataTag}.bed")
    bim = channel.fromPath(
        params.outputDir + "${inputDataTag}/cohortData/" + params.cohortName + ".${inputDataTag}.bim")
    fam = channel.fromPath(
        params.outputDir + "${inputDataTag}/cohortData/" + params.cohortName + ".${inputDataTag}.fam")

    return bed.combine(bim).combine(fam)
}

def getCovariatesReport(inputDataTag) {

    return channel.fromPath(
        params.outputDir + "${inputDataTag}/cohortData/" + params.cohortName + ".${inputDataTag}.cov")
}

process rebuildCovariatesReport {
    label 'plink2'
    label 'mediumMemory'

    input:
        val label
        path covariatesReport
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
    output:
        publishDir "${params.outputDir}/${label}/cohortData", mode: 'copy'
        path "${params.cohortName}.${label}.cov"
    script:
        """
        plink2 \
            --bfile ${cohortBed.getBaseName()} \
            --covar ${covariatesReport} \
            --threads $task.cpus \
            --write-covar cols=fid \
            --out ${params.cohortName}.${label}
        """
}

def printWorkflowExitMessage() {
    if (workflow.success) {
        log.info "Workflow completed without errors".center(60)
    } else {
        log.error "Oops .. something went wrong!".center(60)
    }
    log.info "Check output files in folder:".center(60)
    log.info "${params.outputDir}".center(60)
}

process collectPlotsTogetherAndZip {

    input:
        val label
        path plots

    output:
        publishDir "${params.outputDir}plotArchives/", mode: 'copy'
        path "${label}.tar.gz"

    script:
        """
        mkdir plots
        cp -L *.p?? plots/
        tar -czf ${label}.tar.gz plots
        """
}

def sendWorkflowExitEmail() {
    if (userEmailAddressIsProvided()) {
        sendMail(
            to: "${params.email}",
            subject: getBasicEmailSubject(),
            body: getBasicEmailMessage())
   }
}

def sendWorkflowExitEmailWithPlots(analysisStep) {
    if (userEmailAddressIsProvided()) {
      sendMail(
          to: "${params.email}",
          subject: getBasicEmailSubject(),
          body: getBasicEmailMessage(),
      attach: "${params.outputDir}plotArchives/${analysisStep}.tar.gz")
  }
}