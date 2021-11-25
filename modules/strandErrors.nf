include {
    checkCohortName;
    checkOutputDir;
    checkReferencePanelsDir;
    checkGeneticMapsDir;
    checkEmailAdressProvided;
    checkInputCohortData;
    userEmailAddressIsProvided;
    getBasicEmailSubject;
    getBasicEmailMessage;
    getCohortData;
    getCovariatesReport;
} from "${projectDir}/modules/base.nf"

def checkInputParams() {
    checkCohortName()
    checkOutputDir()
    checkEmailAdressProvided()
    checkInputCohortData('snvFiltered')
    //checkReferencePanelsDir()
    //checkGeneticMapsDir()
}

def getInputChannels() {
    return [
        getCohortData('ancestry'),
        getReferencePanels(),
        getGeneticMapsArchive(),
        getCovariatesReport('ancestry')]
}

process decompressGeneticMapsArchive {
    input:
        path geneticMapsArchive

    output:
        path "*.map"

    script:
        """
        unzip ${geneticMapsArchive}
        """
}

process selectAutosomalGenotypeSet {
    label 'plink2'
    label 'smallMemory'

    tag "cohortData"

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
    output:
        publishDir path: "${params.outputDir}/public-databases/genetic-maps", mode: 'copy'
        path "${cohortBed.getBaseName()}.vcf.gz"

    script:
        """
        plink2 \
            --bfile ${cohortBed.getBaseName()} \
            --chr 1-22 \
            --export vcf-4.2 bgz id-paste='fid' \
            --out ${cohortBed.getBaseName()}
        """
}

process indexReferencePanel {
    label 'htslib'
    label 'mediumMemory'

    tag "chr${chromosome}"

    input:
    tuple val(chromosome), path(vcfFile)
    output:
    tuple val(chromosome), path("${vcfFile}"), path("${vcfFile}.tbi")
    script:
        """
    tabix -p vcf ${vcfFile}
        """
}

process selectBiallelicSnvsWithBcftools {
    label 'bcftools'
    label 'mediumMemory'

    tag "referencePanel_chr${chromosome}"

    input:
        tuple val(chromosome), path(referencePanel), path(referencePanelIndex)
    output:
        publishDir path: "${params.outputDir}/public-databases/reference-panels", mode: 'copy'
        tuple val(chromosome), path("chr${chromosome}.refpanel.biallelic.vcf.gz")
    script:
    """
        bcftools \
            view \
            -M2 \
            -v snps \
            --threads $task.cpus \
            -Oz \
            -o chr${chromosome}.refpanel.biallelic.vcf.gz \
            ${referencePanel}
        """
}

process alignWithConformGt {
    label 'beagle'
    label 'mediumMemory'

    tag "genotypeSet, referencePanel_chr${chromosome}"

    input:
        tuple path(cohortGenotypes), val(chromosome), path(referencePanel)
    output:
        tuple \
            val(chromosome),
            path("${params.cohortName}.${chromosome}.vcf.gz")
        tuple \
            val(chromosome),
            path("${params.cohortName}.${chromosome}.log")
    script:
        """
        conform-gt \
            ref=${referencePanel} \
            gt=${cohortGenotypes} \
            chrom=${chromosome} \
            out=${params.cohortName}.${chromosome}
        """
}

process selectOppositeStrandSnvs {

    tag "${chromosome}"

    input:
        tuple \
            val(chromosome),
            path(alignmentLog)
    output:
        path "opposite-strand-snvs.${chromosome}.txt"
    script:
        """
        awk '{ if(\$10 ~ /OPPOSITE_STRAND/) {print \$3}}' ${alignmentLog} \
            > opposite-strand-snvs.${chromosome}.txt
        """
}

process concatenateFiles {

    input:
        path files
    output:
        path "opposite-strand_snvs.txt"
    script:
        """
        cat *.log > opposite-strand_snvs.txt
        """
}

process rebuildCohortData {
    label 'mediumMemory'
    label 'plink2'

    tag "haplotypeSet, inputFam"

    input:
        path cohortFam
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
    output:
        path "${params.cohortName}.strandFlipped.{bed,bim,fam}"
    script:
        """
        plink2 \
            --vcf ${cohortFam} \
            --fam ${cohortFam} \
            --threads $task.cpus \
            --make-bed \
            --double-id \
            --out ${params.cohortName}.aligned
        """
}

process flipGenotypesOnOppositeStrand {
    label 'mediumMemory'
    label 'plink2'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam), path(snvList)
    output:
        publishDir path: "${params.outputDir}/strandFlipped/cohortData", mode: 'copy'
        path("${params.cohortName}.strandFlipped.{bed,bim,fam}")
    script:
        """
        plink2 \
            --bfile ${cohortBed.getBaseName()} \
            --flip ${snvList} \
            --threads $task.cpus \
            --make-bed \
            --out ${params.cohortName}.strandFlipped
        """
}

process indexWithTabix {
    label 'htslib'
    label 'mediumMemory'

    tag "haplotypeSubset_chr${chromosome}"

    input:
        tuple val(chromosome), path(vcfFile)
    output:
        tuple path("${vcfFile}"), path("${vcfFile}.tbi")
    script:
        """
        tabix -p vcf ${vcfFile}
        """
}

process concatenateWithBcftools {
    label 'bcftools'
    label 'mediumMemory'

    tag "all haplotypeSubsets"

    input:
        path vcfFiles
    output:
        path "${params.cohortName}-phased.vcf.gz"
    script:
        """
        bcftools concat --threads $task.cpus -Oz -o ${params.cohortName}-phased.vcf.gz *.vcf.gz
        """
}

def getReferencePanels() {
    referencePanels = channel
        .fromPath(params.phase.referencePanels)
        .flatten()

    return indexByChromosome(selectAutosomes(referencePanels))
}

def getGeneticMapsArchive() {
    return channel.fromPath(params.phase.geneticMapsArchive)
}

def getGeneticMaps() {
    return channel
        .of(1..22)
        .map{ it -> [
            it,
            file(params.phase.geneticMapsDir + '*chr' + it + '.*.map')[0]]}
}

def selectAutosomes(inputPaths) {
    return inputPaths
        .filter {
            !(it.getName() =~ /chrX/) && !(it.getName() =~ /chrY/)
        }
}

def indexByChromosome(inputPaths) {
    return inputPaths
        .map {
            it ->
            [getChromosomeNumberFromString(it.getName()), it]
        }
}

def getChromosomeNumberFromString(inputString) {
    return (inputString =~ /chr(\d+)/)[0][1].toInteger()
}