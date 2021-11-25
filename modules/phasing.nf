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
        getCohortData('strandFlipped'),
        getReferencePanels(),
	    getGeneticMapsArchive(),
        getCovariatesReport('strandFlipped')]
}

process selectGenotypeSetWithPlink {
    label 'plink2'
    label 'smallMemory'

    tag "cohortData"

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)

    output:
        path "${cohortBed.getBaseName()}.vcf.gz"

    script:
        """
        plink2 \
            --bfile ${cohortBed.getBaseName()} \
            --export vcf-4.2 bgz id-paste='fid' \
            --out ${cohortBed.getBaseName()}
        """
}

process recodeMaleXHaploidAsDiploid {
    label 'smallMemory'
    
    tag 'genotypeSet'
    
    input:
        path cohortGenotypes
    output:
        path "recodedCohortGenotypes.vcf.gz"
    script:
        """
        # ** Note to future developers! **
        # the sed command appears twice on each line (i.e. 4 times
        # in total) on purpose. It is to catch the cases were the 
        # input strings overlap, for example:
        #   ...1/0  1   1   1/1...
        
        zcat ${cohortGenotypes} \
            | sed 's|\t0\t|\t0/0\t|g' | sed 's|\t0\t|\t0/0\t|g' \
            | sed 's|\t1\t|\t1/1\t|g' | sed 's|\t1\t|\t1/1\t|g' \
            > recodedCohortGenotypes.vcf
        gzip recodedCohortGenotypes.vcf
        """
}

process alignWithConformGt {
    label 'beagle'
    label 'mediumMemory'

    tag "genotypeSet, referencePanel_chr${chromosome}"

    input:
        tuple path(cohortGenotypes), val(chromosome), path(referencePanel)
    output:
        tuple val(chromosome), path("${params.cohortName}.${chromosome}.vcf.gz")
    script:
        """
        conform-gt \
            ref=${referencePanel} \
            gt=${cohortGenotypes} \
            chrom=${chromosome} \
            out=${params.cohortName}.${chromosome}
        """
}

process phaseWithBeagle {
    label 'beagle'
    label 'mediumMemory'

    tag "genotypeSubset_chr${chromosome}, referencePanel_chr${chromosome}, geneticMap_chr${chromosome}"

    input:
        tuple val(chromosome), path(alignedGenotypes), path(referencePanel), path(geneticMap)
    output:
        tuple val(chromosome), path("${params.cohortName}.${chromosomeString}.phased.vcf.gz")
    script:
        chromosomeString = ( chromosome < 10 ) ? "0${chromosome}" : "${chromosome}"
        """
        beagle \
            ref=${referencePanel} \
            impute=true \
            map=${geneticMap} \
            gt=${alignedGenotypes} \
            chrom=${chromosome} \
            nthreads=${task.cpus} \
            window=${params.phase.windowSize} \
            out=${params.cohortName}.${chromosomeString}.phased
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

process rebuildCohortDataWithPlink() {
    label 'mediumMemory'
    label 'plink2'

    tag "haplotypeSet, inputFam"

    input:
        path haplotypes
        tuple path(unphasedBed), path(unphasedBim), path(unphasedFam)
    output:
        publishDir path: "${params.outputDir}/phased/cohortData", mode: 'copy'
        path "${params.cohortName}.phased.{bed,bim,fam}"
    script:
        """
        plink2 \
            --vcf ${haplotypes} \
            --fam ${unphasedFam} \
            --threads $task.cpus \
            --make-bed \
            --double-id \
            --out ${params.cohortName}.phased
        """
}

def getReferencePanels() {
    referencePanels = channel
        .fromPath("${params.outputDir}/public-databases/reference-panels/chr*.refpanel.biallelic.vcf.gz")
        .flatten()

    return indexByChromosome(selectAutosomes(referencePanels))
}

def getGeneticMaps() {
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

