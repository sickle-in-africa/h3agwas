include {
    checkOutputDir;
    checkCohortName;
    checkEmailAdressProvided;
    userEmailAddressIsProvided;
    getCohortData;
    getCovariatesReport;
    getBasicEmailSubject;
    getBasicEmailMessage;
    checkReferenceSequence;
    checkInputCohortData;
} from "${projectDir}/modules/base.nf"

def checkInputParams() {
    checkCohortName()
    checkOutputDir()
    checkEmailAdressProvided()
    checkInputCohortData('input')
    checkReferenceSequence()
}

def getInputChannels() {
    return [
        getCohortData('input'),
        getReferenceSequence(),
        getCovariatesReport('input')]
}

def getReferenceSequence() {
    return channel
        .fromPath(params.baseQC.referenceSequence)
}

process selectDuplicatedVariants {
    label 'plink'

    tag "cohortData"

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)

    output:
        path 'plink.dupvar'

    script:
        """
        plink \
            --keep-allele-order \
            --bfile ${cohortBed.getBaseName()} \
            --list-duplicate-vars \
            ids-only \
            suppress-first
        """
}

process removeDuplicatedVariants {
    label 'plink'

    tag "cohortData, duplicatedVariants"

    input:
    tuple path(cohortBed), path(cohortBim), path(cohortFam)
    path duplicatedVariants

    output:
        path("duplicatedVariantsRemoved.{bed,bim,fam}")

    script:
        """
        plink \
            --keep-allele-order \
            --bfile ${cohortBed.getBaseName()} \
            -exclude ${duplicatedVariants} \
            --make-bed \
            --out duplicatedVariantsRemoved
        """
}


process alignGenotypesToReference {
    label 'mediumMemory'
    label 'plink2'

    tag "filteredCohortData, reference"

    cache 'lenient'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
        path referenceSequence
    output:
        path "temporary.vcf.gz"
    script:
        plinkBase = cohortBed.getBaseName()
        """
        plink2 \
            --bfile ${plinkBase} \
            --fa ${referenceSequence} \
            --ref-from-fa force \
            --normalize \
            --threads $task.cpus \
            --export vcf-4.2 id-paste=iid bgz \
            --real-ref-alleles \
            --out temporary
        """
}

process selectBiallelicSnvs {
    label 'mediumMemory'
    label 'bcftools'

    tag "alignedGenotypes"

    input:
        path genotypeSet
    output:
        path "filtered.vcf.gz"
    script:
        """
        bcftools \
            view \
            -m2 \
            -M2 \
            -v snps \
            --threads $task.cpus \
            -Oz \
            -o filtered.vcf.gz \
        ${genotypeSet}
        """
}

process rebuildCohortData {
    label 'mediumMemory'
    label 'plink2'

    tag "alignedGenotypes, filteredCohortFam"

    input:
        path alignedGenotypes
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
    output:
        path "aligned.{bed,bim,fam}"
    script:
        """
        plink2 \
            --vcf ${alignedGenotypes} \
            --fam ${cohortFam} \
            --threads $task.cpus \
            --make-bed \
            --double-id \
            --out aligned
        """
}

process selectCommonSnvs {
    label 'plink'
    label 'smallMemory'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
    output:
        path "common-snvs.{bed,bim,fam}"
    script:
        """
        plink \
            --keep-allele-order \
            --bfile ${cohortBed.getBaseName()} \
            --maf ${params.baseQC.minAlleleFrequency} \
            --make-bed \
            --out common-snvs
        """
}

process removeSamplesWithLowCallRates {
    label 'plink'
    label 'smallMemory'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
    output:
        path "high-sample-call-rate.{bed,bim,fam}"
    script:
        """
        plink \
            --keep-allele-order \
            --bfile ${cohortBed.getBaseName()}  \
            --mind ${params.baseQC.maxMissingnessPerSample} \
            --make-bed \
            --out high-sample-call-rate
        """
}

process removeSnvsWithLowCallRates {
    label 'plink'
    label 'smallMemory'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
    output:
        publishDir "${params.outputDir}/basicFiltered/cohortData", mode: 'copy'
        path "${params.cohortName}.basicFiltered.{bed,bim,fam}"
    script:
        """
        plink \
            --keep-allele-order \
            --bfile ${cohortBed.getBaseName()}  \
            --geno ${params.baseQC.maxMissingnessPerSnv} \
            --make-bed \
            --out ${params.cohortName}.basicFiltered
        """
}
