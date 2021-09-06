include {
    getCohortData;
} from "${projectDir}/modules/base.nf"

def checkInputParams() {
    println 'checking...'
}
def getInputChannels() {
    return getCohortData('phased')
}

process convertCohortDataToEigensoftFormat {
    label 'eigensoft'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
    output:
        tuple path("${cohortBed.getBaseName()}.eigenstratgeno"), path("${cohortBed.getBaseName()}.snp"), path("${cohortBed.getBaseName()}.ind")
    script:
        """
        convertf -p <(printf "genotypename: ${cohortBed}
            snpname: ${cohortBim}
            indivname: ${cohortFam}
            outputformat: EIGENSTRAT
            genotypeoutname: ${cohortBed.getBaseName()}.eigenstratgeno
            snpoutname: ${cohortBed.getBaseName()}.snp
            indivoutname: ${cohortBed.getBaseName()}.ind
            familynames: NO")
        """
}

process selectPrincipalComponents {
    label 'eigensoft'

    input:
        tuple path(cohortGeno), path(cohortSnp), path(cohortInd)
    output:
        tuple path("${cohortGeno.getBaseName()}.evec"), path("${cohortGeno.getBaseName()}.eval")
    script:
        """
        smartpca -p <(printf "genotypename: ${cohortGeno}
            snpname: ${cohortSnp}
            indivname: ${cohortInd}
            evecoutname: ${cohortGeno.getBaseName()}.evec
            evaloutname: ${cohortGeno.getBaseName()}.eval
            numoutlieriter: 0
            numoutlierevec: 100
            numoutevec: 100
            outliersigmathresh: 6
            altnormstyle: NO")
        sed \
            --in-place \
            --expression 's/^[ \t]*//' \
            ${cohortGeno.getBaseName()}.evec
        """
}

process testPrincipalComponentToPhenotypeAssociations {
    label 'tidyverse'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)        
        tuple path(cohortEvec), path(cohortEval)
    output:
        path "${output}"
    script:
        principalComponent = 5
        output = "pc${principalComponent}.csv"
        """
        #!/usr/bin/env Rscript --vanilla

        PHENOTYPE_COLUMN <- 6
        SELECTED_PRINCIPAL_COMPONENT <- ${principalComponent}

        phenotype <- as.data.frame(read.table("${cohortFam}", header=FALSE)[,c(1,PHENOTYPE_COLUMN)])
        colnames(phenotype) <- c('sample_id', 'phenotype')

        significantPrincipalComponents <- as.data.frame(read.table("${cohortEvec}", header=FALSE)[,c(1,2:(SELECTED_PRINCIPAL_COMPONENT+1))])
        colnames(significantPrincipalComponents) <- c('sample_id', paste('pc', 1:SELECTED_PRINCIPAL_COMPONENT, sep=""))

        associationData <- merge(phenotype, significantPrincipalComponents, by='sample_id')
        associationData\$sample_id <- NULL

        associationSummary <- summary(lm(phenotype ~ ., data=associationData))
        pvalue <- associationSummary\$coefficients[SELECTED_PRINCIPAL_COMPONENT+1,4]
        rsquared <- associationSummary\$r.squared

        cat(c(SELECTED_PRINCIPAL_COMPONENT, pvalue, rsquared), file="${output}", sep=",")
        """
}


