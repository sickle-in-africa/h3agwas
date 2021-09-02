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
        """
}
