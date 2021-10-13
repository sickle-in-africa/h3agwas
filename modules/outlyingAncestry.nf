include {
    getCohortData;
} from "${projectDir}/modules/base.nf"

def checkInputParams() {
    println 'checking...'
}
def getInputChannels() {
    return [
        getCohortData('phased'),
        getPrincipalComponentIds()]
}

def getPrincipalComponentIds() {
    return channel.of(1..params.ancestry.numberOfPrincipalComponentsToAnalyze)
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
        tuple path("${cohortGeno.getBaseName()}.evec"), path("${cohortGeno.getBaseName()}.eval"), path("${cohortGeno.getBaseName()}.log")
    script:
        """
        smartpca -p <(printf "genotypename: ${cohortGeno}
            snpname: ${cohortSnp}
            indivname: ${cohortInd}
            evecoutname: ${cohortGeno.getBaseName()}.evec
            evaloutname: ${cohortGeno.getBaseName()}.eval
            numoutlieriter: 0
            numoutlierevec: ${params.ancestry.numberOfPrincipalComponentsToAnalyze}
            numoutevec: ${params.ancestry.numberOfPrincipalComponentsToAnalyze}
            outliersigmathresh: 6
            altnormstyle: NO")  > ${cohortGeno.getBaseName()}.log
        sed \
            --in-place \
            --expression 's/^[ \t]*//' \
            ${cohortGeno.getBaseName()}.evec
        """
}

process testPrincipalComponentToPhenotypeAssociations {
    label 'tidyverse'

    tag "pc ${principalComponent}"

    input:
        tuple val(principalComponent), path(cohortEvec), path(cohortEval), path(cohortLog), path(cohortBed), path(cohortBim), path(cohortFam)
    output:
        path "${output}"
    script:
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

        write(paste(SELECTED_PRINCIPAL_COMPONENT, pvalue, rsquared, sep=","), file="${output}")
        """
}

process concatenateAcrossPrincipalComponents { 

    input:
        path pvalueAndVarianceTuples
    output:
        path "pvalueAndVarianceTable.csv"
    script:
        """
        echo "#PC,pvalue,r-squared" > pvalueAndVarianceTable.csv
        cat ${pvalueAndVarianceTuples} >> pvalueAndVarianceTable.csv
        sort \
            --key=1 \
            --numeric-sort \
            --field-separator="," \
            --output=pvalueAndVarianceTable.csv \
            pvalueAndVarianceTable.csv
        """
}

process calculateVarianceAddedByNewPrincipalComponent {
    input:
        path pvalueAndVarianceTable
    output:
        path "pvalueAndVarianceAndVarianceAddedTable.csv"
    script:
        template 'calculateVarianceAddedByNewPrincipalComponent.py'
}

process getNumberOfSignificantlyAssociatedPrincipalComponents {

    input:
        path pvalueTable
    output:
        env numberOfPrincipalComponents
    script:
        template 'getNumberOfSignificantlyAssociatedPrincipalComponents.sh'
}

process drawScreePlotForPrincipalComponents {
    label 'tidyverse'

    input:
        path pvalueAndVarianceAndVarianceAddedTable  
    output:
        path "scree-plot-for-principal-components.png"
    script:
        """
        #!/usr/bin/env Rscript --vanilla
        library(tidyverse)
        read.csv("${pvalueAndVarianceAndVarianceAddedTable}", header=TRUE) %>% 
            ggplot(aes(x=X.PC,y=r.squared_added)) +
                geom_bar(stat="identity") +
                xlab("principal component") +
                ylab("variance explained")
        ggsave("scree-plot-for-principal-components.png", dpi=300)
        """
}

process removeOutlyingSamples {
    label 'eigensoft'

    input:
        tuple path(cohortGeno), path(cohortSnp), path(cohortInd)
        val numberOfSignificantPrincipalComponents
    output:
        tuple path("${cohortGeno.getBaseName()}.outliers.evec"), path("${cohortGeno.getBaseName()}.outliers.eval"), path("${cohortGeno.getBaseName()}.outliers.log")
    script:
    """
        smartpca -p <(printf "genotypename: ${cohortGeno}
            snpname: ${cohortSnp}
            indivname: ${cohortInd}
            evecoutname: ${cohortGeno.getBaseName()}.outliers.evec
            evaloutname: ${cohortGeno.getBaseName()}.outliers.eval
            numoutlieriter: 0
            numoutlierevec: ${numberOfSignificantPrincipalComponents}
            numoutevec: ${params.ancestry.numberOfPrincipalComponentsToAnalyze}
            outliersigmathresh: 6
            altnormstyle: NO") > ${cohortGeno.getBaseName()}.outliers.log
        sed \
            --in-place \
            --expression 's/^[ \t]*//' \
            --expression 's/:/ /g' \
            ${cohortGeno.getBaseName()}.outliers.evec
    """
}

process drawPrincipalComponentPlotForSamples {
    label 'tidyverse'

    input:
        tuple path(cohortEvec), path(cohortEval), path(cohortLog)
    output:
        path "principal-component-plot-for-samples.pdf"
    script:
        """
        #!/usr/bin/env Rscript --vanilla
        library(tidyverse)

        pcx <- 1
        pcy <- 2

        PCAEVEC <-read.table("${cohortEvec}",head=T)
        colnames(PCAEVEC)[pcx+2] <- paste("PC",pcx,sep="")
        colnames(PCAEVEC)[pcy+2] <- paste("PC",pcy,sep="")
        colnames(PCAEVEC)[ncol(PCAEVEC)] <- "Pheno"
        pdf("principal-component-plot-for-samples.pdf")
        qplot(PCAEVEC[,pcx+2],PCAEVEC[,pcy+2], data=PCAEVEC, color=Pheno) + 
            xlab(paste("PC",(pcx),sep="")) + 
            ylab(paste("PC",(pcy),sep=""))
        dev.off()
        """
}

process drawPrincipalComponentPlotForOutliers {
    label 'tidyverse'

    input:
        tuple path(cohortEvec), path(cohortEval), path(cohortLog)
    output:
        path "principal-component-plot-for-outliers.pdf"
    script:
        """
        #!/usr/bin/env Rscript --vanilla
        library(tidyverse)

        pcx <- 1
        pcy <- 2

        PCAEVEC <-read.table("${cohortEvec}",head=T)
        colnames(PCAEVEC)[pcx+2] <- paste("PC",pcx,sep="")
        colnames(PCAEVEC)[pcy+2] <- paste("PC",pcy,sep="")
        colnames(PCAEVEC)[ncol(PCAEVEC)] <- "Pheno"
        pdf("principal-component-plot-for-outliers.pdf")
        qplot(PCAEVEC[,pcx+2],PCAEVEC[,pcy+2], data=PCAEVEC, color=Pheno) + 
            xlab(paste("PC",(pcx),sep="")) + 
            ylab(paste("PC",(pcy),sep=""))
        dev.off()
        """
}

process extractOutliers {

    input:
        tuple path(cohortEvec), path(cohortEval), path(cohortLog)
    output:
        path "${cohortEvec.getBaseName()}.outliers"
    script:
        """
        awk '/REMOVED/ {print \$3}' ${cohortLog} | sed 's/:/ /g' > ${cohortEvec.getBaseName()}.outliers
        """
}

process rebuildCohortData {
    label 'plink'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam), path(outlyingSamples)
    output:
        publishDir "${params.outputDir}/ancestry/cohortData", mode: 'copy'
        path "${params.cohortName}.ancestry.{bed,bim,fam}"
    script:
        """
        plink \
            --keep-allele-order \
            --bfile ${cohortBed.getBaseName()} \
            --remove ${outlyingSamples} \
            --make-bed \
            --threads $task.cpus \
            --out ${params.cohortName}.ancestry
        """
}