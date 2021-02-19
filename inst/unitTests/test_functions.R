test_K2Tworkflow <- function() {

    data(K2res)

    ## Pre-processing
    K2check <- K2preproc(K2eSet(K2res))
    checkTrue(is(K2check, "K2"), "K2preproc() did not return 'K2' object.\n")

    ## Taxonomies
    K2check <- K2tax(K2check)
    checkTrue(is(K2check, "K2"), "K2tax() did not return 'K2' object.\n")

    ## Dendrograms
    dend <- K2dendro(K2check)
    checkTrue(is(dend, "dendrogram"), "K2dendro() did not return 'dendrogram'
        object.\n")

    ## Meta-variable tests
    K2check <- runTestsMods(K2check,
                            infoClass=c(
                                sex="factor",
                                type="factor1",
                                score="numeric1"))
    checkTrue(is(K2check, "K2"), "runTestsMods() did not return 'K2' object.\n")

    ## Meta-variable results
    mvRes <- getTestsModTable(K2check)
    checkTrue(is(mvRes, "data.frame"), "getTestsModTable() did not return
        'data.frame'.\n")

    ## Differential Analysis
    K2check <- runDGEmods(K2check)
    checkTrue(is(K2check, "K2"), "runDGEmods() did not return 'K2' object.\n")

    ## DGE results
    gRes <- getDGETable(K2check)
    checkTrue(is(gRes, "data.frame"), "getDGETable() did not return
        'data.frame'.\n")

    ## Enrichment analysis
    K2check <- runGSEmods(K2res,
                          genesets=K2genesets(K2res))
    checkTrue(is(K2check, "K2"), "runGSEmods() did not return 'K2' object.\n")

    ## Genes 2 pathways
    gp <- getGenePathways(K2genesets(K2check))
    checkTrue(is(gp, "character"), "runGSEmods() did not return 'character'
        vector.\n")

    ## GSVA
    K2check <- runGSVAmods(K2check,
                         verbose=FALSE)
    checkTrue(is(K2check, "K2"), "runGSVAmods() did not return 'K2' object.\n")

    ## Aggregate GSVA scores
    aggList <- list(c("GS12", "GS1", "GS2"))
    K2check <- aggregateGSVAscores(aggList, K2check)
    checkTrue(is(K2check, "K2"), "aggregateGSVAscores() did not return 'K2'
        object.\n")

    # Differential pathway analysis
    K2check <- runDSSEmods(K2check)
    checkTrue(is(K2check, "K2"), "runDSSEmods() did not return 'K2' object.\n")

    ## Enrichment results
    eRes <- getEnrichmentTable(K2check)
    checkTrue(is(eRes, "data.frame"), "getEnrichmentTable() did not return
        'data.frame'.\n")

    ## Dashboard
    K2dashboard(K2res, output_dir=tempdir())
    K2resDir <- list.files(tempdir(), pattern="K2Taxonomer_",
        full.names=TRUE)[[1]]
    checkTrue(file.exists(file.path(K2resDir, "about.md")),
        "K2dashboard() didn't write files.")
    checkTrue(file.exists(file.path(K2resDir, "K2results.rds")),
        "K2dashboard() didn't write files.")
    checkTrue(file.exists(file.path(K2resDir, "K2Taxonomer.Rmd")),
        "K2dashboard() didn't write files.")

    ## Wrapper
    K2check <- runK2Taxonomer(K2eSet(K2res),
                  genesets=K2genesets(K2res),
                  infoClass=c(
                    sex="factor",
                    type="factor1",
                    score="numeric1")
                )
    checkTrue(is(K2check, "K2"),
        "runK2Taxonomer() did not return 'K2' object.\n")

}
