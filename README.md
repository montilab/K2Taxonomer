# K2Taxonomer
This is an R package under development to perform iterative K=2 clustering, as well as annotate and visualize results.

```r
library(K2Taxonomer)

# Run K2Taxonomer
K2summary <- runK2Taxonomer(eSet,
                           cohorts = "group",
                           covariates = NULL,
                           use = "Z",
                           nFeats = nrow(eSet)/20,
                           nBoots = 200,
                           clustFunc = hclust_wrapper,
                           info = info,
                           genesets = genesets,
                           qthresh = 0.05,
                           oneoff = TRUE,
                           geneURL = NULL,
                           genesetURL = NULL)

# Generate K2 dashboard
K2dashboard(K2summary,
            analysis_name = "example",
            output_dir = "results",
            compile = FALSE)
```