# library(testthat); library(augere.screen); source("test-runScreen.R")

library(SummarizedExperiment)
cd <- DataFrame(group = rep(LETTERS[1:4], each=5), age=1:20, batch=factor(rep(1:5, 4)))
means <- 2^runif(1000, 1, 10)
means[1] <- 0 # to test the filtering.
mat <- matrix(rnbinom(length(means) * nrow(cd), mu=means, size=10), nrow=length(means), ncol=nrow(cd))
se <- SummarizedExperiment(list(counts=mat), colData=cd)
rowData(se)$type <- sample(c("NTC", "NEG", "other"), length(means), replace=TRUE)
rowData(se)$gene <- paste0("GENE-", sample(LETTERS, length(means), replace=TRUE))

test_that("runScreen works with defaults, no genes", {
    output <- tempfile()
    pdf(file=NULL)
    res <- runScreen(se, groups="group", comparisons=c("A", "B"), output.dir=output)
    dev.off()

    expect_identical(length(res$gene), 0L) 
    expect_identical(length(res$barcode), 1L) 
    expect_identical(names(res$barcode), "Increase in `A` over `B`")
    expect_identical(rownames(res$barcode[[1]]), rownames(se))
    expect_identical(rownames(se), rownames(res$normalized))

    expect_identical(rownames(res$barcode[[1]]), rownames(se))
    expect_type(res$barcode[[1]]$PValue, "double")
    expect_type(res$barcode[[1]]$LogFC, "double")
    expect_type(res$barcode[[1]]$AveAb, "double")
    expect_type(res$barcode[[1]]$FDR, "double")

    roundtrip <- augere.core::readResult(file.path(output, "results", "barcode-1"))
    expect_identical(roundtrip$x, res$barcode[[1]])
    expect_match(roundtrip$metadata$title, "Increase")
    expect_identical(roundtrip$metadata$differential_barcode_abundance$contrast$left, list("A"))
    expect_identical(roundtrip$metadata$differential_barcode_abundance$contrast$right, list("B"))
    expect_match(roundtrip$metadata$differential_barcode_abundance$normalization, "TMM")

    roundtrip <- augere.core::readResult(file.path(output, "results", "normalized"))
    expect_identical(dimnames(roundtrip$x), dimnames(res$normalized))
    expect_match(roundtrip$metadata$title, "Normalized")
})

test_that("runScreen works with default genes", {
    output <- tempfile()
    res <- runScreen(se, groups="group", comparisons=c("A", "B"), gene.field="gene", output.dir=output, suppress.plots=TRUE)
    expect_identical(names(res$gene), c("simes", "fry", "holm-min"))

    expected.names <- sort(unique(rowData(se)$gene))
    for (meth in names(res$gene)) {
        expect_match(names(res$gene[[meth]]), "Increase in")
        df <- res$gene[[meth]][[1]]
        expect_identical(rownames(df), expected.names)

        expect_type(df$NumBarcodes, "integer")
        expect_type(df$NumUp, "integer")
        expect_type(df$NumDown, "integer")
        expect_type(df$PValue, "double")
        expect_type(df$LogFC, "double")
        expect_type(df$FDR, "double")
        expect_true(all(df$Direction %in% c("up", "down", NA)))

        roundtrip <- augere.core::readResult(file.path(output, "results", paste0(meth, "-1")))
        expect_identical(roundtrip$metadata$differential_gene_abundance$method, meth)
        expect_match(roundtrip$metadata$differential_gene_abundance$normalization, "TMM")
    }
})

test_that("runScreen works with log-fold changes", {
    tmp <- tempfile()
    out <- runScreen(se, group="group", comparisons=c("C", "A"), lfc.threshold=0.5, output.dir=tmp, suppress.plots=TRUE, save.results=FALSE)

    fname <- readLines(file.path(tmp, "report.Rmd"))
    expect_true(any(grepl("treat\\(.*lfc=0.5", fname)))
})

test_that("runScreen works with custom subsetting", {
    tmp <- tempfile()
    out <- runScreen(se, group="group", comparisons=c("C", "A"), subset.factor="batch", subset.levels=c("1","3","5"), gene.field="gene", output.dir=tmp, suppress.plots=TRUE)

    # Subsetting works correctly.
    expect_lt(ncol(out$normalized), ncol(se))
    expect_identical(unique(out$normalized$batch), factor(c("1", "3", "5"), 1:5))

    meta <- augere.core::readResult(file.path(tmp, "results", "barcode-1"))$metadata
    expect_identical(meta$differential_barcode_abundance$subset, "batch IN (1,3,5)")
    meta <- augere.core::readResult(file.path(tmp, "results", "normalized"))$metadata
    expect_identical(meta$summarized_experiment$subset, "batch IN (1,3,5)")

    for (meth in names(out$gene)) {
        roundtrip <- augere.core::readResult(file.path(tmp, "results", paste0(meth, "-1")))
        expect_identical(roundtrip$metadata$differential_gene_abundance$subset, "batch IN (1,3,5)")
    }
})

test_that("runScreen works with some other filtering parameters", {
    tmp <- tempfile()
    out <- runScreen(se, group="group", comparisons=c("C", "A"), filter.reference.factor="group", filter.reference.levels=c("C", "D"), output.dir=tmp, suppress.plots=TRUE, save.results=FALSE)
    expect_false(rowData(out$normalized)$retained[1])

    out <- runScreen(se, group="group", comparisons=c("C", "A"), filter.default=FALSE, output.dir=tmp, suppress.plots=TRUE, save.results=FALSE)
    expect_true(rowData(out$normalized)$retained[1])
})

test_that("runScreen works with control-based normalization", {
    # Don't suppress plots, we want to check that the normalization-related plots for each strategy are generated correctly.
    tmp <- tempfile()
    pdf(file=NULL)
    out <- runScreen(se, group="group", comparisons=c("C", "A"), norm.control.type.field="type", norm.control.types=c("NTC", "NEG"), output.dir=tmp)
    dev.off()

    meta <- augere.core::readResult(file.path(tmp, "results", "barcode-1"))$metadata
    expect_identical(meta$differential_barcode_abundance$normalization, "TMM on controls (NTC, NEG)")

    meta <- augere.core::readResult(file.path(tmp, "results", "normalized"))$metadata
    expect_identical(meta$summarized_experiment$normalization, "TMM on controls (NTC, NEG)")

    for (meth in names(out$gene)) {
        roundtrip <- augere.core::readResult(file.path(output, "results", paste0(meth, "-1")))
        expect_identical(roundtrip$metadata$differential_gene_abundance$normalization, "TMM on controls (NTC, NEG)")
    }

    tmp <- tempfile()
    pdf(file=NULL)
    out <- runScreen(se, group="group", comparisons=c("C", "A"), norm.control.type.field="type", norm.control.types=c("NTC", "NEG"), norm.tmm=FALSE, output.dir=tmp)
    dev.off()

    meta <- augere.core::readResult(file.path(tmp, "results", "barcode-1"))$metadata
    expect_identical(meta$differential_barcode_abundance$normalization, "library size on controls (NTC, NEG)")
})

test_that("runScreen works with multiple comparisons", {
    # Don't suppress plotting, we want to check that multiple plots are correctly generated.
    tmp <- tempfile()
    pdf(file=NULL)
    out <- runScreen(se, group="group", comparisons=list(foo=c("C", "A"), bar=c("B", "D")), gene.field="gene", output.dir=tmp)
    dev.off()

    expect_identical(names(out$barcode), c("foo", "bar"))
    expect_match(augere.core::readResult(file.path(tmp, "results", "barcode-1"))$meta$title, "foo")
    expect_match(augere.core::readResult(file.path(tmp, "results", "barcode-2"))$meta$title, "bar")
    expect_match(augere.core::readResult(file.path(tmp, "results", "fry-1"))$meta$title, "foo")
    expect_match(augere.core::readResult(file.path(tmp, "results", "simes-2"))$meta$title, "bar")

    # Comparing to the single-comparison results.
    tmp1 <- tempfile()
    alone1 <- runScreen(se, group="group", comparisons=c("C", "A"), gene.field="gene", subset.group=FALSE, output.dir=tmp1, save.results=FALSE, suppress.plots=TRUE)
    expect_identical(out$barcode$foo, alone1$barcode[[1]])
    expect_identical(out$gene$simes$foo, alone1$gene$simes[[1]])

    tmp2 <- tempfile()
    alone2 <- runScreen(se, group="group", comparisons=c("B", "D"), gene.field="gene", subset.group=FALSE, output.dir=tmp2, save.results=FALSE, suppress.plots=TRUE)
    expect_identical(out$barcode$bar, alone2$barcode[[1]])
    expect_identical(out$gene$fry$bar, alone2$gene$fry[[1]])
})

test_that("runScreen works with custom design and contrasts", {
    tmp <- tempfile()
    custom <- runScreen(se, design=~0 + group, contrasts="groupB - groupA", gene.field="gene", output.dir=tmp, suppress.plots=TRUE)
    meta <- augere.de::readResult(file.path(tmp, "results", "barcode-1"))$metadata
    expect_identical(meta$differential_barcode_abundance$contrast$type, "custom")
    meta <- augere.de::readResult(file.path(tmp, "results", "simes-1"))$metadata
    expect_identical(meta$differential_gene_abundance$contrast$type, "custom")

    # Making sure we get the same results with a simple design.
    tmp <- tempfile()
    comp <- runScreen(se, group="group", comparisons=c("B", "A"), gene.field="gene", subset.group=FALSE, output.dir=tmp, suppress.plots=TRUE, save.results=FALSE)
    expect_identical(custom$normalized, comp$normalized)
    expect_identical(unname(custom$barcode), unname(comp$barcode))
    expect_identical(unname(custom$gene$fry), unname(comp$gene$fry))
    expect_identical(unname(custom$gene$simes), unname(comp$gene$simes))
})

test_that("runScreen works with duplicateCorrelation", {
    tmp <- tempfile()
    out <- runScreen(se, group="group", comparisons=c("C", "A"), dc.block="batch", gene.field="gene", output.dir=tmp, suppress.plots=TRUE, save.results=FALSE)

    fname <- readLines(file.path(tmp, "report.Rmd"))
    expect_true(any(grepl("duplicateCorrelation(", fname, fixed=TRUE)))
    expect_gte(sum(grepl("block=dc.block,", fname, fixed=TRUE)), 3) # remember, it's used in fry() as well. 
})

test_that("runScreen works with no consolidation methods", {
    # Attempt saving to make sure it skips the per-gene results.
    output <- tempfile()
    res <- runScreen(se, groups="group", comparisons=c("A", "B"), consolidation=NULL, gene.field="gene", output.dir=output, suppress.plots=TRUE)
    expect_identical(length(res$gene), 0L)

    lines <- readLines(file.path(output, "report.Rmd"))
    expect_false(any(grepl("by.gene", lines, fixed=TRUE)))
    expect_false(any(grepl("extra.gene.data", lines, fixed=TRUE)))
})

test_that("runScreen works with no summary", {
    output <- tempfile()
    res <- runScreen(se, groups="group", comparisons=c("A", "B"), report.summary=FALSE, gene.field="gene", output.dir=output, suppress.plots=TRUE, save.results=FALSE)
    expect_identical(names(res$gene), c("simes", "fry", "holm-min")) 

    expected.names <- sort(unique(rowData(se)$gene))
    for (meth in res$gene) {
        df <- meth[[1]]
        expect_identical(rownames(df), expected.names)
        expect_null(df$NumUp)
        expect_null(df$NumDown)
        expect_null(df$LogFC)
        expect_type(df$NumBarcodes, "integer")
        expect_type(df$PValue, "double")
        expect_type(df$FDR, "double")
    }
})

test_that("runScreen works with extra gene data", {
    output <- tempfile()
    res <- runScreen(se, groups="group", comparisons=c("A", "B"), gene.data="type", gene.field="gene", output.dir=output, suppress.plots=TRUE, save.results=FALSE)
    expect_identical(names(res$gene), c("simes", "fry", "holm-min")) 

    expected.names <- sort(unique(rowData(se)$gene))
    for (meth in res$gene) {
        df <- meth[[1]]
        expect_identical(rownames(df), expected.names)
        expect_type(df$NumUp, "integer")
        expect_type(df$NumDown, "integer")
        expect_type(df$LogFC, "double")
        expect_type(df$type, "character")
        expect_type(df$NumBarcodes, "integer")
        expect_type(df$PValue, "double")
        expect_type(df$FDR, "double")
    }
})

test_that("runScreen works with extra fry arguments", {
    output <- tempfile()
    res <- runScreen(se, groups="group", comparisons=c("A", "B"), fry.args=list(standardize="none"), consolidation="fry", gene.field="gene", output.dir=output, suppress.plots=TRUE, save.results=FALSE)
    expect_identical(names(res$gene), "fry")
    lines <- readLines(file.path(output, "report.Rmd"))
    expect_true(any(grepl("standardize=\"none\"", lines, fixed=TRUE)))
})

test_that("runScreen works with some of the other voom-related options", {
    tmp <- tempfile()
    out <- runScreen(se, group="group", comparisons=c("C", "A"), robust=FALSE, trend=TRUE, quality=FALSE, gene.field="gene", output.dir=tmp, suppress.plots=TRUE, save.results=FALSE)

    fname <- readLines(file.path(tmp, "report.Rmd"))
    expect_false(any(grepl("voomWithQualityWeights(", fname, fixed=TRUE)))
    expect_true(any(grepl("voom(", fname, fixed=TRUE)))
    expect_true(any(grepl("trend=TRUE", fname, fixed=TRUE)))
    expect_false(any(grepl("robust=TRUE", fname, fixed=TRUE)))
})

test_that("runScreen works with dry runs and no saving", {
    tmp <- tempfile()
    out <- runScreen(se, group="group", comparisons=c("C", "A"), gene.field="gene", dry.run=TRUE, output.dir=tmp)
    expect_null(out)
    expect_true(file.exists(file.path(tmp, "report.Rmd")))

    tmp <- tempfile()
    out <- runScreen(se, group="group", comparisons=c("C", "A"), gene.field="gene", output.dir=tmp, suppress.plots=TRUE, save.results=FALSE)
    expect_s4_class(out$normalized, "SummarizedExperiment")
    expect_s4_class(out$barcode[[1]], "DFrame")
    expect_s4_class(out$gene$fry[[1]], "DFrame")
    expect_false(file.exists(file.path(tmp, "results")))
})
