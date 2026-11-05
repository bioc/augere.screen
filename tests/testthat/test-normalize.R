# library(testthat); library(augere.screen); source("test-normalize.R")

library(SummarizedExperiment)
rd <- DataFrame(gene = sample(LETTERS, 1000, replace=TRUE), type = rep(c("NTC", "NEG", "other"), c(50, 50, 900)))
means <- 2^rexp(1000)
means[1] <- 0 # to test the filtering.
mat <- matrix(rnbinom(length(means) * 10, mu=means, size=10), nrow=length(means), ncol=10)
rownames(mat) <- sprintf("BARCODE-%s", seq_len(nrow(mat)))
se <- SummarizedExperiment(list(counts=mat), rowData=rd)

library(edgeR)
y <- DGEList(mat)
y$genes <- data.frame(origin=seq_len(nrow(mat)))
y <- y[seq_len(nrow(y)) %% 2 == 0,] # keeping every even barcode, to check for correct handling of non-trivial filters.

ctrl.barcodes <- list(bctrl=sprintf("BARCODE-%s", seq_len(50)))
ctrl.genes <- list(gctrl=LETTERS[1:3])
ctrl.types <- c("NTC", "NEG")

test_that("control-based normalization works as expected with only barcodes", {
    cmds <- augere.screen:::.normalize(
        norm.control.barcodes=ctrl.barcodes,
        norm.control.gene.field=NULL,
        norm.control.type.field=NULL,
        norm.tmm=TRUE
    )
    expect_true(any(grepl("norm.barcodes <-", cmds)))
    expect_false(any(grepl("norm.genes <-", cmds)))
    expect_false(any(grepl("norm.types <-", cmds)))
    expect_true(any(grepl("normLibSizes", cmds)))

    env <- new.env()
    env$se <- se
    env$y <- y
    tmp <- tempfile()
    dir.create(tmp)
    augere.core::compileReport(file=file.path(tmp, "report.Rmd"), contents=cmds, env=env)

    expect_type(env$y$samples$norm.factors, "double")
    expect_false(all(env$y$samples$norm.factors==1))
    expect_identical(env$use.for.norm, rownames(y) %in% unlist(ctrl.barcodes))
    expect_match(env$norm.meta, "TMM on controls (bctrl)", fixed=TRUE)
})

test_that("control-based normalization works as expected with only genes", {
    cmds <- augere.screen:::.normalize(
        norm.control.barcodes=NULL,
        norm.control.gene.field="gene",
        norm.control.genes=ctrl.genes,
        norm.control.type.field=NULL,
        norm.tmm=TRUE
    )
    expect_false(any(grepl("norm.barcodes <-", cmds)))
    expect_true(any(grepl("norm.genes <-", cmds)))
    expect_false(any(grepl("norm.types <-", cmds)))

    env <- new.env()
    env$se <- se
    env$y <- y
    tmp <- tempfile()
    dir.create(tmp)
    augere.core::compileReport(file=file.path(tmp, "report.Rmd"), contents=cmds, env=env)

    expect_type(env$y$samples$norm.factors, "double")
    expect_false(all(env$y$samples$norm.factors==1))
    expect_identical(env$use.for.norm, rd$gene[y$genes$origin] %in% unlist(ctrl.genes))
    expect_match(env$norm.meta, "TMM on controls (gctrl)", fixed=TRUE)
})

test_that("control-based normalization works as expected with only types", {
    cmds <- augere.screen:::.normalize(
        norm.control.barcodes=NULL,
        norm.control.gene.field=NULL,
        norm.control.type.field="type",
        norm.control.types=ctrl.types,
        norm.tmm=TRUE
    )
    expect_false(any(grepl("norm.barcodes <-", cmds)))
    expect_false(any(grepl("norm.genes <-", cmds)))
    expect_true(any(grepl("norm.types <-", cmds)))

    env <- new.env()
    env$se <- se
    env$y <- y
    tmp <- tempfile()
    dir.create(tmp)
    augere.core::compileReport(file=file.path(tmp, "report.Rmd"), contents=cmds, env=env)

    expect_type(env$y$samples$norm.factors, "double")
    expect_false(all(env$y$samples$norm.factors==1))
    expect_identical(env$use.for.norm, rd$type[y$genes$origin] %in% ctrl.types)
    expect_identical(env$norm.meta, "TMM on controls (NTC, NEG)")
})

test_that("control-based normalization works as expected with everything", {
    cmds <- augere.screen:::.normalize(
        norm.control.barcodes=ctrl.barcodes,
        norm.control.gene.field="gene",
        norm.control.genes=ctrl.genes,
        norm.control.type.field="type",
        norm.control.types=ctrl.types,
        norm.tmm=TRUE
    )
    expect_true(any(grepl("norm.barcodes <-", cmds)))
    expect_true(any(grepl("norm.genes <-", cmds)))
    expect_true(any(grepl("norm.types <-", cmds)))

    env <- new.env()
    env$se <- se
    env$y <- y
    tmp <- tempfile()
    dir.create(tmp)
    augere.core::compileReport(file=file.path(tmp, "report.Rmd"), contents=cmds, env=env)

    expect_type(env$y$samples$norm.factors, "double")
    expect_false(all(env$y$samples$norm.factors==1))
    expect_identical(env$use.for.norm, 
        rownames(y) %in% unlist(ctrl.barcodes) |
        rd$gene[y$genes$origin] %in% unlist(ctrl.genes) |
        rd$type[y$genes$origin] %in% ctrl.types
    )
    expect_identical(env$norm.meta, "TMM on controls (bctrl, gctrl, NTC, NEG)")
})

test_that("control-based normalization works with library sizes", {
    cmds <- augere.screen:::.normalize(
        "control",
        norm.control.barcodes=ctrl.barcodes,
        norm.control.gene.field=NULL,
        norm.control.type.field=NULL,
        norm.tmm=FALSE
    )
    expect_true(any(grepl("norm.barcodes <-", cmds)))
    expect_false(any(grepl("norm.genes <-", cmds)))
    expect_false(any(grepl("norm.types <-", cmds)))
    expect_false(any(grepl("normLibSizes", cmds)))

    env <- new.env()
    env$se <- se
    env$y <- y
    tmp <- tempfile()
    dir.create(tmp)
    augere.core::compileReport(file=file.path(tmp, "report.Rmd"), contents=cmds, env=env)

    expect_type(env$y$samples$norm.factors, "double")
    expect_false(all(env$y$samples$norm.factors==1))
    expect_identical(env$norm.meta, "library size on controls (bctrl)")
})

test_that("normalization on all barcodes works as expected", {
    {
        cmds <- augere.screen:::.normalize(
            norm.control.barcodes=NULL,
            norm.control.gene.field=NULL,
            norm.control.type.field=NULL,
            norm.tmm=TRUE
        )
        expect_false(any(grepl("norm.barcodes <-", cmds)))
        expect_false(any(grepl("norm.genes <-", cmds)))
        expect_false(any(grepl("norm.types <-", cmds)))
        expect_true(any(grepl("normLibSizes", cmds)))

        env <- new.env()
        env$se <- se
        env$y <- y
        tmp <- tempfile()
        dir.create(tmp)
        augere.core::compileReport(file=file.path(tmp, "report.Rmd"), contents=cmds, env=env)

        expect_type(env$y$samples$norm.factors, "double")
        expect_false(all(env$y$samples$norm.factors==1))
        expect_identical(env$norm.meta, "TMM on all barcodes")
    }

    {
        cmds <- augere.screen:::.normalize(
            norm.control.barcodes=NULL,
            norm.control.gene.field=NULL,
            norm.control.type.field=NULL,
            norm.tmm=FALSE
        )
        expect_false(any(grepl("norm.barcodes <-", cmds)))
        expect_false(any(grepl("norm.genes <-", cmds)))
        expect_false(any(grepl("norm.types <-", cmds)))
        expect_false(any(grepl("normLibSizes", cmds)))

        env <- new.env()
        env$se <- se
        env$y <- y
        tmp <- tempfile()
        dir.create(tmp)
        augere.core::compileReport(file=file.path(tmp, "report.Rmd"), contents=cmds, env=env)

        expect_type(env$y$samples$norm.factors, "double")
        expect_true(all(env$y$samples$norm.factors==1))
        expect_identical(env$norm.meta, "library size on all barcodes")
    }
})
