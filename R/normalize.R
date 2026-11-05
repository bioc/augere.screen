#' @import augere.core
.normalize <- function(
    norm.control.barcodes, 
    norm.control.gene.field,
    norm.control.genes,
    norm.control.type.field,
    norm.control.types,
    norm.tmm
) {
    if (
        !is.null(norm.control.barcodes) || 
        (!is.null(norm.control.gene.field) && !is.null(norm.control.genes)) ||
        (!is.null(norm.control.type.field) && !is.null(norm.control.types))
    ) {
        template <- system.file("templates", "normalize-control.Rmd", package="augere.screen", mustWork=TRUE)
        parsed <- parseRmdTemplate(readLines(template))

        if (is.null(norm.control.barcodes)) {
            parsed[["barcodes"]] <- NULL
        } else {
            parsed[["barcodes"]] <- replacePlaceholders(parsed[["barcodes"]], list(NORM_BARCODES=deparseToString(norm.control.barcodes)))
        }

        if (is.null(norm.control.gene.field) || is.null(norm.control.genes)) {
            parsed[["genes"]] <- NULL
        } else {
            parsed[["genes"]] <- replacePlaceholders(
                parsed[["genes"]],
                list(
                    NORM_GENE_FIELD=deparseToString(norm.control.gene.field),
                    NORM_GENES=deparseToString(norm.control.genes)
                )
            )
        }

        if (is.null(norm.control.type.field) || is.null(norm.control.types)) {
            parsed[["types"]] <- NULL
        } else {
            parsed[["types"]] <- replacePlaceholders(
                parsed[["types"]],
                list(
                    NORM_TYPE_FIELD=deparseToString(norm.control.type.field),
                    NORM_TYPES=deparseToString(norm.control.types)
                )
            )
        }

    } else {
        template <- system.file("templates", "normalize-all.Rmd", package="augere.screen", mustWork=TRUE)
        parsed <- parseRmdTemplate(readLines(template))
    }

    if (norm.tmm) {
        parsed[["lib"]] <- NULL
    } else {
        parsed[["tmm"]] <- NULL
    }

    unlist(parsed, use.names=FALSE)
}
