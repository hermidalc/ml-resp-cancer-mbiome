suppressPackageStartupMessages(library("Biobase"))
set.seed(1982)

filterEset <- function(eset, features=NULL, samples=NULL) {
    if (!is.null(features) & !is.null(samples)) {
        return(eset[c(features),c(samples)])
    }
    else if (!is.null(features)) {
        return(eset[c(features),])
    }
    else if (!is.null(samples)) {
        return(eset[,c(samples)])
    }
    else {
        return(eset)
    }
}

esetClassLabels <- function(eset, samples=NULL) {
    if (!is.null(samples)) {
        return(eset$Class[c(samples)])
    }
    else {
        return(eset$Class)
    }
}

esetFeatureAnnot <- function(eset, annot, features=NULL) {
    if (!is.null(features)) {
        symbols <- as.character(featureData(eset)[c(features)][[annot]])
    }
    else {
        symbols <- as.character(featureData(eset)[[annot]])
    }
    symbols[is.na(symbols)] <- ""
    return(symbols)
}

limmaFeatureScore <- function(X, y) {
    suppressPackageStartupMessages(require("limma"))
    design <- model.matrix(~0 + factor(y))
    colnames(design) <- c("Class0", "Class1")
    fit <- lmFit(t(X), design)
    contrast.matrix <- makeContrasts(Class1VsClass0=Class1-Class0, levels=design)
    fit.contrasts <- contrasts.fit(fit, contrast.matrix)
    fit.b <- eBayes(fit.contrasts)
    results <- topTableF(fit.b, number=Inf, adjust.method="BH")
    results <- results[order(as.integer(row.names(results))),]
    return(list(results$F, results$adj.P.Val))
}

limmaFpkmFeatureScore <- function(X, y) {
    suppressPackageStartupMessages(require("limma"))
    design <- model.matrix(~0 + factor(y))
    colnames(design) <- c("Class0", "Class1")
    fit <- lmFit(t(log2(X + 1)), design)
    contrast.matrix <- makeContrasts(Class1VsClass0=Class1-Class0, levels=design)
    fit.contrasts <- contrasts.fit(fit, contrast.matrix)
    fit.b <- eBayes(fit.contrasts, trend=TRUE)
    results <- topTableF(fit.b, number=Inf, adjust.method="BH")
    results <- results[order(as.integer(row.names(results))),]
    return(list(results$F, results$adj.P.Val))
}

fcbfFeatureIdxs <- function(X, y, threshold=0) {
    results <- Biocomb::select.fast.filter(cbind(X, as.factor(y)), disc.method="MDL", threshold=threshold)
    results <- results[order(results$Information.Gain, decreasing=TRUE), , drop=FALSE]
    return(results$NumberFeature - 1)
}

cfsFeatureIdxs <- function(X, y) {
    X <- as.data.frame(X)
    colnames(X) <- seq(1, ncol(X))
    feature_idxs <- FSelector::cfs(as.formula("Class ~ ."), cbind(X, "Class"=as.factor(y)))
    return(as.integer(feature_idxs) - 1)
}

relieffFeatureScore <- function(X, y, num.neighbors=10, sample.size=5) {
    X <- as.data.frame(X)
    colnames(X) <- seq(1, ncol(X))
    results <- FSelector::relief(
        as.formula("Class ~ ."), cbind(X, "Class"=as.factor(y)),
        neighbours.count=num.neighbors, sample.size=sample.size
    )
    results <- results[order(as.integer(row.names(results))), , drop=FALSE]
    return(results$attr_importance)
}