###############################################################################
# GOUtilities.R                                                               #
# Author:  Dario Ghersi                                                       #
# Version: 20161210                                                           #
# Goal:    utility functions for Gene Ontology                                #
###############################################################################

require(igraph)

buildGOGraphs <- function(goFile, regulatoryLinks=FALSE) {
  ## build and return a list of 'igraph' objects, one per namespace
  ## if the 'regulatoryLinks' is on, the edge type 'regulates'
  ## will be included in the graph
  ## N.B. If regulatory links are included, edges that go from one
  ## namespace to another will be present (e.g. some
  ## "biological_process" term will regulate some other
  ## "molecular_function" term)!!

  ## create the edge "lists" (actually, three matrices)
  edgeList <- list("biological_process"=c(), "molecular_function"=c(),
                   "cellular_component"=c())
  
  ## parse the gene ontology file and fill the edge lists
  goCon <- file(goFile, "r")
  goLine <- readLines(goCon, 1)
  children <- c()
  count <- 1
  while (length(goLine) > 0) {
    # reset the children
    if (goLine == "[Term]") {
      children <- c()
    }

    # term id
    if (length(grep(pattern="^id: ", x=goLine))) {
      children <- c(children, strsplit(goLine, split="id: ")[[1]][2])
    }

    # alt id
    if (length(grep(pattern="^alt_id: ", x=goLine))) {
      children <- c(children, strsplit(goLine, split="alt_id: ")[[1]][2])
    }

    # name space
    if (length(grep(pattern="^namespace: ", x=goLine))) {
      namespace <- strsplit(goLine, split="namespace: ")[[1]][2]
    }

    # "is a" relationship
    if (length(grep(pattern="^is_a: ", x=goLine))) {
      parent <- strsplit(goLine, split=" ")[[1]][2]
      for (child in children) {
        edgeList[[namespace]] <- rbind(edgeList[[namespace]],
                                       c(child, parent))
      }
    }

    # "part of" relationship
    if (length(grep(pattern="^relationship: part_of", x=goLine))) {
      parent <- strsplit(goLine,
                         split="relationship: part_of ")[[1]][2]
      parent <- strsplit(parent, " !")[[1]][1]
      for (child in children) {
        edgeList[[namespace]] <- rbind(edgeList[[namespace]],
                                       c(child, parent))
      }
    }

    # regulatory links
    if (regulatoryLinks) {
      if (length(grep(pattern="^relationship: regulates", x=goLine))) {
        parent <- strsplit(goLine,
                           split="relationship: regulates ")[[1]][2]
        parent <- strsplit(parent, " !")[[1]][1]
        for (child in children) {
          edgeList[[namespace]] <- rbind(edgeList[[namespace]],
                                         c(child, parent))
        }
      }
      if (length(grep(pattern="^relationship: positively_regulates",
                      x=goLine))) {
        parent <- strsplit(goLine,
                           split="relationship: positively_regulates ")[[1]][2]
        parent <- strsplit(parent, " !")[[1]][1]
        for (child in children) {
          edgeList[[namespace]] <- rbind(edgeList[[namespace]],
                                         c(child, parent))
        }
      }
      if (length(grep(pattern="^relationship: negatively_regulates",
                      x=goLine))) {
        parent <- strsplit(goLine,
                           split="relationship: negatively_regulates ")[[1]][2]
        parent <- strsplit(parent, " !")[[1]][1]
        for (child in children) {
          edgeList[[namespace]] <- rbind(edgeList[[namespace]],
                                         c(child, parent))
        }
      }
    }
    goLine <- readLines(goCon, 1)
  }
  close(goCon)

  ## create the igraph object
  gBP <- graph.edgelist(edgeList[["biological_process"]])
  gMF <- graph.edgelist(edgeList[["molecular_function"]])
  gCC <- graph.edgelist(edgeList[["cellular_component"]])
  
  return(list("biological_process"=gBP, "molecular_function"=gMF,
              "cellular_component"=gCC))
}

################################################################################

computeFrequencies <- function(annotations, GOGraph) {
  ## build a table with the frequency of each term in the annotations
  ## following the hierarchy
  ## the only terms to look at are the ones that have been used in the subset,
  ## their parents and their children
  ## N.B. the function works on the assumption that the 'annotation'
  ## list is "term-centric"

  ## get the terms and their parents
  nodes <- get.vertex.attribute(GOGraph, "name")
  annTerms <- names(annotations)
  annTermsPos <- unlist(sapply(1:length(annTerms), function(x)
                               which(nodes == annTerms[x])))
  
  parents <- nodes[unlist(neighborhood(GOGraph,
                                       nodes=annTermsPos, order=1E6,
                                       mode="out"))]
  terms <- unique(parents)
  numTerms <- length(terms)

  ## set up the frequency vector
  freq <- rep(0, numTerms)
  names(freq) <- terms

  ## compute the frequency for each term
  for (term in terms) {
    ## get the children terms
    termPos <- which(nodes == term)
    children <-
      unique(nodes[unlist(neighborhood(GOGraph, nodes=termPos,
                                       order=1E6, mode="in"))])
    freq[term] <- length(unique(unlist(annotations[children])))
  }

  return(freq)
}

###############################################################################

computeFunctHomogeneity <- function(terms, GOGraph, IC) {
  ## compute an index of functional homogeneity among terms
  ## as the pairwise functional similarity between them

  if (length(terms) == 1) # if only one term return 1
    return(1)

  ## compute the pairwise semantic similarity between terms
  pairwiseSemSim <- c()
  numTerms <- length(terms)
  for (i in 1:(numTerms - 1)) {
    for (j in (i + 1):numTerms) {
      pairwiseSemSim <- c(pairwiseSemSim,
                          semanticSim(terms[i], terms[j], GOGraph,
                                      IC))
    }
  }
  homogeneity <- mean(pairwiseSemSim)
  
  return(homogeneity)
}

###############################################################################

computeICTerms <- function(annotations, GOGraph) {
  ## compute the information content of each term

  ## get the term-centric annotations
  termCentric <- getTermCentricAnn(annotations)

  ## compute the frequency of each term
  freqs <- computeFrequencies(termCentric, GOGraph)

  ## compute the information content for each term
  numProteins <- length(annotations)
  IC <- -sapply(freqs, function(x) log10(x / numProteins))
  names(IC) <- names(freqs)
  
  return(IC)
}

###############################################################################

createAnnList <- function(annFile) {
  ## create an annotations list by parsing a file with the first
  ## column containing the ID of the gene and the rest the annotations

  ## read the data
  data <- scan(annFile, what="c", sep="\n")

  ## parse the data
  numAnn <- length(data)
  annList <- vector("list", numAnn)
  geneNames <- c()
  for (i in 1:numAnn) {
    fields <- strsplit(data[i], split="\t")[[1]]
    geneNames <- c(geneNames, fields[1])
    annList[[i]] <- fields[2:length(fields)]
  }
  names(annList) <- geneNames

  return(annList)
}

###############################################################################

createEnrichmentTable <- function(subset, termCentricAnn,
                                  enrichment, GOGraph, GODict,
                                  pvalue, outfile) {
  ## create a text file read for parsing with the enriched terms below the
  ## given p-value

  sortedTissueEnrichment <- sort(enrichment)
  termsPvalues <-sortedTissueEnrichment[sortedTissueEnrichment < pvalue]
  genesContrib <- c()
  
  ## add the genes contributing to each enriched term
  for (term in names(termsPvalues)) {
    withTerm <- intersect(haveTerm(termCentricAnn, term, GOGraph),
                          subset)
    withTerm <- paste(withTerm, collapse=",")
    genesContrib <- c(genesContrib, withTerm)
  }
  tableEnrichment <- cbind(names(termsPvalues), GODict[names(termsPvalues)],
                           signif(termsPvalues, 3), genesContrib)
  
  write.table(tableEnrichment, file=outfile, row.names=FALSE,
              col.names=FALSE, quote=FALSE, sep="\t")
}

###############################################################################
# Function by Ryan Hagenson
# Purpose: Provide an alternative createEnrichmentTable function that 
#          returns a data.frame rather than writes to file.
#          Almost entirely identical to original createEnrichmentTable function
###############################################################################

createEnrichmentFrame <- function(subset, termCentricAnn,
                                  enrichment, GOGraph, GODict,
                                  pvalue) {
  ## create a text file read for parsing with the enriched terms below the
  ## given p-value
  
  sortedTissueEnrichment <- sort(enrichment)
  termsPvalues <-sortedTissueEnrichment[sortedTissueEnrichment < pvalue]
  genesContrib <- c()
  
  ## add the genes contributing to each enriched term
  for (term in names(termsPvalues)) {
    withTerm <- intersect(haveTerm(termCentricAnn, term, GOGraph),
                          subset)
    withTerm <- paste(withTerm, collapse=",")
    genesContrib <- c(genesContrib, withTerm)
  }
  tableEnrichment <- cbind(names(termsPvalues), GODict[names(termsPvalues)],
                           signif(termsPvalues, 3), genesContrib)
  
  return(tableEnrichment)
}
###############################################################################
# End function by Ryan Hagenson
###############################################################################

###############################################################################

createGODict <- function(goFile) {
  ## create a dictionary that assignes the definition of a GO term
  ## to its ID

  ## parse the gene ontology file
  goCon <- file(goFile, "r")
  goLine <- readLines(goCon, 1)
  ids <- c()
  GODict <- c()
  while (length(goLine) > 0) {
    if (length(grep(pattern="^id: ", x=goLine))) {
      id <- strsplit(goLine, split="id: ")[[1]][2]
      ids <- c(ids, id)
    }
    if (length(grep(pattern="^name: ", x=goLine))) {
      termName <- strsplit(goLine, split="name: ")[[1]][2]
      GODict <- c(GODict, termName)
    }
    goLine <- readLines(goCon, 1)
  }
  names(GODict) <- ids
  return(GODict)
}

###############################################################################

enrichmentAnalysis <- function(subset, allProteins, annotations,
                               GOGraph, fdr=TRUE,
                               underrepresented=FALSE) {
  ## carry out a GO enrichment analysis for the proteins in the
  ## subset with an optional multiple-testing FDR correction
  ## if the 'underrepresented' flag is on then the terms that are actually
  ## statistically underrepresented are returned

  ## blacklist of terms to ignore (e.g. the root)
  blacklist <- c("GO:0008150", "GO:0005575", "GO:0003674")
  
  ## get the node names in the ontology graph
  GONodes <- get.vertex.attribute(GOGraph, "name")
  numNodes <- length(GONodes)

  ## get all the "term-centric" annotations (removing the
  ## annotations that don't belong to the network, if any)
  annotations <- annotations[intersect(names(annotations),
                                       allProteins)]
  termCentricAllAnn <- getTermCentricAnn(annotations)

  ## compute the background distribution
  freqBackground <- computeFrequencies(termCentricAllAnn, GOGraph)
  numBackground <- length(allProteins)

  ## get the term-centric annotations for the subset
  annSubset <- annotations[subset]
  termCentricAnn <- getTermCentricAnn(annSubset)
    
  ## get the parents of the terms
  subSize <- length(annSubset)
  terms <- names(termCentricAnn)
  numTerms <- length(terms)
  parents <- unlist(sapply(terms, function(x)
                           getAncestors(x, GOGraph)))
  parents <- unique(parents)
  parents <- setdiff(parents, blacklist) # remove root terms
  numParents <- length(parents)
    
  ## calculate the term frequency in the subset
  freq <- computeFrequencies(termCentricAnn, GOGraph)

  ## perform the hypergeometric test on each term
  enrichment <- c()
  for (term in parents) {
    notTermFreq <- numBackground - freqBackground[term]
    p <- phyper(freq[term] - 1, freqBackground[term], notTermFreq,
                subSize, lower.tail=underrepresented)
    enrichment <- c(enrichment, p)
  }
  names(enrichment) <- parents

  if (fdr) {
    enrichment <- p.adjust(enrichment, method="fdr")
  }
  enrichment <- sort(enrichment)
  
  return(enrichment)
}

###############################################################################

filterEnrichmentList <- function(enrichment, GOGraph, threshold) {
  ## filter and enrichment list by extracting the most specific
  ## terms that have
  ## a p-value below threshold

  ## extract all the terms below threshold
  terms <- names(enrichment)[enrichment < threshold]

  ## remove parent-child relationships
  specificTerms <- returnSpecificTerms(terms, GOGraph)

  return(enrichment[specificTerms])
}

###############################################################################

funSim <- function(setA, setB, GOGraph, IC, reference=0) {
  ## (Ref: the approach implemented is an extended version of what
  ## Agarwal et al used in their paper in PLoS Comp Bio (6) 2010)
  ## compute the functional similarity between two sets of annotations
  ## the variable 'reference' controls the denominator of the index:
  ## 0. the sum(IC) of the union of the terms
  ## 1. the IC of the terms in A
  ## 2. the IC of the terms in B

  ## get the parents of setA
  parentsA <- unlist(sapply(setA, function(x)
                            getAncestors(x, GOGraph)))
  termsA <- unique(parentsA)

  ## get the parents of setB
  parentsB <- unlist(sapply(setB, function(x)
                            getAncestors(x, GOGraph)))
  termsB <- unique(parentsB)

  ## compute the functional similarity
  intersectionAB <- intersect(termsA, termsB)
  numerator <- sum(IC[intersectionAB])
  if (reference == 0) {
    unionAB <- union(termsA, termsB)
    denominator <- sum(IC[unionAB])
  }
  else if (reference == 1) {
    denominator <- sum(IC[termsA])
  }
  else if (reference == 2) {
    denominator <- sum(IC[termsB])
  }
  sim <- numerator / denominator

  return(sim)
}

###############################################################################

getAncestors <- function(term, GOGraph) {
  ## get the parents of a term (the term itself is included)

  return(getTermsAboveOrBelow(term, GOGraph, "above"))
}

###############################################################################

getChildren <- function(term, GOGraph) {
  ## get the children of a term (the term itself is included)

  return(getTermsAboveOrBelow(term, GOGraph, "below"))
}

###############################################################################

getInformativeTerms <- function(termCentricAnn, cutoff, GOGraph) {
  ## return a reduced list of terms by removing very specific and
  ## very general terms, as done in Huang et al, Bioinformatics 2007 (23)
  ## The idea is to select terms that annotate more than cutoff proteins but
  ## whose children annotate less than cutoff proteins

  ## compute the frequencies for all the terms
  freq <- computeFrequencies(termCentricAnn, GOGraph)

  ## get the terms that satisfy the requirements
  nodes <- get.vertex.attribute(GOGraph, "name")
  termsWithFreq <- names(freq)
  candidateTerms <- names(freq[freq > cutoff])
  selectedTerms <- c()
  for (term in candidateTerms) {
    ## get the children
    termPos <- which(nodes == term)
    children <- nodes[unlist(neighborhood(GOGraph, nodes=termPos,
                                          order=1, mode="in"))]
    children <- intersect(setdiff(children, term), termsWithFreq)
    if (sum(freq[children] < cutoff) == length(children))
        selectedTerms <- c(selectedTerms, term)
  }

  return(selectedTerms)
}

###############################################################################

getTermCentricAnn <- function(annotations) {
  ## transform a "protein-centric" annotation list into a
  ## "term-centric" one

  ## build the vector with the terms
  numAnnotations <- length(annotations)
  protIds <- names(annotations)
  terms <- unique(unlist(annotations))
  termCentricAnn <-vector("list", length=length(terms))
  names(termCentricAnn) <- terms
  
  ## scan each term
  for (item in terms) {
    mask <- sapply(1:numAnnotations, function(x) item %in%
                   annotations[[x]])
    termCentricAnn[[item]] <- protIds[mask]
  }

  return(termCentricAnn)
}

###############################################################################

getTermsAboveOrBelow <- function(term, GOGraph, whichWay) {
  ## get the ancestors of a term (the term itself is included)

  ## get the nodes of the GO graph
  nodes <- get.vertex.attribute(GOGraph, "name")

  ## get the position of the term in the graph
  termPos <- which(nodes == term)

  ## get the terms
  if (whichWay == "above") {
    direction <- "out"
  }
  else if (whichWay == "below") {
    direction <- "in"
  }
  else {
    return(NA)
  }
  terms <- nodes[unlist(neighborhood(GOGraph, nodes=termPos,
                                     order=1E6, mode=direction))]
  terms <- unique(terms)

  return(terms)
}

###############################################################################

haveTerm <- function(termCentricAnn, term, GOGraph) {
  ## return a vector with the proteins that have the term (or a more
  ## specific one)

  ## get the children of the term (including the term)
  nodes <- get.vertex.attribute(GOGraph, "name")
  termPos <- which(nodes == term)
  children <- nodes[unlist(neighborhood(GOGraph, nodes=termPos,
                                        order=1E6, mode="in"))]

  ## get the proteins with the term
  proteinsWithTerm <- c()
  for (term in children) {
    proteinsWithTerm <- c(proteinsWithTerm, termCentricAnn[[term]])
  }
  proteinsWithTerm <- unique(proteinsWithTerm)

  return(proteinsWithTerm)
}

###############################################################################

map2GOSlim <- function(annotations, GOGraph, GOSlimGraph) {
  ##  map the terms in the annotation list to the slim subset

  ## blacklist of terms to ignore (e.g. the root)
  blacklist <- c("GO:0008150", "GO:0005575", "GO:0003674")
  
  ## create the slim list
  proteins <- names(annotations)
  numProteins <- length(proteins)
  slimAnn <- vector("list", length=numProteins)
  names(slimAnn) <- proteins

  ## populate the slim list
  slimNodes <- get.vertex.attribute(GOSlimGraph, "name")
  toRemove <- c()
  for (i in 1:numProteins) {
    terms <- annotations[[proteins[i]]]
    ## get the ancestors
    parents <- unlist(sapply(terms, function(x)
                             getAncestors(x, GOGraph)))
    parents <- unique(parents)
    slimTerms <- intersect(parents, slimNodes)
    slimTerms <- setdiff(slimTerms, blacklist)
    
    if (length(slimTerms) == 0)
      toRemove <- c(toRemove, i)
    else
      slimAnn[[proteins[i]]] <- slimTerms
  }

  ## get rid of the proteins without a slim ann
  slimAnn <- slimAnn[-toRemove]
  
  return(slimAnn)
}

###############################################################################

returnNumberOfClusters <- function(annotations, GOGraph, IC,
                                   type="Lin", threshold) {
  ## return the number of term clusters per protein at the given
  ## threshold
  ## N.B. the clustering algorithm is single-linkage agglomerative

  ## set up the list
  proteins <- names(annotations)
  numProteins <- length(proteins)
  numClusters <- vector("list", numProteins)
  names(numClusters) <- proteins

  ## carry out the clustering
  for (protein in proteins) {
    if (length(annotations[[protein]]) == 1) {
      numClusters[[protein]] <- 1
    }
    else {
      ## compute the semantic similarity
      terms <- annotations[[protein]]
      numTerms <- length(terms)
      M <- matrix(1, nrow=numTerms, ncol=numTerms)
      for (i in 1:(numTerms - 1)) {
        for (j in (i + 1):numTerms) {
          ss <- semanticSim(terms[i], terms[j], GOGraph, IC, type)
          M[i, j] <- M[j, i] <- ss
        }
      }
      ## perform the clustering
      M <- 1 - M # transform M in a dissimilarity matrix
      cl <- hclust(as.dist(M), method="single")
      ## cut the tree
      numClusters[[protein]] <- max(cutree(cl, h=threshold))
    }
  }

  return(numClusters)
}

###############################################################################

returnSpecificTerms <- function(terms, GOGraph) {
  ## return only the most specific terms by removing parent-child
  ## relationships

  ## get the list of the terms ancestors
  numTerms <- length(terms)
  ancestorList <- vector("list", length=numTerms)
  names(ancestorList) <- terms
  for (term in terms) {
    ancestorList[[term]] <- getAncestors(term, GOGraph)
    }
  
   ## extract the most specific terms
    specificTerms <- terms
    for (i in 1:(numTerms - 1)) {
      for (j in (i + 1):numTerms) {
        if (terms[i] %in% ancestorList[[terms[j]]]) {
          specificTerms <- setdiff(specificTerms, terms[i])
        }
        else if(terms[j] %in% ancestorList[[terms[i]]]) {
          specificTerms <- setdiff(specificTerms, terms[j])
        }
      }
    }

  return(specificTerms)
}

###############################################################################

semanticSim <- function(term1, term2, GOGraph, IC, type="Lin") {
  ## compute the semantic similarity between two terms

  if (!type %in% c("Lin", "Resnik")) {
    print("Type not supported")
    return(NA)
  }

  ## deal with identical terms
  if (term1 == term2) {
    if (type == "Lin") {
      return(1.0)
    }
    if (type == "Resnik") {
      return(IC[term1])
    }
  }

  ## get the parents
  parents1 <- getAncestors(term1, GOGraph)
  parents2 <- getAncestors(term2, GOGraph)

  ## get the common ancestors
  commonAncestors <- intersect(parents1, parents2)

  ## get the IC of the Most Informative Common Ancestor
  ICMICA <- max(IC[commonAncestors])
  names(ICMICA) <- c()

  if (type == "Resnik") {
    return(ICMICA)
  }
  else if (type == "Lin") {
    lin <-  2 * ICMICA / (IC[term1] + IC[term2])
    names(lin) <- c()
    return(lin)
  }
}


###############################################################################

