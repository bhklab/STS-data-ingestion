suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(PharmacoGx)
})

pset <- readRDS("extraction/data/raw/preclinical/CCLE_2019.rds")

cat("\nClass:\n")
print(class(pset))

cat("\nSlots:\n")
print(slotNames(pset))

cat("\nSample columns:\n")
if ("sample" %in% slotNames(pset)) {
  print(colnames(pset@sample))
  print(head(pset@sample))
}

cat("\nTreatment response names:\n")
if ("treatmentResponse" %in% slotNames(pset)) {
  print(names(pset@treatmentResponse))

  if ("profiles" %in% names(pset@treatmentResponse)) {
    cat("\nTreatment profile columns:\n")
    print(colnames(pset@treatmentResponse$profiles))
    print(head(pset@treatmentResponse$profiles))
  }
}

cat("\nMolecular profiles:\n")
if ("molecularProfiles" %in% slotNames(pset)) {
  print(names(pset@molecularProfiles))

  for (profile_name in names(pset@molecularProfiles)) {
    cat("\n---", profile_name, "---\n")
    profile <- pset@molecularProfiles[[profile_name]]

    print(class(profile))

    if (inherits(profile, "SummarizedExperiment")) {
      print(dim(profile))
      print(assayNames(profile))

      cat("First rownames:\n")
      print(head(rownames(profile)))

      cat("First colnames:\n")
      print(head(colnames(profile)))

      cat("rowData columns:\n")
      print(colnames(rowData(profile)))

      cat("colData columns:\n")
      print(colnames(colData(profile)))
    }
  }
}