suppressPackageStartupMessages({
  library(data.table)
  library(SummarizedExperiment)
  library(PharmacoGx)
})

rds_path <- "extraction/data/raw/preclinical/CCLE_2015.rds"

ccle_2015 <- readRDS(rds_path)

ccle_2015 <- updateObject(ccle_2015, verbose = TRUE)

cat("\nClass:\n")
print(class(ccle_2015))

cat("\nSlots:\n")
print(slotNames(ccle_2015))

cat("\nSample columns:\n")
if ("sample" %in% slotNames(ccle_2015)) {
  print(colnames(ccle_2015@sample))
  print(head(ccle_2015@sample))
}

cat("\nTreatment response names:\n")
if ("treatmentResponse" %in% slotNames(ccle_2015)) {
  print(names(ccle_2015@treatmentResponse))

  for (obj_name in names(ccle_2015@treatmentResponse)) {
    cat("\n--- treatmentResponse$", obj_name, " ---\n", sep = "")
    obj <- ccle_2015@treatmentResponse[[obj_name]]
    print(class(obj))

    if (is.data.frame(obj) || is.data.table(obj) || is.matrix(obj)) {
      print(dim(obj))
      print(colnames(obj))
      print(head(obj))
    } else {
      print(str(obj))
    }
  }
}

cat("\nTreatment slot:\n")
if ("treatment" %in% slotNames(ccle_2015)) {
  print(class(ccle_2015@treatment))
  print(colnames(ccle_2015@treatment))
  print(head(ccle_2015@treatment))
}

cat("\nMolecular profiles:\n")
if ("molecularProfiles" %in% slotNames(ccle_2015)) {
  print(names(ccle_2015@molecularProfiles))

  for (profile_name in names(ccle_2015@molecularProfiles)) {
    cat("\n--- ", profile_name, " ---\n", sep = "")

    profile <- ccle_2015@molecularProfiles[[profile_name]]

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

      cat("First rowData rows:\n")
      print(head(as.data.frame(rowData(profile))))

      cat("colData columns:\n")
      print(colnames(colData(profile)))

      cat("First colData rows:\n")
      print(head(as.data.frame(colData(profile))))
    } else {
      print(str(profile))
    }
  }
}