# Load the zen4R package
library(zen4R)

# Initialize ZenodoManager with your API token
zenodo <- ZenodoManager$new(token = "PLRku93CDvT7DtmYddtcpnKNfextYMi7Sp5LQBls7wA7p1VuF1Z6nGmyyEjB")

# Step 1: Create an empty record
record <- zenodo$createEmptyRecord(reserveDOI = TRUE)
record_id <- record$id

# Step 2: Edit the record metadata
record$metadata$title <- "Large Dataset for My R Package"
record$metadata$description <- "This is a large dataset used in my R package. Hosted on Zenodo due to size limitations."
record$metadata$creators <- list(list(name = "MengYan Li", affiliation = "Independent Developer"))

# Step 3: Update the record on Zenodo
zenodo$depositRecord(record = record, publish = FALSE)

# Step 4: Upload the dataset file
zenodo$uploadFile(
  path = "C:/Users/User1/Desktop/R_Package_For_Tianxi/MUGS/data/pairs.rel.CV.rda",
  record = record
)

# Step 5: Publish the record
zenodo$publishRecord(record$id)

# Step 6: Display the record details
cat("Zenodo Record URL:", record$links$record_html, "\n")
cat("Zenodo DOI:", record$metadata$doi, "\n")
