library("dplyr")
library("igraph")

# Load isoform-isoform interaction network
isoform_network_df <- read.csv("data-raw/results_PHYSICAL_Prob_Model_16_02_26.csv",
                               encoding = "utf-8")

# Create vertex table from both sides
vertices_1 <- isoform_network_df %>%
  select(enst = enst_1,
         gene = gene_1,
         protein = ensp_1)

vertices_2 <- isoform_network_df %>%
  select(enst = enst_2,
         gene = gene_2,
         protein = ensp_2)

vertices <- bind_rows(vertices_1,
                      vertices_2) %>%
  distinct(enst,
           .keep_all = TRUE)

# Create edge table
edges <- isoform_network_df %>%
  select(from = enst_1,
         to   = enst_2,
         pi,
         interact)

# Create final igraph object
isoform_network_graph <- graph_from_data_frame(
  d = edges,
  vertices = vertices,
  directed = FALSE
)

# Save isoform interaction network igraph object to internal data
usethis::use_data(isoform_network_graph,
                  internal = TRUE,
                  overwrite = TRUE)

# save raw .csv data
write.csv(isoform_network_df,
          file = "inst/extdata/iso_iso_network_raw.csv")
