# querying relational tables for given taxa
library(here)
library(tidyverse)
library(stringr)
library(janitor)
library(shiny)
library(reactable)
library(htmltools)
library(shinythemes)
library(waiter)

# data preprocess
# read taxonomy
cavs <- read_csv("data/cav_taxonomy.csv")
# read specimen data
specimen_Data <- read_csv("data/master-specimens.csv") %>%
  mutate(specimenID = paste(museumAbbr, specimen)) %>%
  select(sp, specimenID, collection = standardized_name, country = components.country, study)
# read accession data
accession_Data <- read_csv("data/master-accessions.csv")
# read morphology data
morphology_Data <- read_csv("data/master-morph.csv")
# read ecology data
diet_Data <- read_csv("data/master-diet.csv")
habits_Data <- read_csv("data/master-habits.csv")
habitat_Data <- read_csv("data/master-habitats.csv")
# trait definitions
trdefs <- read_csv("data/master-trdefinitions.csv")
# study urls
linkTib <- read_csv("data/study_links.csv")

# focalSp <- "ctenomys pearsoni"
# data sources
master_datasets <- mget(ls(pattern = "Data$"))
# define fn to match scientific names
query_datasets <- function(focalSp) {
  taxDetails <- cavs %>%
    filter(Reduce(`&`, lapply(paste0("^", focalSp, "$"), agrepl, ScientificName, ignore.case = T, max.distance = c(insertions = 0.12))))
  if (nrow(taxDetails) == 0) {
    print("no match in MammalDiversity")
  }
  matches_sp <- map(master_datasets, filter, Reduce(
    `&`,
    lapply(paste0("^", focalSp, "$"), agrepl, sp, ignore.case = T, max.distance = c(insertions = 0.12))
  ))
  matches_sp_hits <- simplify(map(matches_sp, nrow))
  if (sum(matches_sp_hits) == 0) {
    message("no matches or incomplete search term")
  } 
  metadataTibsPre <-
    matches_sp %>%
    map_at(which(matches_sp_hits != 0), remove_empty, "cols")
  hits_w_data <- which(map_lgl(metadataTibsPre, ~ any(names(.x) == "sp")))
  metadataTibs <- map_at(metadataTibsPre, hits_w_data, ~ (rename(.x, scientific_name = sp))) %>% map(clean_names)
  resultsEnf <- map(metadataTibs, nrow) %>%
    simplify() %>%
    enframe() %>%
    rename(Data = 1, `Number of records` = 2) %>%
    mutate(Data = str_remove(Data, "_Data")) %>%
    mutate(Data = str_to_title(Data))
  titlesFixed <- c(
    "MD Taxonomy" = "MD Taxonomy",
    "Data Matches" = "Data Matches"
  )
  titlesDyn <- c(
    "accession_Data" = "Gene Sequences",
    "diet_Data" = "Diet Data",
    "habitat_Data" = "Habitat Data",
    "habits_Data" = "Habits Data",
    "morphology_Data" = "Morphology Data",
    "specimen_Data" = "Specimen Data",
    "trait_defs" = "Trait definitions"
  )
  captionsFixed <- c(
    "MD Taxonomy" = "Matches in MammalDiversity",
    "Data Matches" = "Number of records found for each data type"
  )
  captionsDyn <- c(
    "accession_Data" = "Gene sequence identifiers",
    "diet_Data" = "Feeding and diet habits",
    "habitat_Data" = "Habitat preferences",
    "habits_Data" = "Substrate use and locomotion mode",
    "morphology_Data" = "Morphology data as shared in each publication",
    "specimen_Data" = "Specimen identifiers and collection information",
    "trait_defs" = "Trait descriptions as defined in each study. Includes landmark configurations for geometric morphometric approaches."
  )


  getTRdefs <- function(listOut) {
    if (resultsEnf[5, 2] != 0) {
      studsDefs <- listOut$specimen_Data$study
      tib_tr <- trdefs %>%
        filter(study %in% studsDefs) %>%
        remove_empty("cols")
      append(metadataTibs, list(trait_defs = tib_tr))
    } else {
      listOut
    }
  }

  metadataTibs_tr <- getTRdefs(metadataTibs)
  n_trdefs <- tibble(Data = "trait_defs", `Number of records` = NA)
  if (resultsEnf[5, 2] != 0) {
    n_trdefs[, 2] <- length(metadataTibs_tr$specimen_Data$study)
  } else {
    n_trdefs[, 2] <- 0
  }
  # <-
  matches_sp_hits_new <- c(matches_sp_hits, deframe(n_trdefs))
  metadataTibsF <- metadataTibs_tr[names(matches_sp_hits_new[matches_sp_hits_new != 0])]

  titlesFinal <- c(titlesFixed, titlesDyn[names(matches_sp_hits_new[matches_sp_hits_new != 0])])
  captionsFinal <- c(captionsFixed, captionsDyn[names(matches_sp_hits_new[matches_sp_hits_new != 0])])
  return(
    list(simplify(append(list(taxonomy = taxDetails, matches = resultsEnf), metadataTibsF)), titlesFinal, captionsFinal)
  )
}

query_datasets("ctenomys pearsoni")
ui <- fluidPage(theme=shinytheme("paper"),
  use_butler(),
  br(),
  fluidRow(
  column(11,
       titlePanel("Patterns in research and data-sharing for the study of form and function in caviomorph rodents:  
                  Interactive Web App", windowTitle="Caviomorph Ecomorphology App")
     ),
   column(1,
              actionButton("github",
                          label = "View Code",
                          width = "95px",
                          onclick ="window.open(`https://github.com/luisDVA/Caviomorph-Ecomorphology-Resources-App`, '_blank')",
                          style="color: #fff; background-color: #4c9ed9; border-color: black"),
          style="position:absolute;right:4em"
   )
  ),
  h3("Luis D. Verde Arregoitia, Pablo Teta & Guillermo D'Elía"),
  br(),
  br(),
  p("Enter a taxonomic name below to search.",style= "font-size: 22px"),
  p("Output is split into interactive tables for each data type, which can be sorted by clicking on the column names, or exported separately using the 'Download' button. Entries in the 'study' column provide live links to each publication.",style= "font-size: 22px"),
  fluidRow(column(4, align="center",offset=4,
       textInput("textstring", h4("Query Term"), NULL),
       tags$style("#textstring {font-size:23px;font-style: italic;color: black}")
        )),
  br(),
  p("To see the source code for this application and for offline use, click on the 'View Code' button at the top of this page.",style= "font-size: 20px"),
  br(),
  p("Note: This app uses approximate string matching. The spelling of scientific names and level of resolution may vary between data sources (e.g. searching for \"Kannabateomys\" will produce more matches than searching for \"Kannabateomys amblyonyx\"). Allow for latency as the results update.",style= "font-size: 20px"),
  mainPanel(
    # UI output
    uiOutput("dtss")
  )
)




server <- function(input, output) {
  observeEvent(input$textstring,{
    show_butler()
  output$dtss <- renderUI({
      if (input$textstring != "") {
      data <- query_datasets(input$textstring)
      titles <- data[[2]]
      captions <- data[[3]]

      tibble(
        tbl = data[[1]],
        tle = titles,
        cpts = captions
      ) %>%
        pmap(function(tbl, tle, cpts) {
          tags$div(
            tags$hr(),
            tags$h1(tle),
            fluidRow(
              column(
                8,
                tags$p(cpts)
              ),
              column(
                3,
                downloadHandler(
                  filename = function() {
                    paste0(gsub(" ", "-", tle), "-", gsub(" ", "-", input$textstring), ".csv")
                  },
                  content = function(file) {
                    write.csv(tbl, file)
                  }
                )
              )
            ),
            tags$br(),
            reactable(tbl,
              if ("study" %in% names(tbl)) {
                columns <- list(
                  study = colDef(cell = function(value, index) {
                    # Render as a link
                    suppressMessages(urlpre <- left_join(tbl, linkTib) %>% select(doiurl))
                    url <- urlpre$doiurl[index]
                    htmltools::tags$a(href = url, target = "_blank", as.character(value))
                  })
                )
              },
              defaultColDef = colDef(
                headerStyle = list(background = "#0E4D2D", color = "white")
              ),
              highlight = TRUE, compact = TRUE, sortable = TRUE
            ),
            tags$br()
          )
        })
    }
  })
  hide_butler()
  })
}

shinyApp(ui, server)
