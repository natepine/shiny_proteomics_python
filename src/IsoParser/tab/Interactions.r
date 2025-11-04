# Tab constants
bioplexData <- readRDS(file = paste0(dataPath, "bioplex.rds"))
#ID lookup table
bioplexMap <- bind_rows(
      bioplexData %>% transmute(GeneSymbol = SymbolA, UniprotID = UniprotA, BioplexID = GeneA),
      bioplexData %>% transmute(GeneSymbol = SymbolB, UniprotID = UniprotB, BioplexID = GeneB)
   ) %>% distinct(GeneSymbol, .keep_all = TRUE)

output$ui_select_protein_bioplex <- proteinSelect$addInput("protein_select_bioplex", label = paste("Search for", dataset$idLabel))

output$plot_bioplex <- renderVisNetwork({
   nearest <- dataset$proteins$info[values$nearestIDs, "GeneSymbol"] %>%
      str_remove("__.*$")# Remove customId portion of GeneSymbol
   current <- nearest[[1]]

   df_gs <- select(bioplexData, SymbolA, SymbolB)

   # Filter out irrelevant edges:
   if(input$bioplex_graph_type == "current"){
      df_gs_current <- filter(df_gs, SymbolA == current | SymbolB == current)
      symbols <- unique(c(df_gs_current$SymbolA, df_gs_current$SymbolB))
      df_gs <- filter(df_gs, SymbolA %in% symbols & SymbolB %in% symbols)
      required_nodes <- current
   }else if(input$bioplex_graph_type == "nearest"){
      df_gs <- filter(df_gs, SymbolA %in% nearest & SymbolB %in% nearest)
      required_nodes <- nearest
   }else{
      df_gs <- filter(df_gs, SymbolA %in% nearest | SymbolB %in% nearest)
      # Set all yellows to 2 nearest at least:
      symbols <- unique(c(df_gs$SymbolA, df_gs$SymbolB))
      symbols <- sapply(symbols, function(sym){
         if(nrow(filter(df_gs, SymbolA == sym | SymbolB == sym)) >= 2){
            return(sym)
         }
         return(NA)
      })
      symbols <- symbols[!is.na(symbols)]
      df_gs <- filter(df_gs, SymbolA %in% symbols & SymbolB %in% symbols)
      if(input$bioplex_graph_type == "all"){
         required_nodes <- nearest
      }else{
         # Set all other nodes to at least n connections:
         symbols <- unique(c(df_gs$SymbolA, df_gs$SymbolB))
         symbols <- sapply(symbols, function(sym){
            if(nrow(filter(df_gs, SymbolA == sym | SymbolB == sym)) >= as.numeric(input$bioplex_graph_type)){
               return(sym)
            }
            return(NA)
         })
         symbols <- symbols[!is.na(symbols)]
         df_gs <- filter(df_gs, SymbolA %in% symbols | SymbolB %in% symbols)
         required_nodes <- nearest
      }
   }

   geneSymbol <- sort(unique(c(df_gs$SymbolA, df_gs$SymbolB, required_nodes)))
   nodes <- data.frame(id = geneSymbol, label = geneSymbol)
   ids <- bioplexMap[match(nodes$id, bioplexMap$GeneSymbol), ]
   nodes$UniprotID <- ids$UniprotID
   nodes$BioplexID <- ids$BioplexID
   nodes <- nodes %>% mutate(title = 
      paste0(
         bold(paste0("Gene Symbol: ", id, "<br>")),
         ifelse(is.na(BioplexID),
            color("Bioplex: unavailable<br>", "#E2E7EB"),
            paste0("<a href=\"https://bioplex.hms.harvard.edu/explorer/externalQuery.php?geneQuery=", BioplexID, "\"  target=\"_blank\" style=\"color:#FF9100\">BioplexID: ", BioplexID, "</a><br>")),
         "<a href=\"https://www.uniprot.org/uniprot/", UniprotID, "\"  target=\"_blank\" style=\"color:#019FFF\">UniprotID: ", UniprotID, "</a>"
      ))
   nodes <- mutate(nodes, group = sapply(1:nrow(nodes), function(i){
      sym <- nodes$id[i]
      if(is.na(nodes$BioplexID[i])){
         return("na")
      }else if(sym == current){
         return("current")
      }else if(sym %in% nearest){
         if(sym %in% bioplexData$SymbolA){
            return("bait_nearest")
         }else{
            return("nearest")
         }
      }else if(sym %in% bioplexData$SymbolA){
         return("bait_else")
      }else{
         return("else")
      }
   }))
   
   edges <- data.frame(from = df_gs$SymbolA, to = df_gs$SymbolB)
   edges <- edges %>% mutate(arrows = "to", shadow = TRUE, smooth = TRUE)

   v <- visNetwork(nodes, edges) %>%
      visGroups(groupname = "na", color = input$bioplex_na_color, shape = "diamond", shadow = list(enabled = TRUE)) %>%
      visGroups(groupname = "current", color = input$bioplex_current_color, shape = "square", shadow = list(enabled = TRUE)) %>%
      visGroups(groupname = "nearest", color = input$bioplex_nearest_color, shape = "diamond", shadow = list(enabled = TRUE)) %>%
      visGroups(groupname = "bait_nearest", color = input$bioplex_nearest_color, shadow = list(enabled = TRUE)) %>%
      visGroups(groupname = "else", color = input$bioplex_other_color, shape = "diamond", shadow = list(enabled = TRUE)) %>%
      visGroups(groupname = "bait_else", color = input$bioplex_bait_color, shadow = list(enabled = TRUE)) %>%
      visOptions(highlightNearest = TRUE) %>%
      visPhysics(stabilization = FALSE) %>%
      visEdges(smooth = FALSE) %>%
      visIgraphLayout(layout = input$bioplex_graph_layout) %>%
      visExport(type = tolower(input$bioplex_format), label = paste0("Download Graph (as ", input$bioplex_format, ")"), name = paste0(dataset$ID, "_", input$bioplex_graph_type, "_", values$activeMosaicID), 
            style="width:250px; border-color:#2ECC71; color:#000000; background-color:#45FF42")
   

   return(v)
})
