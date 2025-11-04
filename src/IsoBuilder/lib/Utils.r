# Authentication
check_gfy_credentials <- function(gfy_server) {
   if(is.null(gfy_server) && length(config$server_choices) == 1){
      #Single subdomain configurations don't have a server input
      gfy_server <- config$server_choices
   }
   verifySSL <- config$SERVERS[[gfy_server]]$VerifySSL
   url <- paste0("https://", gfy_server, "/gfy/www/modules/api/v1")
   function(username, password) {
      gfy.obj <- gfy.login(url, username, password, verifySSL)
      return(list(
         result = !is.null(gfy.obj),
         user_info = list(
            gfy.obj = gfy.obj,
            gfy_server = gfy_server,
            name = username)
      ))
   }
}


# Input validation
validateQID <- function(QID) {
   if(is.null(QID) || is.na(QID)){
      showModal(modalDialog(
         title = "Error | Blank Field",
         "One or more of the fields is blank. Please fill them and try submitting again."))
   }else if(is.na(as.numeric(QID))){
      showModal(modalDialog(
         title = "Error | Bad Params Input",
         "That Quant ID is not valid. Please check and try submitting again."))
   }else{
      return(TRUE)
   }
   return(FALSE)
}

validateCSVColumns <- function(columns, patterns, desc) {
   found <- sapply(patterns, function(pattern) {any(str_detect(columns, pattern))})

   if(!all(found)) {
      showModal(modalDialog(
         title = "Error | Missing Required Column(s)",
         paste(desc, "does not contain columns matching column(s):",
               paste(patterns[!found], collapse = ", "))))
      return(FALSE)
   }
   return (TRUE)
}

validateTables <- function(col_df, group_df, class_df, max_name_length){
   PASS <- "Passed"
   status <- c(
      "unique" = PASS,
      "length" = PASS
   )

   if(!any(col_df$Use)) {
      # Abort tests
      status <- setNames(
         rep("No selection: At least one sample must check 'Use'", length(status)),
         names(status)
      )
      return(status)
   }

   used <- which(col_df$Use)
   col_df <- col_df[used, ]

   all_names <- list(col_df$Name, group_df$`Group Name`, class_df$`Class Name`)

   duplicates <- sapply(all_names, duplicated) %>% unlist
   all_names <- unlist(all_names)
   if(any(duplicates)) {
      status[["unique"]] <- paste0(
         "Change duplicate names:</br>",
         paste0(unique(all_names[duplicates]), collapse = "</br>")
      )
   }

   too_long <- nchar(all_names) > max_name_length
   if(any(too_long)) {
      status[["length"]] <- paste0(
         "Limit names to ", max_name_length, " characters:</br>", 
         paste(all_names[too_long], collapse = "</br>")
      )
   }

   return(status)
}

validateNotes <- function(name, notes){
   if(name == "" || notes == ""){
      showModal(modalDialog(
         title = "Error | Empty Fields",
         "Please fill out all of the fields (dataset name and notes) before continuing.",
         easyClose = TRUE))
   }else{
      return(TRUE)
   }
   return(FALSE)
}

validateFinal <- function(metadata){
   requiredValues <- c('Username', 'Dataset', 'Notes', 'NumSamples',
      'NumGroups', 'PlexLevel', 'NumPlexes', 'Server', 'InitialID',
      'AreReps', 'IsSiteQuant', 'ColumnClasses', 'GroupNames', 'GroupColors',
      'GroupIDs', 'ColumnNames', 'ColumnIDs')

   isEmpty <- function(value) {
      any(is.null(value) | is.na(value) | value == "")
   }

   invalid <- sapply(requiredValues, function(key) isEmpty(metadata[[key]]))
   if(any(invalid)){
      showModal(modalDialog(
         title = "Error | Empty Fields",
         paste(
            "One or more required value is missing. Please review the following input(s):",
            paste(names(invalid[invalid]), collapse = ", ")
         ),
         easyClose = TRUE))
   }else{
      return(TRUE)
   }
   return(FALSE)
}

# Helper Functions
symbolIDs <- function(symbols, length.out) {
   n <- length.out
   nSymbols <- length(symbols)
   out <- c()

   cycle <- 0
   while(n > 0) {
      indices <- 0:n %% nSymbols + 1
      prefix <- rep(symbols[indices], each = nSymbols ^ cycle, length.out = length.out)

      out <- paste0(prefix, out)
      n <- n %/% nSymbols
      cycle <- cycle + 1
   }
   return(out)
}

# Database Functions
autoIncrID <- function(dbPath) {
   db <- dbConnect(SQLite(), dbPath)
   on.exit(dbDisconnect(db))

   if (!dbExistsTable(db, "uploadID")) {
      dbExecute(db,
         "CREATE TABLE uploadID (`id` INTEGER PRIMARY KEY AUTOINCREMENT, `key` TEXT NOT NULL)")
   }

   tmpKey <- generateStateKey()
   dbExecute(db, "INSERT INTO uploadID (key) VALUES (?);", tmpKey)
   id <- dbGetQuery(db, "SELECT id FROM uploadID WHERE key=?;", tmpKey)[1, ]
   dbExecute(db, "DELETE FROM uploadID WHERE key=?", tmpKey)

   return(id)
}

getViewerEntry <- function(db, QID, IsSiteQuant, Server) {
   # Construct the fetching query
   searchKey <- paste0("SELECT * FROM metadata WHERE QID = ? AND IsSiteQuant = ? AND Server = ?;")
   viewer <- dbGetQuery(db, searchKey, params = list(QID, as.numeric(IsSiteQuant), Server))

   if(length(viewer) == 0 || nrow(viewer) == 0) {
      return(NULL)
   }else if (nrow(viewer) > 1) {
      # Should never execute
      warning(paste0("Query '", searchKey, "' returned more than one Key:\n",
         paste0(viewer$Key, collapse = "\n"), "\nDefaulting to the first entry."))
      viewer <- viewer[1, ]
   }
   return(dbParseRow(viewer))
}

buildSchema <- function(db) {
   buildMetadata <- paste0(
      "CREATE TABLE  `metadata` (",
      "`QID`          INT NOT NULL,",  # Protein or SiteQuant ID
      "`Key`          TEXT NOT NULL,", # Unique key passed in url to query a viewer
      "`Username`     TEXT,",          # CORE username of viewer's creator
      "`Dataset`      TEXT NOT NULL,", # Dataset's name
      "`Notes`        TEXT NOT NULL,", # Dataset's notes
      "`NumSamples`   INT NOT NULL,",  # Total number of columns
      "`NumGroups`    INT NOT NULL,",  # Total number of groups
      "`PlexLevel`    INT NOT NULL,",  # Number of columns per plex
      "`NumPlexes`    INT NOT NULL,",  # Total number of plexes
      "`Server`       TEXT,",          # Full domain of server to query for data
      "`Species`      TEXT,",          # Dataset's species
      "`InitialID`    TEXT NOT NULL,", # Initial site/protein to display (saved as mosaicID)
      "`AreReps`      INT NOT NULL,",  # Whether plexes can be summarized
      "`IsSiteQuant`  INT NOT NULL,",  # Whether QID describes a SiteQuant ID or a ProteinQuant ID
      "`Date`         TEXT NOT NULL,", # Date this entry was added

# Serialized by as.character(jsonlite::toJSON(x)), unserialized by jsonlite::fromJSON(x)
      "`GroupNames`    TEXT NOT NULL,", # User-defined group name, accessed with GroupIDs
      "`GroupColors`   TEXT NOT NULL,", # User-defined group colors, accessed with GroupIDs
      "`GroupIDs`      TEXT NOT NULL,", # Group index of each column
      "`ColumnClasses` TEXT NOT NULL,", # Class membership of column (user-defined class names)
      "`ColumnNames`   TEXT NOT NULL,", # Column names to be displayed
      "`ColumnIDs`     TEXT NOT NULL,", # Column index of each column to display in display order
      "PRIMARY KEY(`Key`)",
   ");")

   dbExecute(db, buildMetadata)
}

findKey <- function(db, QID, IsSiteQuant, Server) {
   searchKey <- paste0("SELECT Key FROM metadata WHERE QID = ? AND IsSiteQuant = ? AND Server = ?;")
   key <- dbGetQuery(db, searchKey, params = c(QID, as.numeric(IsSiteQuant), Server))$Key

   if(length(key) == 0) {
      return(NULL)
   }else if (length(key) > 1) {
      # Should never execute
      warning(paste0("Query '", searchKey, "' returned more than one Key:\n",
         paste0(key, collapse = "\n"), "\nDefaulting to the first entry."))
      key <- key[[1]]
   }
   return(key)
}

clearKey <- function(db, key) {
   if (is.null(key) || !nzchar(key)) {
      cat(paste("Error clearing NULL or empty key"), file = stderr())
   }else{
      deleteMetadata <- paste0("DELETE FROM metadata WHERE Key = ?;")
      dbExecute(db, deleteMetadata, params = key)
   }
}

generateStateKey <- function(){
   digits <- 0:9
   paste0(c(
      sample(LETTERS, 5, replace = TRUE),
      sample(digits, 4, replace = TRUE),
      sample(LETTERS, 1, replace = TRUE),
      sample(digits, 2, replace = TRUE),
      sample(LETTERS, 5, replace = TRUE),
      sample(digits, 1, replace = TRUE),
      sample(LETTERS, 1, replace = TRUE),
      sample(digits, 2, replace = TRUE),
      sample(LETTERS, 4, replace = TRUE),
      sample(digits, 5, replace = TRUE),
      sample(LETTERS, 3, replace = TRUE)), collapse = "")
}

asJSON <- function(data, method) {
   for(target in names(method)){
      as.type <- method[[target]]
      data[[target]] <- data[[target]] %>% as.type %>% toJSON %>% as.character
   }
   return(data)
}

dbParseRow <- function(row) {
   # Logicals are saved as numerics
   logicals <- c("IsSiteQuant", "AreReps")
   row[logicals] <- row[logicals] %>% as.logical

   # Restore JSON columns
   row <- as.list(row)
   for(target in names(JSON_Vars)){
      as.type <- JSON_Vars[[target]]
      row[[target]] <- row[[target]] %>% fromJSON %>% as.type
   }
   return(row)
}

saveToDatabase <- function(metadata, dbPath) {
   # Connect to the database
   db <- dbConnect(SQLite(), dbPath)
   on.exit(dbDisconnect(db))

   if (!dbExistsTable(db, "metadata")) {
      buildSchema(db)
   }

   QID <- metadata$QID
   key <- metadata$Key
   if(!is.null(key)){
      # Overwrite an old key
      message("Reusing key for ", ifelse(metadata$IsSiteQuant, "S", "P"), "QID: ", QID, '\n', key)
      # Clear all info
      clearKey(db, key)
   }else{
      # Write to a new key
      key <- generateStateKey()
      message("Wrote new key: ", key)
   }

   if(!is.null(key) && nzchar(key)){
      metadata$Key <- key
      metadata <- asJSON(metadata, JSON_Vars)

      # Save to database
      saveMetadata(metadata, db)

      return(key)
   }else{
      msg <- paste0("Failed to generate key for Quant ID: ", QID)
      showNotification(msg, duration = NULL)
      cat(paste0(
         msg, "\n",
         keyValueStr(metadata, collapse = "\n")),
         file = stderr())
   }
   return(NULL)
}

#' @param params A dataframe or list with names identical to: c(QID, Key, Username, Dataset, Notes, NumSamples, NumGroups, PlexLevel,
#'    NumPlexes, Server, Species, InitialID, AreReps, IsSiteQuant, ColumnClasses, GroupNames, GroupColors, GroupIDs, ColumnNames, ColumnIDs)
#' @param db SQLite database connection
saveMetadata <- function(params, db) {
   # Construct a query to Insert a new row
   insert <- paste0("INSERT INTO metadata ",
      "(QID, Key, Username, Dataset, Notes, NumSamples, NumGroups, PlexLevel, ",
      "NumPlexes, Server, Species, InitialID, AreReps, IsSiteQuant, Date, ",
      "ColumnClasses, GroupNames, GroupColors, GroupIDs, ColumnNames, ColumnIDs) ",
      "VALUES ",
      "(@QID, @Key, @Username, @Dataset, @Notes, @NumSamples, @NumGroups, @PlexLevel, ",
      "@NumPlexes, @Server, @Species, @InitialID, @AreReps, @IsSiteQuant, date('now'), ",
      "@ColumnClasses, @GroupNames, @GroupColors, @GroupIDs, @ColumnNames, @ColumnIDs);")
   
   message("Save Metadata:\n")
   message(keyValueStr(params, collapse = "\n"))

   dbExecute(db, insert, params = params)
}

getAdminTableData <- function(dbPath) {
   if(!file.exists(dbPath)) return(NULL)

   getTableData <- paste("SELECT Date, QID, Username, Dataset, Key, Notes, Server, Species, NumSamples,",
      "NumGroups, PlexLevel, AreReps, IsSiteQuant, InitialID FROM metadata;")

   db <- dbConnect(SQLite(), dbPath, flags = SQLITE_RO)
   on.exit(dbDisconnect(db))

   response <- dbGetQuery(db, getTableData)

   if(length(response) == 0 || nrow(response) == 0){
      return(NULL)
   }
   return(response)
}

deleteViewer <- function(key, dbPath) {
   if(file.exists(dbPath)) {
      deleteKeyRow <- paste("DELETE FROM metadata WHERE Key = ?;")

      db <- dbConnect(SQLite(), dbPath, flags = SQLITE_RW)
      on.exit(dbDisconnect(db))

      dbExecute(db, deleteKeyRow, params = key)
   }
}
