## DB version handling
suppressPackageStartupMessages({
   library(RSQLite)
   library(dplyr)
})

#' Convert from Jeremy's database to a keyed database for faster access and
#' redefine boolean columns as INT
#'
#' @param db Database connection to upgrade
version2to3 <- function(db){
   build <- paste0(
         "CREATE TABLE  `tmp` (
         `PQID`         INT NOT NULL,
         `KEY`          TEXT NOT NULL,
         `USERNAME`     TEXT NOT NULL,
         `DATANAME`     TEXT NOT NULL,
         `NOTES`        TEXT NOT NULL,
         `NUMSAMPLES`   INT NOT NULL,
         `NUMGROUPS`    INT NOT NULL,
         `PLEXLEVEL`    INT NOT NULL,
         `NUMPLEX`      INT NOT NULL,
         `SERVER`       TEXT NOT NULL,
         `SPECIES`      TEXT NOT NULL,
         `STARTPROT`    TEXT NOT NULL,
         `AREREPS`      INT NOT NULL,
         `PHOSPHO`      INT NOT NULL,
         `DATE`         DATETIME,
         PRIMARY KEY(`PQID`)
      );")

   names <- "PQID, KEY, USERNAME, DATANAME, NOTES, NUMSAMPLES, NUMGROUPS, PLEXLEVEL, NUMPLEX, SERVER, SPECIES, STARTPROT, DATE"
   copy <- paste0("INSERT INTO 'tmp' (", names, ", AREREPS, PHOSPHO, DATE) SELECT ", names, ", 0, 0, DATE FROM 'datainfo';")

   dbExecute(db, build)
   dbExecute(db, copy)
   dbExecute(db, "UPDATE tmp SET AREREPS = 1 WHERE KEY IN (SELECT KEY FROM datainfo WHERE AREREPS = 'TRUE')")
   dbExecute(db, "UPDATE tmp SET PHOSPHO = 1 WHERE KEY IN (SELECT KEY FROM datainfo WHERE PHOSPHO = 'TRUE')")
   dbExecute(db, "DROP TABLE datainfo;")
   dbExecute(db, "ALTER TABLE tmp RENAME TO datainfo;")
}

#' Convert to new database schema and rename cached data for compatibility with v3.3
#'
#' @param db Database connection to upgrade
#' @param cachePath Path to cache folder to upgrade. '/' will be appended.
#' @param confPath Path to conf.yml for appending API_DOMAIN to the SERVER column.
#' @param domain Ignore confPath and explicitly supply the domain name to append.
version3to3.3 <- function(db, cachePath = "temp", confPath = "data/conf.yml", domain = NULL){
   if(is.null(domain)){
      domain <- ""
      suppressPackageStartupMessages({
         library(yaml)
      })
      conf <- read_yaml(confPath)
      if(file.exists(conf)) {
         if("API_DOMAIN" %in% names(conf))
            domain <- conf$API_DOMAIN
         else
            warning("Defaulting to empty top-level domain, MANUALLY verify resulting server names in database.db")
      }else {
         stop("Error finding API_DOMAIN, no file at supplied path: ", normalizePath(confPath))
      }
   }

   #Build schema
   query <- paste0(
      "CREATE TABLE  `colMetadata` (",
      "`Key`          TEXT NOT NULL,", # Viewer key passed in url to query a viewer
      "`Class`        TEXT NOT NULL,", # Class membership of column (user-defined class names)
      "`Name`         TEXT NOT NULL,", # Column names to be displayed
      "`Color`        TEXT NOT NULL,", # Column colors
      "`GroupName`    TEXT NOT NULL,", # User-defined group name
      "`GroupNum`     INT NOT NULL,",  # Group membership of column
      "`Position`     INT NOT NULL",   # Custom column order
   ");")

   dbExecute(db, query)

   query <- paste0(
      "CREATE TABLE  `metadata` (",
      "`QID`          INT NOT NULL,",  # Protein or SiteQuant ID
      "`Key`          TEXT NOT NULL,", # Unique key passed in url to query a viewer
      "`Username`     TEXT NOT NULL,", # CORE username of viewer's creator
      "`Dataset`      TEXT NOT NULL,", # Dataset's name
      "`Notes`        TEXT NOT NULL,", # Dataset's notes
      "`NumSamples`   INT NOT NULL,",  # Total number of columns
      "`NumGroups`    INT NOT NULL,",  # Total number of groups
      "`PlexLevel`    INT NOT NULL,",  # Number of columns per plex
      "`NumPlexes`    INT NOT NULL,",  # Total number of plexes
      "`Server`       TEXT NOT NULL,", # Full domain of server to query for data
      "`Species`      TEXT,",          # Dataset's species
      "`InitialID`    TEXT NOT NULL,", # Initial site/protein to display (saved as mosaicID)
      "`AreReps`      INT NOT NULL,",  # Whether plexes can be summarized
      "`IsSiteQuant`  INT NOT NULL,",  # Whether QID describes a SiteQuant ID or a ProteinQuant ID
      "`Date`         TEXT NOT NULL,", # Date this entry was added
      "PRIMARY KEY(`Key`)",
   ");")

   dbExecute(db, query)

   insertMetadata <- paste0("INSERT INTO metadata ",
      "(QID, Key, Username, Dataset, Notes, NumSamples, NumGroups, PlexLevel, ",
      "NumPlexes, Server, Species, InitialID, AreReps, IsSiteQuant, Date) ",
      "VALUES ",
      "(@PQID, @KEY, @USERNAME, @DATANAME, @NOTES, @NUMSAMPLES, @NUMGROUPS, @PLEXLEVEL, ",
      "@NUMPLEX, @SERVER, @SPECIES, @STARTPROT, @AREREPS, @PHOSPHO, @DATE);")
   insertColMetadata <- paste0("INSERT INTO colMetadata ",
      "(Key, Class, Name, Color, GroupName, GroupNum, Position) ",
      "VALUES ",
      "(@KEY, @CLASS, @NAME, @COLOR, @GROUPNAME, @GROUPNUM, @POSITION);")

   # Load tables into R for processing
   colMetadata <- dbReadTable(db, "columninfo")
   metadata <- dbReadTable(db, "datainfo")

   # Handle missing or blank values, enforce data types
   metadata <- metadata %>% transmute(
      PQID = as.numeric(PQID),
      KEY = as.character(KEY),
      USERNAME = as.character(USERNAME),
      DATANAME = as.character(DATANAME),
      NOTES = as.character(NOTES),
      NUMSAMPLES = as.numeric(NUMSAMPLES),
      NUMGROUPS = as.numeric(NUMGROUPS),
      PLEXLEVEL = as.numeric(PLEXLEVEL),
      NUMPLEX = as.numeric(NUMPLEX),
      SERVER = as.character(ifelse(SERVER == "", domain, paste0(SERVER, ".", domain))),
      SPECIES = as.character(replace(SPECIES, SPECIES == "", NA)),
      STARTPROT = as.character(STARTPROT),
      AREREPS = as.logical(AREREPS),
      PHOSPHO = as.logical(PHOSPHO),
      DATE = as.character(replace(DATE, is.na(DATE), " - "))
   )

   keyMap <- setNames(metadata$KEY, metadata$PQID)
   colMetadata <- colMetadata %>% transmute(
      KEY = as.character(keyMap[as.character(PQID)]),
      CLASS = as.character(replace(CLASS, is.na(CLASS) | CLASS == "", "default")),
      NAME = as.character(ifelse(NAME == "", paste0("default_", POSITION), NAME)),
      COLOR = as.character(replace(COLOR, COLOR == "", "#777777")),
      GROUPNAME = as.character(paste("Group", GROUPNUM)),
      GROUPNUM = as.numeric(GROUPNUM),
      POSITION = as.numeric(POSITION)
   )

   # Save to database
   dbExecute(db, insertMetadata, params = metadata)
   dbExecute(db, insertColMetadata, params = colMetadata)

   dbExecute(db, "DROP TABLE datainfo;")
   dbExecute(db, "DROP TABLE columninfo;")

   # Rename cache
   renameFiles(cachePath, "proteins", metadata$PQID, metadata$PHOSPHO)
   renameFiles(cachePath, "peptides", metadata$PQID, metadata$PHOSPHO)
}

#' Rename cached files for v3.3. For example:
#' temp/proteins/123_proteins.rds -> temp/proteins/sq_123.rds
#'
#' @param path The cacheFolder
#' @param folder Which subdir to convert
#' @param qids Vector of Quant ids to search for
#' @param isSiteQuant Matching vector of booleans describing Quant id type
renameFiles <- function(path, folder, qids, isSiteQuant){
   library(fs)

   oldPaths <- paste0(path, "/", folder, "/", qids, "_", folder, ".rds")
   exists <- file_exists(oldPaths)
   if(any(exists)){
      newPaths <- paste0(path, "/", folder, "/", ifelse(isSiteQuant[exists], "s", "p"), "q_", qids[exists], ".rds")
      file_move(oldPaths[exists], newPaths)
   }
}

as.json <- function(x){
   as.character(toJSON(x))
}

#' Convert a v3.3 database to v3.4 schema
#' Merge colMetadata entries into metadata table
#'
#' @param db Database connection for 3.3 database
version3.3to3.4 <- function(db){
   library(jsonlite)

   #Build schema
   query <- paste0(
      "CREATE TABLE  `tmp` (",
      "`QID`          INT NOT NULL,",  # Protein or SiteQuant ID
      "`Key`          TEXT NOT NULL,", # Unique key passed in url to query a viewer
      "`Username`     TEXT NOT NULL,", # CORE username of viewer's creator
      "`Dataset`      TEXT NOT NULL,", # Dataset's name
      "`Notes`        TEXT NOT NULL,", # Dataset's notes
      "`NumSamples`   INT NOT NULL,",  # Total number of columns
      "`NumGroups`    INT NOT NULL,",  # Total number of groups
      "`PlexLevel`    INT NOT NULL,",  # Number of columns per plex
      "`NumPlexes`    INT NOT NULL,",  # Total number of plexes
      "`Server`       TEXT NOT NULL,", # Full domain of server to query for data
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
      "`ColumnIDs`     TEXT NOT NULL,", # Column index of each column in display order
      "PRIMARY KEY(`Key`)",
   ");")

   dbExecute(db, query)

   insert <- paste0("INSERT INTO tmp ",
      "(QID, Key, Username, Dataset, Notes, NumSamples, NumGroups, PlexLevel, ",
      "NumPlexes, Server, Species, InitialID, AreReps, IsSiteQuant, Date, ",
      "ColumnClasses, GroupNames, GroupColors, GroupIDs, ColumnNames, ColumnIDs) ",
      "VALUES ",
      "(@QID, @Key, @Username, @Dataset, @Notes, @NumSamples, @NumGroups, @PlexLevel, ",
      "@NumPlexes, @Server, @Species, @InitialID, @AreReps, @IsSiteQuant, @Date, ",
      "@ColumnClasses, @GroupNames, @GroupColors, @GroupIDs, @ColumnNames, @ColumnIDs);")

   # Load tables into R for processing
   colMetadata <- dbReadTable(db, "colMetadata")
   metadata <- dbReadTable(db, "metadata")

   # Just in case
   keys <- intersect(metadata$Key, unique(colMetadata$Key))
   if(length(keys) > 1){
      # There is data
      metadata <- metadata %>% filter(Key %in% keys)
      colMetadata <- colMetadata %>% filter(Key %in% keys)

      # Enforce types
      metadata <- metadata %>% transmute(
         QID = as.numeric(QID),
         Key = as.character(Key),
         Username = as.character(Username),
         Dataset = as.character(Dataset),
         Notes = as.character(Notes),
         NumSamples = as.numeric(NumSamples),
         NumGroups = as.numeric(NumGroups),
         PlexLevel = as.numeric(PlexLevel),
         NumPlexes = as.numeric(NumPlexes),
         Server = as.character(Server),
         Species = as.character(replace(Species, Species == "", NA)),
         InitialID = as.character(InitialID),
         AreReps = as.logical(AreReps),
         IsSiteQuant = as.logical(IsSiteQuant),
         Date = as.character(replace(Date, is.na(Date), " - ")),
      )

      flattenedCols <- data.frame()
      for(key in keys) {
         colKeyed <- colMetadata %>% filter(Key == key)

         groups <- distinct(colKeyed, GroupName, GroupNum, Color) %>%
            arrange(GroupNum)

         flattenedCols <- rbind(flattenedCols, 
            with(colKeyed, data.frame(
               Key = key,
               ColumnClasses = as.json(as.character(Class)),
               GroupNames = as.json(as.character(groups$GroupName)),
               GroupColors = as.json(as.character(groups$Color)),
               GroupIDs = as.json(as.numeric(GroupNum)),
               ColumnNames = as.json(as.character(Name)),
               ColumnIDs = as.json(as.numeric(Position))
            ))
         )
      }

      metadata <- metadata %>% inner_join(flattenedCols, by = "Key")
      dbExecute(db, insert, params = metadata)
   }

   dbExecute(db, "DROP TABLE metadata;")
   dbExecute(db, "DROP TABLE colMetadata;")
   dbExecute(db, "ALTER TABLE tmp RENAME TO metadata;")
}

#' Convert 3.4 'temp' cache to v3.5 server caches
#'
#' @param db Database connection for 3.4 database
version3.4to3.5 <- function(db){
   library(fs)

   viewers <- dbGetQuery(db, "SELECT QID, IsSiteQuant, Server FROM metadata")
   if(length(viewers) != 0) {
      moveCache <- function(viewers, subfolder) {
         dirs <- paste0(unique(viewers$Server), "/", subfolder)
         lapply(dirs, dir_create)
         paths <- paste0("/", subfolder, "/",
            ifelse(viewers$IsSiteQuant, "sq", "pq"), "_", viewers$QID, ".rds")

         oldPaths <- paste0("temp", paths)
         newPaths <- paste0(viewers$Server, paths)

         exist <- file_exists(oldPaths)
         if(any(exist)) {
            file_move(oldPaths[exist], newPaths[exist])
         }
      }
      moveCache(viewers, "proteins")
      moveCache(viewers, "peptides")
      if(dir.exists("temp")){
         dir_copy("temp",
            "temp~",
            overwrite = TRUE)
        dir_delete("temp")
      }
   }
}

withDB <- function(path, FUN) {
   db <- dbConnect(SQLite(), path)
   on.exit(dbDisconnect(db))
   FUN(db)
}

#             1  2    3    4    5
versions <- c(2, 3, 3.3, 3.4, 3.5)
upToDate <- length(versions)

#' Version checks
#'
#' @param db Database connection to version
#' @return Version index
getDBVersion <- function(db) {
   table_info <- dbGetQuery(db, "PRAGMA table_info(datainfo);")
   if(dbExistsTable(db, "datainfo")){
      areRepsType <- table_info[table_info$name == "AREREPS", "type"]
      if(areRepsType == "TEXT"){
         return(1)
      }
      return(2)
   }
   if(dbExistsTable(db, "colMetadata")){
      return(3)
   }
   if(dir.exists("temp")) {
      return(4)
   }
   return(upToDate)
}

# Automatic docker migration
if(!interactive()){
   dbPath <- "data/database.db"

   if(file.exists(dbPath)){
      # Only update existing data
      dbVersion <- withDB(dbPath, getDBVersion)
      message(paste("Detected database version", versions[[dbVersion]]))

      if(dbVersion != upToDate){
         # Only update out-of-date data
         tmpDB <- paste0(dbPath, ".tmp")
         file.copy(dbPath, tmpDB, overwrite = TRUE)

         db <- dbConnect(SQLite(), tmpDB)
         on.exit(dbDisconnect(db))

         if(dbVersion <= 1) version2to3(db)
         if(dbVersion <= 2) version3to3.3(db)
         if(dbVersion <= 3) version3.3to3.4(db)
         if(dbVersion <= 4) version3.4to3.5(db)

         #Sanity check
         currentVersion <- getDBVersion(db)
         if(currentVersion != upToDate){
            stop(
               paste0("Failed to migrate DB from version ", versions[[currentVersion]], " to ",
                  versions[[currentVersion + 1]], ". Check bad migration results at ", getwd(), "/",  tmpDB)
            )
         }else{
            # Free connection
            dbDisconnect(db)
            on.exit()

            # Replace db
            file.rename(dbPath, paste0(dbPath, ".bak"))
            file.rename(tmpDB, dbPath)
            message("Successfully migrated from db version ", versions[[dbVersion]], " to ", versions[[currentVersion]])
         }
      }else{
         message("The Database is already up-to-date")
      }

   }
}
