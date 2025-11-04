# Summary

The viewer app analyzes data and produces web-accessible plots. It is packaged as a Docker container to enforce permissions and ensure reproducible builds.

Categories of data include: Masspike API Data, Uniprot Annotations, and Bioplex Interactions.

## Masspike API

Protein, protein info, and peptide data are queried and cached as received. The API formats are then converted to a common format for use in the app.

## Uniprot Annotations

Annotations are stored in the `annotations.db`. Administrators may updated and configure this database via the `annotations` script.

## Bioplex Interactions

Bioplex data is automatically managed by the migration and validation scripts. It is stored in `bioplex.rds` and accessed for Human datasets.

# Data Flow

A user creates a viewer by accessing the TMT Editor webpage.

1. User authenticates with Masspike Credentials
2. Inputs viewer details
3. Protein-level data is queried
4. User finalizes viewer creation
5. Viewer configuration is stored in `database.db`
6. Protein data is cached for faster viewer loading

A user accesses a viewer in TMT Viewer.

1. The query string parameter `sessid` is used to retrieve viewer configuration
2. Protein, protein info, and peptide data are loaded by first checking the cache and then querying missing data
3. Data is formated to the common format and transformed by applying the specified configuration
4. Annotation data is retrieved and computed, if available
5. The web page is populated with interactive plots from the `tabs` folder

# Data Structures

The following is a reference to explain the fields of important data structures.

## Viewer Configuration

Database Columns are:
- QID:           Protein or SiteQuant ID
- Key:           Unique key passed in url to query a viewer
- Username:      Masspike username of viewer's creator
- Dataset:       Dataset's name
- Notes:         Dataset's notes
- NumSamples:    Total number of columns
- NumGroups:     Total number of groups
- PlexLevel:     Number of columns per plex
- NumPlexes:     Total number of plexes
- Server:        Full domain of server to query for data
- Species:       Dataset's species
- InitialID:     Initial site/protein to display (saved as mosaicID)
- AreReps:       Whether plexes can be summarized
- IsSiteQuant:   Whether QID describes a SiteQuant ID or a ProteinQuant ID
- Date:          Date this entry was added
- GroupNames:    User-defined group name, accessed with GroupIDs
- GroupColors:   User-defined group colors, accessed with GroupIDs
- GroupIDs:      Group index of each column
- ColumnClasses: Class membership of column (user-defined class names)
- ColumnNames:   Column names to be displayed
- ColumnIDs:     Column index of each column to display in display order

## Mosaic Dataset

The global variable `dataset` in `IsoParser` server code is a list that includes all of the fields in a `viewer configuration`, plus fields named `proteins`, `peptides`, and `normalization`.

Fields from the database include:

- ID          : numeric   [QID]
- key         : character [Key]
- username    : character [Username]
- name        : character [Dataset]
- notes       : character [Notes]
- numSamples  : numeric   <| Computed; number of ColumnIDs present
- numGroups   : numeric   [NumGroups]
- numPlexes   : numeric   [NumPlexes]
- plexLevel** : numeric   [PlexLevel]
- server      : character [Server]
- species     : character [Species]
- initialID   : character [InitialID]
- areReps     : logical   [AreReps]
- isSiteQuant : logical   [IsSiteQuant]
- date        : character [Date]
- groupNames  : character [GroupNames]    fromJSON / group
- groupColors : character [GroupColors]   fromJSON / group
- groups*     : numeric   [GroupIDs]      fromJSON / column
- classes*    : character [ColumnClasses] fromJSON / column
- names*      : character [ColumnNames]   fromJSON / column
- order       : numeric   [ColumnIDs]     fromJSON / column

Notes:

- \* Columns are reordered according to `ColumnIDs` to reflect user-specified column order
- `plexLevel` is deprecated and should not be used in new code
- See `src/IsoParser/lib/InitializeData.R` @ `parseMetadata for details

## Proteins

The `proteins` field of `dataset` is a list with `info` and `columns` fields. Info contains
qualitative data for each Site/Protein with the following columns:

- MosaicID: A unique identifier for use in TMT Mosaic
- UniprotID: Any database ID. If it is a UniprotID, Annotations are enabled
- GeneSymbol: Gene symbol associated with this protein
- Description: Short description of the protein
- Peptides: Number of peptides detected
- Site: Position of this site in the protein sequence [Site Quant Only]
- Sequence: AA sequence of fragment [Site Quant Only]

Quantitative data is stored as a matrix in the `columns` field with `MosaicID` rownames in the same
row order as the `info` data.frame.

## Peptides

The `peptides` field of `dataset` is also a list with `info` and `columns` fields.

- MosaicID: A unique identifier for use in TMT Mosaic
- UniprotID: Any database ID. If it is a UniprotID, Annotations are enabled
- GeneSymbol: Gene symbol associated with this protein
- Class: An abstract grouping mechanism
- SearchID: Identifier linking this peptide to a search
- PeptideID: Unique identifier of this peptide
- PeptideSequence: AA sequence of this peptide
- FirstScanNumber: Number of first scan with this peptide
- RunLoadPath: Path of data for the run
- Site: Position of this site in the protein sequence [Site Quant Only]

Quantitative data is stored as a list of matrices in the `columns` field with `MosaicID` rownames.
For example, a dataset with classes "A" and "B" may access quantitation data with this code:

```
dataset$peptides$columns[["A"]][<MosaicID_of_interest>, ]
```

## Normalization

The `normalization` field of `dataset` describes the kind of normalization applied and records the
exact normalization factors.

- `factors` - A vector of multiplicative normalization factors in the same order as `proteins$columns`
- `normBy` - What kind of normalization was applied. Options include {all, row, none}
  - all - The factors are proportional to column-wise sums
  - row - The factors are proportional to a protein quant row
  - none - No normalization was performed
- `normRow` - Is `NA` unless `normBy == "row"`, in which case, it is the `UniprotID` of the
              protein `factors` is derived from

## Attributes

`dataset` also stores data in attributes:
```
removed_contaminant_count <- attr(dataset, "contam")
removed_decoy_count <- attr(dataset, "decoy")
```
Where a
- `contaminant` is any entry whose `UniprotID` contains "contaminant" as a substring
- `decoy` contains "##" in either `UniprotID` or `GeneSymbol`

# Adminiatration

The web app is configured with `conf.yml`:

- Enumerates Masspike Servers with credentials that allow using `all_users=1` in API queries
- A web app admin password that allows access to all users' viewers in TMT Editor
- Enumerates species available for selection during viewer creation and guides the `annotations` utility

See `conf.tmp` for a more detailed description of the format.

Several scripts facilitate common administrative tasks:

- `annotations`: Download species annotation data from KEGG, GO, Reactome, and Interpro servers
- `genConf`: Interactively populate a `conf.yml` from the command line
- `migrateDB`: Check for new DB version and runs migration if detected
- `rename` : Easily edit species and server names in `database.db`
- `validateBuild`: Runs `migrateDB`, then checks that
    - `database.db` and `conf.yml` agree
    - `bioplex.rds` exists and is up-to-date
    - `annotations.db` exists with tables for each configured species

# Deployment

The `tests` folder contains scripts to verify backwards compatibility with existing viewers and
perform basic checks on the data loading phase.

## Unit Tests

Always run the unit tests before deploying to production.
```
tests/runTests
```
Add tests as bugs are discovered to detect future regressions.

## Adding Tests

Use the `snapshot` utility to automatically store a summarized form of the cache data and its loaded dataset
as a test. Note that the codebase MUST produce the correct output at the time of the snapshot, as it will be
saved into the test suite. The test name must be unique and a description of the test is required as well.

For example, on a working installation where "123MyKeY" is a key in `database.db`:
```
tests/snapshot example_test load 123MyKeY "a simple example dataset loads"
```

Creates
- Test entry in `tests/store/test_graph.csv`
- Input database entry at `tests/store/metadata/example_test.db.csv`*
- Input data at `tests/store/gfyRaw/example_test.csv`
- Output datafiles `proteins.csv` and `peptides.csv` under `tests/store/mosaic/example_test/`

\* `metadata/` files should be migrated if the DB format changes

See `tests/snapshot --help` for details.

Since the tests are all in csv format in `tests/store/`, it is also possible to manually create or modify
test input and output if necessary.

## DB Checks

The `healthcheck` utility displays warnings, errors and crash information for loading the selected
viewers.

It has several options. For example,
```
tests/healthcheck before=2020-10-01
```
loads all viewers created before October 10, 2020.

```
tests/healthcheck --help
```

displays all available options.

## Publish Changes

In the app folder, run `docker compose build` and `docker compose up -d` to update and start the app.

# Coding Standards

- Handle all exceptions to avoid crashing the app
- Write analyses in the `tabs` directory, with one tab per file
- Make general, reusable libraries if multiple parts of the app need the same functionality (never copy code)
- Always ensure backwards compatibility with the database and cached files via 1) migrations or 2) non-breaking changes
- If a major bug is reported, use `snapshot` and write a test to avoid future regressions
- Run `tests` before deploying new changes

