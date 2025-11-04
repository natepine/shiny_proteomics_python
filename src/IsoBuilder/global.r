source("../common/Globals.r")

libraries(c(
#    "sendmailR",
    "shinymanager"
))

navTabs = c(
#   "Tab Name"      = "tab_id",
    "Data Source"   = "source",
    "Display Names/Order" = "name",
    "Colors"        = "color",
    "Notes"         = "notes",
    "Finalize"      = "review"
)
navBttns = setNames(c(
#   "Button Label"
    "Check Data",
    "Save Names",
    "Save Colors",
    "Save Notes",
    "Create Viewer"
), navTabs)
