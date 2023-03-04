/* ANSI color Codes
---------------------------------------------------------------------------- */
// Used for determining the color of text output to the console
ANSI_RESET  = "\u001B[0m"
ANSI_GREEN  = "\u001B[32m"
ANSI_RED    = "\u001B[31m"
ANSI_YELLOW = "\u001B[33m"
ANSI_BLUE   = "\u001B[34m"
ANSI_GRAY   = "\033[0;90m"
ANSI_GREY   = ANSI_GRAY

/* Table logging function
---------------------------------------------------------------------------- */

// Master param logging function
def log_params(help_messages, parameters) {
    // Get the key names for all parameters and then print them out using a for loop
    // Recursively get all the keys listed in the help. We will output these values
    keys_in_help = [:]
    help_messages.values().toList().forEach(keys_in_help::putAll);

    // Begin
    log.info "Run parameters:"
    println ""

    // Table header
    table_log("parameter",  "value", false)
    table_log("=================", "===============================", false)

    // Table body
    // Get the key names for all parameters and then print them out using a for loop
    for (param in parameters.keySet()) {
        if (keys_in_help.containsKey(param)) {
            param_log(param)
        }
    }

    // Output the table footer
    println "----------------------------------------------------------------------------"
    table_log("User", workflow.userName, false)
    table_log("config", workflow.configFiles, false)
    table_log("Resume", workflow.resume, false)

    // Just for some space
    println ""
}

// Lower level helper for logging tables
def table_log(col1, col2, color = true) {
    tbl_format = color ? "${ANSI_YELLOW}%-35s%s${ANSI_RESET}" : "%-35s%s"
    log.info String.format(tbl_format, col1, col2)
}

// Lower level helper for logging params in a neat table
def param_log(param_name) {
    table_log("--" + param_name, params."$param_name", true)
}
