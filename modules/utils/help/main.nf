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

/* Help Functions
---------------------------------------------------------------------------- */
def show_help(messages) {
    messages.forEach({key, value ->
        if (value instanceof LinkedHashMap) {
            println ""
            println "${ANSI_GREEN}${key}:${ANSI_RESET}"
            show_help(value)
        } else {
            println String.format("    ${ANSI_YELLOW}--%-35s%s${ANSI_RESET}", key, value)
        }
    })
}
