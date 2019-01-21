.inputUI <- function(cwl, inputList){
    ilist <- inputs(cwl)
    dList <- lapply(ilist, function(x){
        if(x@id %in% names(inputList)){
            d <- selectInput(inputId = x@id, label = paste0(x@id, " (", x@type, ")"), choices = inputList[[x@id]])
        }else{
            if(class(x@type) == "InputArrayParam"){
                d <- textAreaInput(inputId = x@id, label = paste0(x@id, " (array)"), value = x@default)
            }else if(x@type == "boolean"){
                d <- selectInput(inputId = x@id, label = paste0(x@id, " (", x@type, ")"), choices = list("TRUE" = TRUE, "FALSE" = FALSE))
            }else if(x@type %in% c("string", "int", "File", "Directory")){
                d <- textInput(inputId = x@id, label = paste0(x@id, " (", x@type, ")"), value = x@default)
            }else if(grepl("\\[", x@type)){
                d <- textAreaInput(inputId = x@id, label = paste0(x@id, " (", x@type, ")"), value = x@default)
            }else{
                d <- textInput(inputId = x@id, label = paste0(x@id, " (", x@type, ")"), value = x @default)
            }
        }
        
        if(length(x@doc) > 0){
            list(d, helpText(x@doc))
        }else{
            d
        }
    })
    dList
}

#' cwlShiny
#'
#' Function to generate shiny app automaticlly for a `cwlParam` object.
#' @param cwl A cwlParam object.
#' @param inputList a list of choices for the inputs of cwl object. The name of the list must match the inputs of the cwl object.
#' @param ... More options for `runCWL`.
#' @import Rcwl
#' @import shiny
#' @export
#' @examples
#' input1 <- InputParam(id = "sth")
#' echo <- cwlParam(baseCommand = "echo", inputs = InputParamList(input1))
#' echoApp <- cwlShiny(echo)

cwlShiny <- function(cwl, inputList = list(), ...){
    stopifnot(class(cwl) == "cwlParam")
    tList <- titlePanel(cwl@id)
    if(length(cwl@label) > 0) tList <- list(tList, h3(cwl@label))
    divList <- .inputUI(cwl, inputList)

    ui <- fluidPage(                
        tList,
        sidebarLayout(
            sidebarPanel(divList),
            mainPanel(
                actionButton("run", "Run"),
                tabsetPanel(
                    tabPanel("Output", verbatimTextOutput("outPath")),
                    tabPanel("Command", verbatimTextOutput("command")),
                    tabPanel("Log", verbatimTextOutput("log"))
                )
            )
        )
    )

    ilist <- inputs(cwl)
    server <- function(input, output){
        res <- eventReactive(input$run, {
            for(x in names(ilist)){
                if(class(ilist[[x]]@type) == "InputArrayParam" ||
                   grepl("\\[", ilist[[x]]@type)){
                    v <- unlist(strsplit(input[[x]], split = ","))
                    eval(parse(text=paste0("cwl$",x," <- v")))
                }else{    
                    eval(parse(text=paste0("cwl$",x," <- input$", x)))
                }
            }
            print("Running...")
            runCWL(cwl, ...)
        })
        output$outPath <- renderText({
            paste(res()$output, collapse = "\n")
        })
        output$command <- renderText({
            paste(res()$command, collapse = "\n")
        })
        output$log <- renderText({
            paste(res()$log, collapse = "\n")
        })}

    shinyApp(ui, server)
}
