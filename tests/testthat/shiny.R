
library(shiny)
library(Rcwl)
input1 <- InputParam(id = "sth")
echo <- cwlParam(baseCommand = "echo", inputs = InputParamList(input1))

e1 <- InputParam(id = "flag", type = "boolean", prefix = "-f", doc = "boolean flag")
e2 <- InputParam(id = "string", type = "string", prefix = "-s")
e3 <- InputParam(id = "int", type = "int", prefix = "-i", default = 123)
e4 <- InputParam(id = "file", type = "File", prefix = "--file=", separate = FALSE)
e5 <- InputParam(id = "array", type = "string[]", prefix = "-A", doc = "separated by comma")
echoA <- cwlParam(baseCommand = "echo", id = "mulEcho", label = "Test parameter types",
                  inputs = InputParamList(e1, e2, e3, e4, e5),
                  stdout = "output.txt")

.inputUI <- function(cwl, inputList, upload=FALSE){
    ilist <- inputs(cwl)
    dList <- lapply(ilist, function(x){
        if(x@id %in% names(inputList)){
            d <- selectInput(inputId = x@id,
                             label = paste0(x@id, " (", x@type, ")"),
                             choices = inputList[[x@id]])
        }else{
            if(class(x@type) == "InputArrayParam"){
                d <- textAreaInput(inputId = x@id,
                                   label = paste0(x@id, " (array)"),
                                   value = x@default)
            }else if(x@type == "boolean"){
                d <- selectInput(inputId = x@id,
                                 label = paste0(x@id, " (", x@type, ")"),
                                 choices = list("TRUE" = TRUE, "FALSE" = FALSE))
            }else if(x@type %in% c("string", "int", "File", "Directory")){
                if(x@type == "File" & upload){
                    d <- fileInput(inputId = x@id,
                                   label = paste0(x@id, " (", x@type, ")"))
                }else{
                    d <- textInput(inputId = x@id,
                                   label = paste0(x@id, " (", x@type, ")"),
                                   value = x@default)
                }
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

cwlShiny <- function(cwl, inputList = list(), upload = FALSE, ...){
    tList <- titlePanel(cwl@id)
    if(length(cwl@label) > 0) tList <- list(tList, h3(cwl@label))
    divList <- .inputUI(cwl, inputList, upload)

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

tmp1 <- tempfile()
tmp2 <- tempfile()
file.create(tmp1)
file.create(tmp2)
inputList <- list(file = c(tmp1, tmp2))
echoS <- cwlShiny(echoA, inputList)
runApp(echoS)

inputList <- list()
