server <- function(input, output, session) {
  observeEvent(input$help, introjs(session, options = list("showBullets" = "false", "showProgress" = "true", 
                                                          "showStepNumbers" = "false","nextLabel" = "Next",
                                                          "prevLabel" = "Prev", "skipLabel" = "Skip")))
  react.input <- reactive({
    req(input$exp_file)
    exp.file <- input$exp_file
    ext <- tools::file_ext(exp.file$name)
    shiny::validate(
      need(ext == "tsv",
           "Please upload a valid UCSC Xena gene expression TSV file!")
    )
    read_delim(exp.file$datapath, "\t", escape_double = FALSE, trim_ws = TRUE)
  })
  
  react.dataset <- reactive({
    shiny::validate(
      need(input$dataset %in% react.input()$`cancer type abbreviation`,
           "Selected TCGA dataset not found in input file!")
    )
    sub.pancan(react.input(), input$dataset)
  })

  react.surv <- reactive({
    set.surv(react.dataset()$`cancer type abbreviation`[1])
  })
  
  observe({
    updateSliderInput(session, "n_genes", max = length(3:(ncol(react.dataset())-2)), step = 1)
  })
  
  gene.choices <- reactiveValues(gene_vec = c())
  
  observe({
    if(input$selectall > 0) {
      if(input$selectall %% 2 == 0 & input$n_genes == length(3:(ncol(react.dataset())-2) )) {
        updateSelectizeInput(session, "genes", selected = gene.choices$gene_vec)
      } else{
        updateSelectizeInput(session, "genes", selected = NULL)
      }
    } 
  })

  observe({
    gene.choices$gene_vec <- set.genes(react.dataset())
    updateSelectizeInput(session, "genes", choices = gene.choices$gene_vec, options = list(maxItems = input$n_genes))
  })
  
  observe({
    updateSelectizeInput(session, "fold_genes", choices = setdiff(gene.choices$gene_vec, input$genes))
  })

  react.exp <- reactive({
    shiny::validate(
      need(input$genes %in% colnames(react.dataset()),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!any(colSums(is.na(react.dataset()[, input$genes])) >= nrow(react.dataset()[, input$genes])*99/100),
           paste("Selected gene:", names(which(colSums(is.na(react.dataset()[, input$genes])) >= nrow(react.dataset()[, input$genes])*99/100)), "contain missing expressions for this dataset!"))
    )
    shiny::validate(
      need(!(is.null(input$genes)),
           "No genes selected!")
    )

    set.exp(react.dataset(), input$genes)
  })
  react.fold <- reactive({
    shiny::validate(
      need(input$genes %in% colnames(react.dataset()),
           "Selected genes not found in input dataset!")
    )
    set.exp(react.dataset(), input$fold_genes)
  })

  react.exp.fold <- reactive({
    req(input$genes, input$fold_genes)
    shiny::validate(
      need(!any(colSums(is.na(react.dataset()[, input$fold_genes])) >= nrow(react.dataset()[, input$fold_genes])*99/100),
           paste("Selected fold gene(s):", names(which(colSums(is.na(react.dataset()[, input$fold_genes])) >= 
                                                         nrow(react.dataset()[, input$fold_genes])*99/100)), "contain missing expressions for this dataset!"))
    )
    fold.exprs <- cbind(react.fold()$sample, react.exp()[2:ncol(react.exp())] - rowSums(react.fold()[2:ncol(react.fold())]))
    colnames(fold.exprs) <- "sample"
    set.foldchange(react.dataset(), input$genes, fold.exprs)
  })

##############################################
  react.cox.os <- reactive({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("OS", "OS.time")])))) >= nrow(react.surv()[, c("OS", "OS.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    if(!is.null(input$fold_genes)) {
      multi.cox.model.os(react.surv(), react.exp.fold())
    } else {
      multi.cox.model.os(react.surv(), react.exp())
    }
  })

  react.cox.dss <- reactive({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("DSS", "DSS.time")])))) >= nrow(react.surv()[, c("DSS", "DSS.time")])*(99/100),
           "Insufficient data for disease-specific survival!")
    )
    if(!is.null(input$fold_genes)) {
      multi.cox.model.dss(react.surv(), react.exp.fold())
    } else {
      multi.cox.model.dss(react.surv(), react.exp())
    }
  })

  react.cox.dfi <- reactive({
    req(input$genes)
    shiny::validate(
       need(!length(which(is.na(rowSums(react.surv()[, c("DFI", "DFI.time")])))) >= nrow(react.surv()[, c("DFI", "DFI.time")])*(99/100),
            "Insufficient data for disease-free survival!")
    )
    if(!is.null(input$fold_genes)) {
      multi.cox.model.dfi(react.surv(), react.exp.fold())
    } else {
      multi.cox.model.dfi(react.surv(), react.exp())
    }
  })

  react.cox.pfi <- reactive({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("PFI", "PFI.time")])))) >= nrow(react.surv()[, c("PFI", "PFI.time")])*(99/100),
           "Insufficient data for progress-free survival")
    )
    if(!is.null(input$fold_genes)) {
      multi.cox.model.pfi(react.surv(), react.exp.fold())
    } else {
      multi.cox.model.pfi(react.surv(), react.exp())
    }
  })
########################################
  react.cox.os.data <- reactive({
    req(input$genes)
    if(!is.null(input$fold_genes)) {
      exp.data <- react.exp.fold()
      shiny::validate(
        need(any(colSums(is.na(exp.data[2:ncol(exp.data)])) < nrow(exp.data[2:ncol(exp.data)])*99/100),
             "Missing expressions for one or more gene(s)!")
      )
      data.os <- merge(react.surv(), exp.data)
      data.os[, c(1, 3:4, 12:ncol(data.os))]
    } else {
      exp.data <- react.exp()
      shiny::validate(
        need(any(colSums(is.na(exp.data[2:ncol(exp.data)])) < nrow(exp.data[2:ncol(exp.data)])*99/100),
             "Missing expressions for one or more gene(s)!")
      )
      data.os <- merge(react.surv(), exp.data)
      data.os[, c(1, 3:4, 12:ncol(data.os))]
    }
  })

  react.cox.dss.data <- reactive({
    req(input$genes)
    if(!is.null(input$fold_genes)) {
      exp.data <- react.exp.fold()
      shiny::validate(
        need(any(colSums(is.na(react.exp.fold()[2:ncol(react.exp.fold())])) < nrow(react.exp.fold()[2:ncol(react.exp.fold())])*99/100),
             "Missing expressions for one or more gene(s)!")
      )
      data.dss <- merge(react.surv(), exp.data)
      data.dss[, c(1, 5:6, 12:ncol(data.dss))]
    } else {
      exp.data <- react.exp()
      shiny::validate(
        need(any(colSums(is.na(react.exp()[2:ncol(react.exp())])) < nrow(react.exp()[2:ncol(react.exp())])*99/100),
             "Missing expressions for one or more gene(s)!")
      )
      data.dss <- merge(react.surv(), exp.data)
      data.dss[, c(1, 5:6, 12:ncol(data.dss))]
    }
  })

  react.cox.dfi.data <- reactive({
    req(input$genes)
    if(!is.null(input$fold_genes)) {
      exp.data <- react.exp.fold()
      shiny::validate(
        need(any(colSums(is.na(react.exp.fold()[2:ncol(react.exp.fold())])) < nrow(react.exp.fold()[2:ncol(react.exp.fold())])*99/100),
             "Missing expressions for one or more gene(s)!")
      )
      data.dfi <- merge(react.surv(), exp.data)
      data.dfi[, c(1, 7:8, 12:ncol(data.dfi))]
    } else {
      exp.data <- react.exp()
      shiny::validate(
        need(any(colSums(is.na(react.exp()[2:ncol(react.exp())])) < nrow(react.exp()[2:ncol(react.exp())])*99/100),
             "Missing expressions for one or more gene(s)!")
      )
      data.dfi <- merge(react.surv(), exp.data)
      data.dfi[, c(1, 7:8, 12:ncol(data.dfi))]
    }
  })
  
  react.cox.pfi.data <- reactive({
    req(input$genes)
    if(!is.null(input$fold_genes)) {
      exp.data <- react.exp.fold()
      shiny::validate(
        need(any(colSums(is.na(react.exp.fold()[2:ncol(react.exp.fold())])) < nrow(react.exp.fold()[2:ncol(react.exp.fold())])*99/100),
             "Missing expressions for one or more gene(s)!")
      )
      data.pfi <- merge(react.surv(), exp.data)
      data.pfi[, c(1, 9:10, 12:ncol(data.pfi))]
    } else {
      exp.data <- react.exp()
      shiny::validate(
        need(any(colSums(is.na(react.exp()[2:ncol(react.exp())])) < nrow(react.exp()[2:ncol(react.exp())])*99/100),
             "Missing expressions for one or more gene(s)!")
      )
      data.pfi <- merge(react.surv(), exp.data)
      data.pfi[, c(1, 9:10, 12:ncol(data.pfi))]
    }
  })

#################################
  react.cph.os <- reactive({
    cph(react.cox.os()$formula, data = react.cox.os.data(), x = TRUE, y = TRUE, surv = TRUE)
  })

  react.cph.dss <- reactive({
    cph(react.cox.dss()$formula, data = react.cox.dss.data(), x = TRUE, y = TRUE, surv = TRUE)
  })

  react.cph.dfi <- reactive({
    cph(react.cox.dfi()$formula, data = react.cox.dfi.data(), x = TRUE, y = TRUE, surv = TRUE)
  })

  react.cph.pfi <- reactive({
    cph(react.cox.pfi()$formula, data = react.cox.pfi.data(), x = TRUE, y = TRUE, surv = TRUE)
  })

#########################################################
  output$data_os <- renderTable({
    surv.data <- react.cox.os()$model
    surv.data[2:ncol(surv.data)]
  })
  
  output$data_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_genes", ".tsv", sep = "")},
    content = function(file) {
      write.table(react.cox.os()$model, file)}
  )

  output$data_dss <- renderTable({
    surv.data <- react.cox.dss()$model
    surv.data[2:ncol(surv.data)]
  })

  output$data_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_genes", ".tsv", sep = "")},
    content = function(file) {
      write.table(react.cox.dss()$model, file)}
  )
  
  output$data_dfi <- renderTable({
    surv.data <- react.cox.dfi()$model
    surv.data[2:ncol(surv.data)]
  })
  
  output$data_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_genes", ".tsv", sep = "")},
    content = function(file) {
      write.table(react.cox.dfi()$model, file)}
  )
  
  output$data_pfi <- renderTable({
    surv.data <- react.cox.pfi()$model
    surv.data[2:ncol(surv.data)]
  })
  
  output$data_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_genes", ".tsv", sep = "")},
    content = function(file) {
      write.table(react.cox.pfi()$model, file)}
  )
  
  output$data_box_os <- renderPlot({
    boxplot(react.cox.os.data()[4:ncol(react.cox.os.data())], col = "pink", pch = 2, las = 2)
  })
  
  output$data_box_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_genes", "_box", ".png", sep = "")},
    content = function(file) {
      png(filename = file)
      boxplot(react.cox.os.data()[4:ncol(react.cox.os.data())], col = "pink", pch = 2, las = 2)
      dev.off()
    }
  )
  
  output$data_box_dss <- renderPlot({
    boxplot(react.cox.dss.data()[4:ncol(react.cox.dss.data())], col = "pink", pch = 2, las = 2)
  })
  
  output$data_box_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_genes", "_box", ".png", sep = "")},
    content = function(file) {
      png(filename = file)
      boxplot(react.cox.dss.data()[4:ncol(react.cox.dss.data())], col = "pink", pch = 2, las = 2)
      dev.off()
    }
  )
  
  output$data_box_dfi <- renderPlot({
    boxplot(react.cox.dfi.data()[4:ncol(react.cox.dfi.data())], col = "pink", pch = 2, las = 2)
  })
  
  output$data_box_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_genes", "_box", ".png", sep = "")},
    content = function(file) {
      png(filename = file)
      boxplot(react.cox.dfi.data()[4:ncol(react.cox.dfi.data())], col = "pink", pch = 2, las = 2)
      dev.off()
    }
  )
  
  output$data_box_pfi <- renderPlot({
    boxplot(react.cox.pfi.data()[4:ncol(react.cox.pfi.data())], col = "pink", pch = 2, las = 2)
  })
  
  output$data_box_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_genes", "_box", ".png", sep = "")},
    content = function(file) {
      png(filename = file)
      boxplot(react.cox.pfi.data()[4:ncol(react.cox.pfi.data())], col = "pink", pch = 2, las = 2)
      dev.off()
    }
  )
  
  output$forest_os <- renderPlot({
    shiny::validate(
      need(react.cox.os()$iter <= 20,
           "Forest plot not drawn due to unconverging model!")
    )
    shiny::validate(
      need(all(!is.infinite(react.cox.os()$coefficients)),
           "Forest plot not drawn due to infinite coeff(s)!")
    )
    ggforest(react.cox.os(), data = react.cox.os.data())
  })
  
  output$forest_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_genes", "_forest", ".png", sep = "")},
    content = function(file) {
      shiny::validate(
        need(react.cox.os()$iter <= 20,
             "Forest plot not drawn due to unconverging model!")
      )
      shiny::validate(
        need(all(!is.infinite(react.cox.os()$coefficients)),
             "Forest plot not drawn due to infinite coeff(s)!")
      )
      ggsave(file, ggforest(react.cox.os(), data = react.cox.os.data()), dpi = 320)
    }
  )
    
  output$forest_dss <- renderPlot({
    shiny::validate(
      need(react.cox.dss()$iter <= 20,
           "Forest plot not drawn due to unconverging model!")
    )
    shiny::validate(
      need(all(!is.infinite(react.cox.dss()$coefficients)),
           "Forest plot not drawn due to infinite coeff(s)!")
    )
    ggforest(react.cox.dss(), data = react.cox.dss.data())
  })
  
  output$forest_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_genes", "_forest", ".png", sep = "")},
    content = function(file) {
      shiny::validate(
        need(react.cox.dss()$iter <= 20,
             "Forest plot not drawn due to unconverging model!")
      )
      shiny::validate(
        need(all(!is.infinite(react.cox.dss()$coefficients)),
             "Forest plot not drawn due to infinite coeff(s)!")
      )
      ggsave(file, ggforest(react.cox.dss(), data = react.cox.dss.data()), dpi = 320)
    }
  )
  
  output$forest_dfi <- renderPlot({
    shiny::validate(
      need(react.cox.dfi()$iter <= 20,
           "Forest plot not drawn due to unconverging model!")
    )
    shiny::validate(
      need(all(!is.infinite(react.cox.dfi()$coefficients)),
           "Forest plot not drawn due to infinite coeff(s)!")
    )
    ggforest(react.cox.dfi(), data = react.cox.dfi.data())
  })

  output$forest_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_genes", "_forest", ".png", sep = "")},
    content = function(file) {
      shiny::validate(
        need(react.cox.dfi()$iter <= 20,
             "Forest plot not drawn due to unconverging model!")
      )
      shiny::validate(
        need(all(!is.infinite(react.cox.dfi()$coefficients)),
             "Forest plot not drawn due to infinite coeff(s)!")
      )
      ggsave(file, ggforest(react.cox.dfi(), data = react.cox.dfi.data()), dpi = 320)
    }
  )
  
  output$forest_pfi <- renderPlot({
    shiny::validate (
      need(react.cox.pfi()$iter <= 20,
           "Forest plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.pfi()$coefficients)),
           "Forest plot not drawn due to infinite coeff(s)!")
    )
    ggforest(react.cox.pfi(), data = react.cox.pfi.data())
  })
  
  output$forest_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_genes", "_forest", ".png", sep = "")},
    content = function(file) {
      shiny::validate(
        need(react.cox.pfi()$iter <= 20,
             "Forest plot not drawn due to unconverging model!")
      )
      shiny::validate(
        need(all(!is.infinite(react.cox.pfi()$coefficients)),
             "Forest plot not drawn due to infinite coeff(s)!")
      )
      ggsave(file, ggforest(react.cox.pfi(), data = react.cox.pfi.data()), dpi = 320)
    }
  )
  
  output$schoen_os <- renderPlot({
    shiny::validate (
      need(react.cox.os()$iter <= 20,
           "Schoenfeld plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.os()$coefficients)),
           "Schoenfeld plot not drawn due to infinite coeff(s)!")
    )
    ggcoxzph(cox.zph(react.cox.os()), legend.labs = names(react.cox.os()$coefficients),
              font.main = 9, font.submain = 9, font.legend = 12, font.tickslab = 7, font.x = 10, font.y = 10, point.size = 0.3)
  })
  
  output$schoen_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_genes", "_schoen", ".png", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.os()$iter <= 20,
             "Schoenfeld plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.os()$coefficients)),
             "Schoenfeld plot not drawn due to infinite coeff(s)!")
      )
      ggsave(file, print(ggcoxzph(cox.zph(react.cox.os()), legend.labs = names(react.cox.os()$coefficients), 
               font.main = 9, font.submain = 9, font.legend = 12, font.tickslab = 7,
               font.x = 10, font.y = 10, point.size = 0.3)), width = 20, units = "cm", dpi = 320)
      }
  )
  
  output$schoen_dss <- renderPlot({
    shiny::validate (
      need(react.cox.dss()$iter <= 20,
           "Schoenfeld plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.dss()$coefficients)),
           "Schoenfeld plot not drawn due to infinite coeff(s)!")
    )
    ggcoxzph(cox.zph(react.cox.dss()), legend.labs = names(react.cox.dss()$coefficients), font.main = 9, font.submain = 9, 
             font.legend = 12, font.tickslab = 7, font.x = 10, font.y = 10, point.size = 0.3)
  })
  
  output$schoen_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_genes", "_schoen", ".png", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.dss()$iter <= 20,
             "Schoenfeld plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.dss()$coefficients)),
             "Schoenfeld plot not drawn due to infinite coeff(s)!")
      )
      ggsave(file, print(ggcoxzph(cox.zph(react.cox.dss()), legend.labs = names(react.cox.dss()$coefficients), 
                            font.main = 9, font.submain = 9, font.legend = 12, font.tickslab = 7, font.x = 10, 
                            font.y = 10, point.size = 0.3)), width = 20, units = "cm", dpi = 320)
    }
  )
  
  output$schoen_dfi <- renderPlot({
    shiny::validate (
      need(react.cox.dfi()$iter <= 20,
           "Schoenfeld plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.dfi()$coefficients)),
           "Schoenfeld plot not drawn due to infinite coeff(s)!")
    )
    ggcoxzph(cox.zph(react.cox.dfi()), legend.labs = names(react.cox.dfi()$coefficients), 
             font.main = 9, font.submain = 9, font.legend = 12, font.tickslab = 7,
             font.x = 10, font.y = 10, point.size = 0.3)
  })
  
  output$schoen_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_genes", "_schoen", ".png", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.dfi()$iter <= 20,
             "Schoenfeld plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.dfi()$coefficients)),
             "Schoenfeld plot not drawn due to infinite coeff(s)!")
      )
      ggsave(file, print(ggcoxzph(cox.zph(react.cox.dfi()), legend.labs = names(react.cox.dfi()$coefficients), 
                            font.main = 9, font.submain = 9, font.legend = 12, font.tickslab = 7,
                            font.x = 10, font.y = 10, point.size = 0.3)), width = 20, units = "cm", dpi = 320)
    }
  )
  
  output$schoen_pfi <- renderPlot({
    shiny::validate (
      need(react.cox.pfi()$iter <= 20,
           "Schoenfeld plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.pfi()$coefficients)),
           "Schoenfeld plot not drawn due to infinite coeff(s)!")
    )
    ggcoxzph(cox.zph(react.cox.pfi()), legend.labs = names(react.cox.pfi()$coefficients), 
             font.main = 9, font.submain = 9, font.legend = 12, font.tickslab = 7,
             font.x = 10, font.y = 10, point.size = 0.3)
  })

  output$schoen_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_genes", "_schoen", ".png", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.pfi()$iter <= 20,
             "Schoenfeld plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.pfi()$coefficients)),
             "Schoenfeld plot not drawn due to infinite coeff(s)!")
      )
      ggsave(file, print(ggcoxzph(cox.zph(react.cox.pfi()), legend.labs = names(react.cox.pfi()$coefficients), 
                            font.main = 9, font.submain = 9, font.legend = 12, font.tickslab = 7,
                            font.x = 10, font.y = 10, point.size = 0.3)), width = 20, units = "cm", dpi = 320)
    }
  )

  output$roc_os <- renderPlot({
    shiny::validate (
      need(react.cox.os()$iter <= 20,
           "ROC plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.os()$coefficients)),
           "ROC plot not drawn due to infinite coeff(s)!")
    )
    data.os <- na.omit(react.cox.os.data())
    coxph.fit <- coxph(react.cox.os()$formula, data = data.os, model = TRUE, x = TRUE, y = TRUE)
    roc.score <- Score(list("COX(coxph.fit$formula[[3]])" = coxph.fit), formula = Surv(OS.time, OS) ~ 1,
                       data = data.os, plots = "ROC", metrics = "AUC")
    plotROC(roc.score)
  })
  
  output$roc_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_genes", "_roc", ".png", sep = "")},
    content = function(file) {
      shiny::validate(
        need(nrow(react.dataset()) != 0,
             "Please choose other dataset, as there is no primary tumor info for this!")
      )
      shiny::validate (
        need(react.cox.os()$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.os()$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
      data.os <- na.omit(react.cox.os.data())
      coxph.fit <- coxph(react.cox.os()$formula, data = data.os, model = TRUE, x = TRUE, y = TRUE)
      roc.score <- Score(list("COX(coxph.fit$formula[[3]])" = coxph.fit), formula = Surv(OS.time, OS) ~ 1,
                         data = data.os, plots = "ROC", metrics = "AUC")
      ggsave(file, print(plotROC(roc.score)), dpi = 320)
    }
  )
  
  output$roc_dss <- renderPlot({
    shiny::validate (
      need(react.cox.dss()$iter <= 20,
           "ROC plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.dss()$coefficients)),
           "ROC plot not drawn due to infinite coeff(s)!")
    )
    data.dss <- na.omit(react.cox.dss.data())
    coxph.fit <- coxph(react.cox.dss()$formula, data = data.dss, model = TRUE, x = TRUE, y = TRUE)
    roc.score <- Score(list("COX(coxph.fit$formula[[3]])" = coxph.fit), formula = Surv(DSS.time, DSS) ~ 1,
                       data = data.dss, plots = "ROC", metrics = "AUC")
    plotROC(roc.score)
  })
  
  output$roc_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_genes", "_roc", ".png", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.dss()$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.dss()$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
      data.dss <- na.omit(react.cox.dss.data())
      coxph.fit <- coxph(react.cox.dss()$formula, data = data.dss, model = TRUE, x = TRUE, y = TRUE)
      roc.score <- Score(list("COX(coxph.fit$formula[[3]])" = coxph.fit), formula = Surv(DSS.time, DSS) ~ 1,
                         data = data.dss, plots = "ROC", metrics = "AUC")
      ggsave(file, print(plotROC(roc.score)), dpi = 320)
    }
  )
  
  output$roc_dfi <- renderPlot({
    shiny::validate (
      need(react.cox.dfi()$iter <= 20,
           "ROC plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.dfi()$coefficients)),
           "ROC plot not drawn due to infinite coeff(s)!")
    )
    data.dfi <- na.omit(react.cox.dfi.data())
    coxph.fit <- coxph(react.cox.dfi()$formula, data = data.dfi, model = TRUE, x = TRUE, y = TRUE)
    roc.score <- Score(list("COX(coxph.fit$formula[[3]])" = coxph.fit), formula = Surv(DFI.time, DFI) ~ 1,
                       data = data.dfi, plots = "ROC", metrics = "AUC")
    plotROC(roc.score)
  })
  
  output$roc_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_genes", "_roc", ".png", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.dfi()$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.dfi()$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
        data.dfi <- na.omit(react.cox.dfi.data())
        coxph.fit <- coxph(react.cox.dfi()$formula, data = data.dfi, model = TRUE, x = TRUE, y = TRUE)
        roc.score <- Score(list("COX(coxph.fit$formula[[3]])" = coxph.fit), formula = Surv(DFI.time, DFI) ~ 1,
                           data = data.dfi, plots = "ROC", metrics = "AUC")
        ggsave(file, print(plotROC(roc.score)), dpi = 320)
    }
  )
  
  output$roc_pfi <- renderPlot({
    shiny::validate (
      need(react.cox.pfi()$iter <= 20,
           "ROC plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.pfi()$coefficients)),
           "ROC plot not drawn due to infinite coeff(s)!")
    )
      data.pfi <- na.omit(react.cox.pfi.data())
      coxph.fit <- coxph(react.cox.pfi()$formula, data = data.pfi, model = TRUE, x = TRUE, y = TRUE)
      roc.score <- Score(list("COX(coxph.fit$formula[[3]])" = coxph.fit), formula = Surv(PFI.time, PFI) ~ 1,
                         data = data.pfi, plots = "ROC", metrics = "AUC")
      plotROC(roc.score)
  })

  output$roc_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_genes", "_roc", ".png", sep = "")},
    content = function(file) {
      shiny::validate(
        need(nrow(react.dataset()) != 0,
             "Please choose other dataset, as there is no primary tumor info for this!")
      )
      shiny::validate (
        need(react.cox.pfi()$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.pfi()$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
      data.pfi <- na.omit(react.cox.pfi.data())
      coxph.fit <- coxph(react.cox.pfi()$formula, data = data.pfi, model = TRUE, x = TRUE, y = TRUE)
      roc.score <- Score(list("COX(coxph.fit$formula[[3]])" = coxph.fit), formula = Surv(PFI.time, PFI) ~ 1,
                         data = data.pfi, plots = "ROC", metrics = "AUC")
      ggsave(file, print(plotROC(roc.score)), dpi = 320)
    }
  )
  
  output$roc_data_os <- renderTable({
    shiny::validate (
      need(react.cox.os()$iter <= 20,
           "ROC plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.os()$coefficients)),
           "ROC plot not drawn due to infinite coeff(s)!")
    )
    data.os <- na.omit(react.cox.os.data())
    coxph.fit <- coxph(react.cox.os()$formula, data = data.os, model = TRUE, x = TRUE, y = TRUE)
    roc.score <- Score(list("COX(coxph.fit$formula[[3]])" = coxph.fit), formula = Surv(OS.time, OS) ~ 1,
                       data = data.os, plots = "ROC", metrics = "AUC")
    roc.data <- cbind(roc.score$AUC$score$AUC*100, roc.score$AUC$score$lower*100, roc.score$AUC$score$upper*100, roc.score$AUC$score$se*100)
    colnames(roc.data) <- c("lower AUC%", "AUC%", "upper AUC%", "standard error%")
    roc.data
  })
  
  output$roc_data_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_genes", "_rocdata", ".tsv", sep = "")},
    content = function(file) {
      shiny::validate(
        need(nrow(react.dataset()) != 0,
             "Please choose other dataset, as there is no primary tumor info for this!")
      )
      shiny::validate (
        need(react.cox.os()$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.os()$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
      data.os <- na.omit(react.cox.os.data())
      coxph.fit <- coxph(react.cox.os()$formula, data = data.os, model = TRUE, x = TRUE, y = TRUE)
      roc.score <- Score(list("COX(coxph.fit$formula[[3]])" = coxph.fit), formula = Surv(OS.time, OS) ~ 1,
                         data = data.os, plots = "ROC", metrics = "AUC")
      roc.data <- cbind(roc.score$AUC$score$AUC*100, roc.score$AUC$score$lower*100, roc.score$AUC$score$upper*100, roc.score$AUC$score$se*100)
      colnames(roc.data) <- c("lower AUC%", "AUC%", "upper AUC%", "standard error%")
      write.table(roc.data, file)}
  )
  
  output$roc_data_dss <- renderTable({
    shiny::validate (
      need(react.cox.dss()$iter <= 20,
           "ROC plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.dss()$coefficients)),
           "ROC plot not drawn due to infinite coeff(s)!")
    )
    data.dss <- na.omit(react.cox.dss.data())
    coxph.fit <- coxph(react.cox.dss()$formula, data = data.dss, model = TRUE, x = TRUE, y = TRUE)
    roc.score <- Score(list("COX(coxph.fit$formula[[3]])" = coxph.fit), formula = Surv(DSS.time, DSS) ~ 1,
                       data = data.dss, plots = "ROC", metrics = "AUC")
    roc.data <- cbind(roc.score$AUC$score$AUC*100, roc.score$AUC$score$lower*100, roc.score$AUC$score$upper*100, roc.score$AUC$score$se*100)
    colnames(roc.data) <- c("lower AUC%", "AUC%", "upper AUC%", "standard error%")
    roc.data
  })
  
  output$roc_data_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_genes", "_rocdata", ".tsv", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.dss()$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.dss()$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
      data.dss <- na.omit(react.cox.dss.data())
      coxph.fit <- coxph(react.cox.dss()$formula, data = data.dss, model = TRUE, x = TRUE, y = TRUE)
      roc.score <- Score(list("COX(coxph.fit$formula[[3]])" = coxph.fit), formula = Surv(DSS.time, DSS) ~ 1,
                         data = data.dss, plots = "ROC", metrics = "AUC")
      roc.data <- cbind(roc.score$AUC$score$AUC*100, roc.score$AUC$score$lower*100, roc.score$AUC$score$upper*100, roc.score$AUC$score$se*100)
      colnames(roc.data) <- c("lower AUC%", "AUC%", "upper AUC%", "standard error%")
      write.table(roc.data, file)}
  )
  
  output$roc_data_dfi <- renderTable({
    shiny::validate (
      need(react.cox.dfi()$iter <= 20,
           "ROC plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.dfi()$coefficients)),
           "ROC plot not drawn due to infinite coeff(s)!")
    )
    data.dfi <- na.omit(react.cox.dfi.data())
    coxph.fit <- coxph(react.cox.dfi()$formula, data = data.dfi, model = TRUE, x = TRUE, y = TRUE)
    roc.score <- Score(list("COX(coxph.fit$formula[[3]])" = coxph.fit), formula = Surv(DFI.time, DFI) ~ 1,
                       data = data.dfi, plots = "ROC", metrics = "AUC")
    roc.data <- cbind(roc.score$AUC$score$AUC*100, roc.score$AUC$score$lower*100, roc.score$AUC$score$upper*100, roc.score$AUC$score$se*100)
    colnames(roc.data) <- c("lower AUC%", "AUC%", "upper AUC%", "standard error%")
    roc.data
  })
  
  output$roc_data_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_genes", "_rocdata", ".tsv", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.dfi()$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.dfi()$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
      data.dfi <- na.omit(react.cox.dfi.data())
      coxph.fit <- coxph(react.cox.dfi()$formula, data = data.dfi, model = TRUE, x = TRUE, y = TRUE)
      roc.score <- Score(list("COX(coxph.fit$formula[[3]])" = coxph.fit), formula = Surv(DFI.time, DFI) ~ 1,
                         data = data.dfi, plots = "ROC", metrics = "AUC")
      roc.data <- cbind(roc.score$AUC$score$AUC*100, roc.score$AUC$score$lower*100, roc.score$AUC$score$upper*100, roc.score$AUC$score$se*100)
      colnames(roc.data) <- c("lower AUC%", "AUC%", "upper AUC%", "standard error%")
      write.table(roc.data, file)}
  )
  
  output$roc_data_pfi <- renderTable({
    shiny::validate (
      need(react.cox.pfi()$iter <= 20,
           "ROC plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.pfi()$coefficients)),
           "ROC plot not drawn due to infinite coeff(s)!")
    )
    data.pfi <- na.omit(react.cox.pfi.data())
    coxph.fit <- coxph(react.cox.pfi()$formula, data = data.pfi, model = TRUE, x = TRUE, y = TRUE)
    roc.score <- Score(list("COX(coxph.fit$formula[[3]])" = coxph.fit), formula = Surv(PFI.time, PFI) ~ 1,
                       data = data.pfi, plots = "ROC", metrics = "AUC")
    roc.data <- cbind(roc.score$AUC$score$AUC*100, roc.score$AUC$score$lower*100, roc.score$AUC$score$upper*100, roc.score$AUC$score$se*100)
    colnames(roc.data) <- c("lower AUC%", "AUC%", "upper AUC%", "standard error%")
    roc.data
  })
  
  output$roc_data_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_genes", "_rocdata", ".tsv", sep = "")},
    content = function(file) {
      shiny::validate(
        need(nrow(react.dataset()) != 0,
             "Please choose other dataset, as there is no primary tumor info for this!")
      )
      shiny::validate (
        need(react.cox.pfi()$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.pfi()$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
      data.pfi <- na.omit(react.cox.pfi.data())
      coxph.fit <- coxph(react.cox.pfi()$formula, data = data.pfi, model = TRUE, x = TRUE, y = TRUE)
      roc.score <- Score(list("COX(coxph.fit$formula[[3]])" = coxph.fit), formula = Surv(PFI.time, PFI) ~ 1,
                         data = data.pfi, plots = "ROC", metrics = "AUC")
      roc.data <- cbind(roc.score$AUC$score$AUC*100, roc.score$AUC$score$lower*100, roc.score$AUC$score$upper*100, roc.score$AUC$score$se*100)
      colnames(roc.data) <- c("lower AUC%", "AUC%", "upper AUC%", "standard error%")
      write.table(roc.data, file)}
  )
  
  
  
  output$anova_data_os <- renderTable({
    shiny::validate (
      need(react.cox.os()$iter <= 20,
           "ANOVA CPH table not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.os()$coefficients)),
           "ANOVA CPH table not drawn due to infinite coeff(s)!")
    )
    as.table(anova(react.cph.os()))
  })
  
  output$anova_data_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_genes", "_cph_anova", ".tsv", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.os()$iter <= 20,
             "CPH ANOVA table not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.os()$coefficients)),
             "CPH ANOVA table not drawn due to infinite coeff(s)!")
      )
      write.table(anova(react.cph.os()), file)
    }
  )
  
  output$anova_data_dss <- renderTable({
    shiny::validate (
      need(react.cox.dss()$iter <= 20,
           "CPH ANOVA table not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.dss()$coefficients)),
           "CPH ANOVA table not drawn due to infinite coeff(s)!")
    )
    as.table(anova(react.cph.dss()))
  })
  
  output$anova_data_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_genes", "_cph_anova", ".tsv", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.dss()$iter <= 20,
             "CPH ANOVA table not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.dss()$coefficients)),
             "CPH ANOVA table not drawn due to infinite coeff(s)!")
      )
      write.table(anova(react.cph.dss()), file)
    }
  )
  
  output$anova_data_dfi <- renderTable({
    shiny::validate (
      need(react.cox.dfi()$iter <= 20,
           "CPH ANOVA table not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.dfi()$coefficients)),
           "CPH ANOVA table not drawn due to infinite coeff(s)!")
    )
    as.table(anova(react.cph.dfi()))
  })
  
  output$anova_data_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_genes", "_cph_anova", ".tsv", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.dfi()$iter <= 20,
             "CPH ANOVA table not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.dfi()$coefficients)),
             "CPH ANOVA table not drawn due to infinite coeff(s)!")
      )
      write.table(aanova(react.cph.dfi()), file)
    }
  )
  
  output$anova_data_pfi <- renderTable({
    shiny::validate (
      need(react.cox.pfi()$iter <= 20,
           "CPH ANOVA table not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.pfi()$coefficients)),
           "CPH ANOVA table not drawn due to infinite coeff(s)!")
    )
    as.table(anova(react.cph.pfi()))
  })
  
  output$anova_data_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_genes", "_cph_anova", ".tsv", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.pfi()$iter <= 20,
             "CPH ANOVA table not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.pfi()$coefficients)),
             "CPH ANOVA table not drawn due to infinite coeff(s)!")
      )
      write.table(anova(react.cph.pfi()), file)
    }
  )
  
  output$anova_plot_os <- renderPlot({
    shiny::validate (
      need(react.cox.os()$iter <= 20,
           "CPH ANOVA plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.os()$coefficients)),
           "CPH ANOVA plot not drawn due to infinite coeff(s)!")
    )
    anova.cph <- as.data.frame(anova(react.cph.os()))[1:length(input$genes), ]
    ggplot(anova.cph, aes(x = anova.cph[, 1] - anova.cph[, 2], y = anova.cph[, 3], label = rownames(anova.cph))) + 
      geom_point(size = 0.5, col = "darkgreen") + labs(x = expression(paste("\U03C7"^2, "- df")), y = "p-value") + 
      geom_label_repel(fill = "white", size = 3)
  })
  
  output$anova_plot_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_genes", "_cph_anova", ".png", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.os()$iter <= 20,
             "CPH ANOVA plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.os()$coefficients)),
             "CPH ANOVA plot not drawn due to infinite coeff(s)!")
      )
      anova.cph <- as.data.frame(anova(react.cph.os()))[1:length(input$genes), ]
      ggsave(file, ggplot(anova.cph, aes(x = anova.cph[, 1] - anova.cph[, 2], y = anova.cph[, 3], label = rownames(anova.cph))) + 
               geom_point(size = 0.5, col = "darkgreen") + labs(x = expression(paste("\U03C7"^2, "- df")), y = "p-value") + 
               geom_label_repel(fill = "white", size = 3), dpi = 320)
    }
  )
  
  output$anova_plot_dss <- renderPlot({
    shiny::validate (
      need(react.cox.dss()$iter <= 20,
           "CPH ANOVA plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.dss()$coefficients)),
           "CPH ANOVA plot not drawn due to infinite coeff(s)!")
    )
    anova.cph <- as.data.frame(anova(react.cph.dss()))[1:length(input$genes), ]
    ggplot(anova.cph, aes(x = anova.cph[, 1] - anova.cph[, 2], y = anova.cph[, 3], label = rownames(anova.cph))) + 
      geom_point(size = 0.5, col = "darkgreen") + labs(x = expression(paste("\U03C7"^2, "- df")), y = "p-value") + 
      geom_label_repel(fill = "white", size = 3)
    
  })
  
  output$anova_plot_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_genes", "_cph_anova", ".png", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.dss()$iter <= 20,
             "CPH ANOVA plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.dss()$coefficients)),
             "CPH ANOVA plot not drawn due to infinite coeff(s)!")
      )
      anova.cph <- as.data.frame(anova(react.cph.dss()))[1:length(input$genes), ]
      ggsave(file, ggplot(anova.cph, aes(x = anova.cph[, 1] - anova.cph[, 2], y = anova.cph[, 3], label = rownames(anova.cph))) + 
               geom_point(size = 0.5, col = "darkgreen") + labs(x = expression(paste("\U03C7"^2, "- df")), y = "p-value") + 
               geom_label_repel(fill = "white", size = 3), dpi = 320)
    }
  )
  
  output$anova_plot_dfi <- renderPlot({
    shiny::validate (
      need(react.cox.dfi()$iter <= 20,
           "CPH ANOVA plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.dfi()$coefficients)),
           "CPH ANOVA plot not drawn due to infinite coeff(s)!")
    )
    anova.cph <- as.data.frame(anova(react.cph.dfi()))[1:length(input$genes), ]
    ggplot(anova.cph, aes(x = anova.cph[, 1] - anova.cph[, 2], y = anova.cph[, 3], label = rownames(anova.cph))) + 
      geom_point(size = 0.5, col = "darkgreen") + labs(x = expression(paste("\U03C7"^2, "- df")), y = "p-value") + 
      geom_label_repel(fill = "white",size = 3)
  })
  
  output$anova_plot_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_genes", "_cph_anova", ".png", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.dfi()$iter <= 20,
             "CPH ANOVA plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.dfi()$coefficients)),
             "CPH ANOVA plot not drawn due to infinite coeff(s)!")
      )
      anova.cph <- as.data.frame(anova(react.cph.dfi()))[1:length(input$genes), ]
      ggsave(file, ggplot(anova.cph, aes(x = anova.cph[, 1] - anova.cph[, 2], y = anova.cph[, 3], label = rownames(anova.cph))) + 
               geom_point(size = 0.5, col = "darkgreen") + labs(x = expression(paste("\U03C7"^2, "- df")), y = "p-value") + 
               geom_label_repel(fill = "white", size = 3), dpi = 320)
    }
  )
  
  output$anova_plot_pfi <- renderPlot({
    shiny::validate (
      need(react.cox.pfi()$iter <= 20,
           "CPH ANOVA plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.pfi()$coefficients)),
           "CPH ANOVA plot not drawn due to infinite coeff(s)!")
    )
    anova.cph <- as.data.frame(anova(react.cph.pfi()))[1:length(input$genes), ]
    ggplot(anova.cph, aes(x = anova.cph[, 1] - anova.cph[, 2], y = anova.cph[, 3], label = rownames(anova.cph))) + 
      geom_point(size = 0.5, col = "darkgreen") + labs(x = expression(paste("\U03C7"^2, "- df")), y = "p-value") + 
      geom_label_repel(fill = "white", size = 3)
  })
  
  output$anova_plot_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_genes", "_cph_anova", ".png", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.pfi()$iter <= 20,
             "CPH ANOVA plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.pfi()$coefficients)),
             "CPH ANOVA plot not drawn due to infinite coeff(s)!")
      )
      anova.cph <- as.data.frame(anova(react.cph.pfi()))[1:length(input$genes), ]
      ggsave(file, ggplot(anova.cph, aes(x = anova.cph[, 1] - anova.cph[, 2], y = anova.cph[, 3], label = rownames(anova.cph))) + 
               geom_point(size = 0.5, col = "darkgreen") + labs(x = expression(paste("\U03C7"^2, "- df")), y = "p-value") + 
               geom_label_repel(fill = "white", size = 3), dpi = 320)
    }
  )
  
  output$valid_data_os <- renderTable({
    shiny::validate (
      need(react.cox.os()$iter <= 20,
           "CPH validation table not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.os()$coefficients)),
           "CPH validation table not drawn due to infinite coeff(s)!")
    )
    valid <- rms::validate(react.cph.os())
    shiny::validate(
      need(!all(is.na(valid[,2])), "NA results in training set!")
    )
    as.table(valid)
  })
  
  output$valid_data_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_genes", "_cph_validate", ".tsv", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.os()$iter <= 20,
             "CPH validation table not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.os()$coefficients)),
             "CPH validation table not drawn due to infinite coeff(s)!")
      )
      valid <- rms::validate(react.cph.os())
      shiny::validate(
        need(!all(is.na(valid[,2])), "NA results in training set!")
      )
      write.table(valid, file)
    }
  )
  
  output$valid_data_dss <- renderTable({
    shiny::validate (
      need(react.cox.dss()$iter <= 20,
           "CPH validation table not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.dss()$coefficients)),
           "CPH validation table not drawn due to infinite coeff(s)!")
    )
    valid <- rms::validate(react.cph.dss())
    shiny::validate(
      need(!all(is.na(valid[,2])), "NA results in training set!")
    )
    as.table(valid)
  })
  
  output$valid_data_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_genes", "_cph_validate", ".tsv", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.dss()$iter <= 20,
             "CPH validation table not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.dss()$coefficients)),
             "CPH validation table not drawn due to infinite coeff(s)!")
      )
      valid <- rms::validate(react.cph.dss())
      shiny::validate(
        need(!all(is.na(valid[,2])), "NA results in training set!")
      )
      write.table(valid, file)
    }
  )
  
  output$valid_data_dfi <- renderTable({
    shiny::validate (
      need(react.cox.dfi()$iter <= 20,
           "CPH validation table not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.dfi()$coefficients)),
           "CPH validation table not drawn due to infinite coeff(s)!")
    )
    valid <- rms::validate(react.cph.dfi())
    shiny::validate(
      need(!all(is.na(valid[,2])), "NA results in training set!")
    )
    as.table(valid)
  })
  
  output$valid_data_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_genes", "_cph_validate", ".tsv", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.dfi()$iter <= 20,
             "CPH validation table not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.dfi()$coefficients)),
             "CPH validation table not drawn due to infinite coeff(s)!")
      )
      valid <- rms::validate(react.cph.dfi())
      shiny::validate(
        need(!all(is.na(valid[,2])), "NA results in training set!")
      )
      write.table(valid, file)
    }
  )
  
  output$valid_data_pfi <- renderTable({
    shiny::validate (
      need(react.cox.pfi()$iter <= 20,
           "CPH validation table not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.pfi()$coefficients)),
           "CPH validation table not drawn due to infinite coeff(s)!")
    )
    valid <- rms::validate(react.cph.pfi())
    shiny::validate(
      need(!all(is.na(valid[,2])), "NA results in training set!")
    )
    as.table(valid)
  })
  
  output$valid_data_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_genes", "_cph_validate", ".tsv", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.pfi()$iter <= 20,
             "CPH validation table not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.pfi()$coefficients)),
             "CPH validation table not drawn due to infinite coeff(s)!")
      )
      valid <- rms::validate(react.cph.pfi())
      shiny::validate(
        need(!all(is.na(valid[,2])), "NA results in training set!")
      )
      write.table(valid, file)
    }
  )
  
  output$valid_plot_os <- renderPlot({
    shiny::validate (
      need(react.cox.os()$iter <= 20,
           "CPH validation plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.os()$coefficients)),
           "CPH validation plot not drawn due to infinite coeff(s)!")
    )
    valid <- rms::validate(react.cph.os())
    shiny::validate(
      need(!all(is.na(valid[,2])), "NA results in training set!")
    )
    plot(valid, type = "p", col = "darkblue", pch = 18)
  })
  
  output$valid_plot_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_genes", "_cph_validate", ".png", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.os()$iter <= 20,
             "CPH validation plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.os()$coefficients)),
             "CPH validation plot not drawn due to infinite coeff(s)!")
      )
      valid <- rms::validate(react.cph.os())
      shiny::validate(
        need(!all(is.na(valid[,2])), "NA results in training set!")
      )
      ggsave(file, print(plot(valid, type = "p", col = "darkblue", pch = 18)), dpi = 320)
    }
  )
  
  output$valid_plot_dss <- renderPlot({
    shiny::validate (
      need(react.cox.dss()$iter <= 20,
           "CPH validation plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.dss()$coefficients)),
           "CPH validation plot not drawn due to infinite coeff(s)!")
    )
    valid <- rms::validate(react.cph.dss())
    shiny::validate(
      need(!all(is.na(valid[,2])), "NA results in training set!")
    )
    plot(valid, type = "p", col = "darkblue", pch = 18)
  })
  
  output$valid_plot_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_genes", "_cph_validate", ".png", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.dss()$iter <= 20,
             "CPH validation plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.dss()$coefficients)),
             "CPH validation plot not drawn due to infinite coeff(s)!")
      )
      valid <- rms::validate(react.cph.dss())
      shiny::validate(
        need(!all(is.na(valid[,2])), "NA results in training set!")
      )
      ggsave(file, print(plot(valid, type = "p", col = "darkblue", pch = 18)), dpi = 320)
    }
  )
  
  output$valid_plot_dfi <- renderPlot({
    shiny::validate (
      need(react.cox.dfi()$iter <= 20,
           "CPH validation plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.dfi()$coefficients)),
           "CPH validation plot not drawn due to infinite coeff(s)!")
    )
    valid <- rms::validate(react.cph.dfi())
    shiny::validate(
      need(!all(is.na(valid[,2])), "NA results in training set!")
    )
    plot(valid, type = "p", col = "darkblue", pch = 18)
  })
  
  output$valid_plot_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_genes", "_cph_validate", ".png", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.dfi()$iter <= 20,
             "CPH validation plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.dfi()$coefficients)),
             "CPH validation plot not drawn due to infinite coeff(s)!")
      )
      valid <- rms::validate(react.cph.dfi())
      shiny::validate(
        need(!all(is.na(valid[,2])), "NA results in training set!")
      )
      ggsave(file, print(plot(valid, type = "p", col = "darkblue", pch = 18)), dpi = 320)
    }
  )
  
  output$valid_plot_pfi <- renderPlot({
    shiny::validate (
      need(react.cox.pfi()$iter <= 20,
           "CPH validation plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(react.cox.pfi()$coefficients)),
           "CPH validation plot not drawn due to infinite coeff(s)!")
    )
    valid <- rms::validate(react.cph.pfi())
    shiny::validate(
      need(!all(is.na(valid[,2])), "NA results in training set!")
    )
    plot(valid, type = "p", col = "darkblue", pch = 18, cex = 1)
  })
  
  output$valid_plot_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_genes", "_cph_validate", ".png", sep = "")},
    content = function(file) {
      shiny::validate (
        need(react.cox.pfi()$iter <= 20,
             "CPH validation plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(react.cox.pfi()$coefficients)),
             "CPH validation plot not drawn due to infinite coeff(s)!")
      )
      valid <- rms::validate(react.cph.pfi())
      shiny::validate(
        need(!all(is.na(valid[,2])), "NA results in training set!")
      )
      ggsave(file, print(plot(valid, type = "p", col = "darkblue", pch = 18)), dpi = 320)
    }
  )

###########################################################
  output$glm_pred_os <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.os.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    glm <- glmnet.os(data, input$n_folds, input$alph)[[1]]
    pr <- as.data.frame(glm)
    colnames(pr) <- "coef"
    shiny::validate(
      need(!is.na(pr$coef),
           "No gene (coefficient) found in model!")
    )
    ggplot(pr, aes(rownames(pr), coef)) + geom_point(col = "darkred", size = 0.5) + 
      xlab("genes") + ylab("coefficients")
  })

  output$glm_pred_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_glm_coeffs", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.os.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      glm <- glmnet.os(data, input$n_folds, input$alph)[[1]]
      pr <- as.data.frame(glm)
      colnames(pr) <- "coef"
      shiny::validate(
        need(!is.na(pr$coef),
             "No gene (coefficient) found in model!")
      )
      ggsave(file, ggplot(pr, aes(rownames(pr), coef)) + geom_point(col = "darkred", size = 0.5) + 
               xlab("genes") + ylab("coefficients"), dpi = 320)
    }
  )

  output$glm_pred_dss <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dss.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
           "Insufficient data for disease-specific survival")
    )
    glm <- glmnet.dss(data, input$n_folds, input$alph)[[1]]
    pr <- as.data.frame(glm)
    colnames(pr) <- "coef"
    shiny::validate(
      need(!is.na(pr$coef),
           "No gene (coefficient) found in model!")
    )
    ggplot(pr, aes(rownames(pr), coef)) + geom_point(col = "darkred", size = 0.5) + xlab("genes") + ylab("coefficients")
  })

  output$glm_pred_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_glm_coeffs", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dss.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
             "Insufficient data for disease-specific survival")
      )
      glm <- glmnet.dss(data, input$n_folds, input$alph)[[1]]
      pr <- as.data.frame(glm)
      colnames(pr) <- "coef"
      shiny::validate(
      need(!is.na(pr$coef),
           "No gene (coefficient) found in model!")
      )
      ggsave(file, ggplot(pr, aes(rownames(pr), coef)) + geom_point(col = "darkred", size = 0.5) + 
               xlab("genes") + ylab("coefficients"), dpi = 320)
    }
  )

  output$glm_pred_dfi <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    glm <- glmnet.dfi(data, input$n_folds, input$alph)[[1]]
    pr <- as.data.frame(glm)
    colnames(pr) <- "coef"
    shiny::validate(
      need(!is.na(pr$coef),
           "No gene (coefficient) found in model!")
    )
    ggplot(pr, aes(rownames(pr), coef)) + geom_point(col = "darkred", size = 0.5) + xlab("genes") + ylab("coefficients")
  })

  output$glm_pred_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_glm_coeffs", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      glm <- glmnet.dfi(data, input$n_folds, input$alph)[[1]]
      pr <- as.data.frame(glm)
      colnames(pr) <- "coef"
      shiny::validate(
        need(!is.na(pr$coef),
             "No gene (coefficient) found in model!")
      )
      ggsave(file, ggplot(pr, aes(rownames(pr), coef)) + geom_point(col = "darkred", size = 0.5) + 
               xlab("genes") + ylab("coefficients"), dpi = 320)
    }
  )

  output$glm_pred_pfi <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.pfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
           "Insufficient data for progression-free survival")
    )
    glm <- glmnet.pfi(data, input$n_folds, input$alph)[[1]]
    pr <- as.data.frame(glm)
    colnames(pr) <- "coef"
    shiny::validate(
      need(!is.na(pr$coef),
           "No gene (coefficient) found in model!")
    )
    ggplot(pr, aes(rownames(pr), coef)) + geom_point(col = "darkred", size = 0.5) + xlab("genes") + ylab("coefficients")
  })

  output$glm_pred_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_glm_coeffs", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.pfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
             "Insufficient data for progression-free survival")
      )
      glm <- glmnet.pfi(data, input$n_folds, input$alph)[[1]]
      pr <- as.data.frame(glm)
      colnames(pr) <- "coef"
      shiny::validate(
        need(!is.na(pr$coef),
             "No gene (coefficient) found in model!")
      )
      ggsave(file, ggplot(pr, aes(rownames(pr), coef)) + geom_point(col = "darkred", size = 0.5) + 
               xlab("genes") + ylab("coefficients"), dpi = 320)
    }
  )

  output$glm_coef_tab_os <- renderTable({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.os.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    glm <- glmnet.os(data, input$n_folds, input$alph)[[1]]
    pr <- as.data.frame(glm)
    colnames(pr) <- "coef"
    shiny::validate(
      need(!is.na(pr$coef),
           "No gene (coefficient) found in model!")
    )
    pr <- cbind(rownames(pr), pr)
    colnames(pr)[1] <- "gene"
    pr
  })
  
  output$glm_coef_tab_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_glm_coeffs", ".tsv", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.os.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      glm <- glmnet.os(data, input$n_folds, input$alph)[[1]]
      pr <- as.data.frame(glm)
      colnames(pr) <- "coef"
      shiny::validate(
        need(!is.na(pr$coef),
             "No gene (coefficient) found in model!")
      )
      pr <- cbind(rownames(pr), pr)
      colnames(pr)[1] <- "gene"
      write.table(pr, file)
    }
  )
  
  output$glm_coef_tab_dss <- renderTable({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dss.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
           "Insufficient data for disease-specific survival")
    )
    glm <- glmnet.dss(data, input$n_folds, input$alph)[[1]]
    pr <- as.data.frame(glm)
    colnames(pr) <- "coef"
    shiny::validate(
      need(!is.na(pr$coef),
           "No gene (coefficient) found in model!")
    )
    pr <- cbind(rownames(pr), pr)
    colnames(pr)[1] <- "gene"
    pr
  })
  
  output$glm_coef_tab_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_glm_coeffs", ".tsv", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dss.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
             "Insufficient data for disease-specific survival")
      )
      glm <- glmnet.dss(data, input$n_folds, input$alph)[[1]]
      pr <- as.data.frame(glm)
      colnames(pr) <- "coef"
      shiny::validate(
        need(!is.na(pr$coef),
             "No gene (coefficient) found in model!")
      )
      pr <- cbind(rownames(pr), pr)
      colnames(pr)[1] <- "gene"
      write.table(pr, file)
    }
  )
  
  output$glm_coef_tab_dfi <- renderTable({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    glm <- glmnet.dfi(data, input$n_folds, input$alph)[[1]]
    pr <- as.data.frame(glm)
    colnames(pr) <- "coef"
    shiny::validate(
      need(!is.na(pr$coef),
           "No gene (coefficient) found in model!")
    )
    pr <- cbind(rownames(pr), pr)
    colnames(pr)[1] <- "gene"
    pr
  })
  
  output$glm_coef_tab_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_glm_coeffs", ".tsv", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      glm <- glmnet.dfi(data, input$n_folds, input$alph)[[1]]
      pr <- as.data.frame(glm)
      colnames(pr) <- "coef"
      shiny::validate(
        need(!is.na(pr$coef),
             "No gene (coefficient) found in model!")
      )
      pr <- cbind(rownames(pr), pr)
      colnames(pr)[1] <- "gene"
      write.table(pr, file)
    }
  )
  
  output$glm_coef_tab_pfi <- renderTable({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.pfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
           "Insufficient data for progression-free survival")
    )
    glm <- glmnet.pfi(data, input$n_folds, input$alph)[[1]]
    pr <- as.data.frame(glm)
    colnames(pr) <- "coef"
    shiny::validate(
      need(!is.na(pr$coef),
           "No gene (coefficient) found in model!")
    )
    pr <- cbind(rownames(pr), pr)
    colnames(pr)[1] <- "gene"
    pr
  })
  
  output$glm_coef_tab_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_glm_coeffs", ".tsv", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.pfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
             "Insufficient data for progression-free survival")
      )
      glm <- glmnet.pfi(data, input$n_folds, input$alph)[[1]]
      pr <- as.data.frame(glm)
      colnames(pr) <- "coef"
      shiny::validate(
        need(!is.na(pr$coef),
             "No gene (coefficient) found in model!")
      )
      pr <- cbind(rownames(pr), pr)
      colnames(pr)[1] <- "gene"
      write.table(pr, file)
    }
  )
  
  output$glm_km_os <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.os.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    glm.data <- glmnet.os(data, input$n_folds, input$alph)[[3]]
    shiny::validate(
      need(!is.null(glm.data),
           "No gene (coefficient) found in model!")
    )
    surv.fit.os <- survfit(Surv(OS.time, OS) ~ Risk, data = glm.data)
    ggsurvplot(surv.fit.os, data = glm.data, palette = c("#FF9E29", "#86AA00"), pval = TRUE,
               risk.table = TRUE, conf.int = TRUE, risk.table.col = "strata") + xlab("days")
  })

  output$glm_km_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_glmPI", "_km", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.os.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      glm.data <- glmnet.os(data, input$n_folds, input$alph)[[3]]
      shiny::validate(
        need(!is.na(glm.data),
             "No gene (coefficient) found in model!")
      )
      surv.fit.os <- survfit(Surv(OS.time, OS) ~ Risk, data = glm.data)
      ggsave(file, print(ggsurvplot(surv.fit.os, data = glm.data, palette = c("#FF9E29", "#86AA00"), pval = TRUE,
                                    risk.table = TRUE, conf.int = TRUE, risk.table.col = "strata") + xlab("days")), dpi = 320)
    }
  )

  output$glm_km_dss <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dss.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    glm.data <- glmnet.dss(data, input$n_folds, input$alph)[[3]]
    shiny::validate(
      need(!is.null(glm.data),
           "No gene (coefficient) found in model!")
    )
    surv.fit.dss <- survfit(Surv(DSS.time, DSS) ~ Risk, data = glm.data)
    ggsurvplot(surv.fit.dss, data = glm.data, palette = c("#FF9E29", "#86AA00"), pval = TRUE,
               risk.table = TRUE, conf.int = TRUE, risk.table.col = "strata") + xlab("days")
  })

  output$glm_km_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_glmPI", "_km", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dss.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      glm.data <- glmnet.dss(data, input$n_folds, input$alph)[[3]]
      shiny::validate(
        need(!is.null(glm.data),
             "No gene (coefficient) found in model!")
      )
      surv.fit.dss <- survfit(Surv(DSS.time, DSS) ~ Risk, data = glm.data)
      ggsave(file, print(ggsurvplot(surv.fit.dss, data = glm.data, palette = c("#FF9E29", "#86AA00"), pval = TRUE,
                                    risk.table = TRUE, conf.int = TRUE, risk.table.col = "strata") + xlab("days")), dpi = 320)
    }
  )

  output$glm_km_dfi <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
          "Insufficient data for disease-free survival!")
    )
    glm.data <- glmnet.dfi(data, input$n_folds, input$alph)[[3]]
    shiny::validate(
      need(!is.null(glm.data),
           "No gene (coefficient) found in model!")
    )
    surv.fit.dfi <- survfit(Surv(DFI.time, DFI) ~ Risk, data = glm.data)
    ggsurvplot(surv.fit.dfi, data = glm.data, palette = c("#FF9E29", "#86AA00"), pval = TRUE,
               risk.table = TRUE, conf.int = TRUE, risk.table.col = "strata") + xlab("days")
  })

  output$glm_km_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_glmPI", "_km", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
             "Insufficient data for disease-free survival!")
      )
      glm.data <- glmnet.dfi(data, input$n_folds, input$alph)[[3]]
      shiny::validate(
        need(!is.null(glm.data),
             "No gene (coefficient) found in model!")
      )
      surv.fit.dfi <- survfit(Surv(DFI.time, DFI) ~ Risk, data = glm.data)
      ggsave(file, print(ggsurvplot(surv.fit.dfi, data = glm.data, palette = c("#FF9E29", "#86AA00"), pval = TRUE,
                                    risk.table = TRUE, conf.int = TRUE, risk.table.col = "strata") + xlab("days")), dpi = 320)
    }
  )

  output$glm_km_pfi <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.pfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
           "Insufficient data for disease-free survival!")
    )
    glm.data <- glmnet.pfi(data, input$n_folds, input$alph)[[3]]
    shiny::validate(
      need(!is.null(glm.data),
           "No gene (coefficient) found in model!")
    )
    surv.fit.pfi <- survfit(Surv(PFI.time, PFI) ~ Risk, data = glm.data)
    ggsurvplot(surv.fit.pfi, data = glm.data, palette = c("#FF9E29", "#86AA00"), pval = TRUE,
               risk.table = TRUE, conf.int = TRUE, risk.table.col = "strata") + xlab("days")
  })

  output$glm_km_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_glmPI", "_km", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.pfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
             "Insufficient data for disease-free survival!")
      )
      glm.data <- glmnet.pfi(data, input$n_folds, input$alph)[[3]]
      shiny::validate(
        need(!is.null(glm.data),
             "No gene (coefficient) found in model!")
      )
      surv.fit.pfi <- survfit(Surv(PFI.time, PFI) ~ Risk, data = glm.data)
      ggsave(file, print(ggsurvplot(surv.fit.pfi, data = glm.data, palette = c("#FF9E29", "#86AA00"), pval = TRUE,
                                    risk.table = TRUE, conf.int = TRUE, risk.table.col = "strata") + xlab("days")), dpi = 320)
    }
  )

  output$glm_data_os <- renderTable({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.os.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    glm.data <- glmnet.os(data, input$n_folds, input$alph)[[3]]
    shiny::validate(
      need(!is.null(glm.data),
           "No gene (coefficient) found in model!")
    )
    glm.data
  })

  output$glm_data_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_glmPI", "_data", ".tsv", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.os.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      glm.data <- glmnet.os(data, input$n_folds, input$alph)[[3]]
      shiny::validate(
        need(!is.null(glm.data),
             "No gene (coefficient) found in model!")
      )
      glm.data$prog_ind <- glm.data$prog_ind[[1]]
      write.table(glm.data, file)}
  )

  output$glm_data_dss <- renderTable({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dss.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    glm.data <- glmnet.dss(data, input$n_folds, input$alph)[[3]]
    shiny::validate(
      need(!is.null(glm.data),
           "No gene (coefficient) found in model!")
    )
    glm.data
  })

  output$glm_data_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_glmPI", "_data", ".tsv", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dss.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      glm.data <- glmnet.dss(data, input$n_folds, input$alph)[[3]]
      shiny::validate(
        need(!is.null(glm.data),
             "No gene (coefficient) found in model!")
      )
      glm.data$prog_ind <- glm.data$prog_ind[[1]]
      write.table(glm.data, file)}
  )

  output$glm_data_dfi <- renderTable({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    glm.data <- glmnet.dfi(data, input$n_folds, input$alph)[[3]]
    shiny::validate(
      need(!is.null(glm.data),
           "No gene (coefficient) found in model!")
    )
    glm.data
  })

  output$glm_data_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_glmPI", "_data", ".tsv", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      glm.data <- glmnet.dfi(data, input$n_folds, input$alph)[[3]]
      shiny::validate(
        need(!is.null(glm.data),
             "No gene (coefficient) found in model!")
      )
      glm.data$prog_ind <- glm.data$prog_ind[[1]]
      write.table(glm.data, file)}
  )

  output$glm_data_pfi <- renderTable({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.pfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    glm.data <- glmnet.pfi(data, input$n_folds, input$alph)[[3]]
    shiny::validate(
      need(!is.null(glm.data),
           "No gene (coefficient) found in model!")
    )
    glm.data
  })

  output$glm_data_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_glmPI", "_data", ".tsv", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.pfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
             "Insufficient data for progress-free survival")
      )
      glm.data <- glmnet.pfi(data, input$n_folds, input$alph)[[3]]
      shiny::validate(
        need(!is.null(glm.data),
             "No gene (coefficient) found in model!")
      )
      glm.data$prog_ind <- glm.data$prog_ind[[1]]
      write.table(glm.data, file)}
  )

  output$glm_cumhz_os <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.os.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    glm.data <- glmnet.os(data, input$n_folds, input$alph)[[3]]
    shiny::validate(
      need(!is.null(glm.data),
           "No gene (coefficient) found in model!")
    )
    surv.fit.os <- survfit(Surv(OS.time, OS) ~ Risk, data = glm.data)
    ggsurvplot(surv.fit.os, data = glm.data, palette = c("#FF007F", "#3399FF"), pval = TRUE, fun = "cumhaz",
               risk.table = TRUE, conf.int = TRUE, risk.table.col = "strata") + xlab("days")
  })

  output$glm_cumhz_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_glmPI", "_cumhaz", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.os.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      glm.data <- glmnet.os(data, input$n_folds, input$alph)[[3]]
      shiny::validate(
        need(!is.null(glm.data),
             "No gene (coefficient) found in model!")
      )
      surv.fit.os <- survfit(Surv(OS.time, OS) ~ Risk, data = glm.data)
      ggsave(file, print(ggsurvplot(surv.fit.os, data = glm.data, palette = c("#FF007F", "#3399FF"), pval = TRUE, fun = "cumhaz",
                                    risk.table = TRUE, conf.int = TRUE, risk.table.col = "strata") + xlab("days")),
                                      dpi = 320)
    }
  )

  output$glm_cumhz_dss <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dss.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    glm.data <- glmnet.dss(data, input$n_folds, input$alph)[[3]]
    shiny::validate(
      need(!is.null(glm.data),
           "No gene (coefficient) found in model!")
    )
    surv.fit.dss <- survfit(Surv(DSS.time, DSS) ~ Risk, data = glm.data)
    ggsurvplot(surv.fit.dss, data = glm.data, palette = c("#FF007F", "#3399FF"), pval = TRUE, fun = "cumhaz",
               risk.table = TRUE, conf.int = TRUE, risk.table.col = "strata") + xlab("days")
  })

  output$glm_cumhz_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_glmPI", "_cumhz", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dss.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      glm.data <- glmnet.dss(data, input$n_folds, input$alph)[[3]]
      shiny::validate(
        need(!is.null(glm.data),
             "No gene (coefficient) found in model!")
      )
      surv.fit.dss <- survfit(Surv(DSS.time, DSS) ~ Risk, data = glm.data)
      ggsave(file, print(ggsurvplot(surv.fit.dss, data = glm.data, palette = c("#FF007F", "#3399FF"), pval = TRUE, fun = "cumhaz",
                                    risk.table = TRUE, conf.int = TRUE, risk.table.col = "strata") + xlab("days")),
                                      dpi = 320)
    }
  )

  output$glm_cumhz_dfi <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    glm.data <- glmnet.dfi(data, input$n_folds, input$alph)[[3]]
    shiny::validate(
      need(!is.null(glm.data),
           "No gene (coefficient) found in model!")
    )
    surv.fit.dfi <- survfit(Surv(DFI.time, DFI) ~ Risk, data = glm.data)
    ggsurvplot(surv.fit.dfi, data = glm.data, palette = c("#FF007F", "#3399FF"), pval = TRUE, fun = "cumhaz",
               risk.table = TRUE, conf.int = TRUE, risk.table.col = "strata") + xlab("days")
  })

  output$glm_cumhz_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_glmPI", "_cumhz", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      glm.data <- glmnet.dfi(data, input$n_folds, input$alph)[[3]]
      shiny::validate(
        need(!is.null(glm.data),
             "No gene (coefficient) found in model!")
      )
      surv.fit.dfi <- survfit(Surv(DFI.time, DFI) ~ Risk, data = glm.data)
      ggsave(file, print(ggsurvplot(surv.fit.dfi, data = glm.data, palette = c("#FF007F", "#3399FF"), pval = TRUE, fun = "cumhaz",
                                    risk.table = TRUE, conf.int = TRUE, risk.table.col = "strata") + xlab("days")),
                                      dpi = 320)
    }
  )

  output$glm_cumhz_pfi <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.pfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    glm.data <- glmnet.pfi(data, input$n_folds, input$alph)[[3]]
    shiny::validate(
      need(!is.null(glm.data),
           "No gene (coefficient) found in model!")
    )
    surv.fit.pfi <- survfit(Surv(PFI.time, PFI) ~ Risk, data = glm.data)
    ggsurvplot(surv.fit.pfi, data = glm.data, palette = c("#FF007F", "#3399FF"), pval = TRUE, fun = "cumhaz",
               risk.table = TRUE, conf.int = TRUE, risk.table.col = "strata") + xlab("days")
  })

  output$glm_cumhz_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_glmPI", "_cumhz", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.pfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      glm.data <- glmnet.pfi(data, input$n_folds, input$alph)[[3]]
      shiny::validate(
        need(!is.null(glm.data),
             "No gene (coefficient) found in model!")
      )
      surv.fit.pfi <- survfit(Surv(PFI.time, PFI) ~ Risk, data = glm.data)
      ggsave(file, print(ggsurvplot(surv.fit.pfi, data = glm.data, palette = c("#FF007F", "#3399FF"), pval = TRUE,
                                    fun = "cumhaz", risk.table = TRUE, conf.int = TRUE, risk.table.col = "strata") + xlab("days")),
                                      dpi = 320)
    }
  )
  
  output$glm_cv_os <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.os.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    cv.glm <- glmnet.os(data, input$n_folds, input$alph)[[2]]
    plot(cv.glm, xlab = expression(paste(log, "\U03BB")))
  })

  output$glm_cv_os_down <- downloadHandler(
     filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_glm", "_cv", ".png", sep = "")},
     content = function(file) {
       req(length(input$genes) >= 2)
       data <- isolate(react.cox.os.data())
       shiny::validate(
         need(input$genes %in% colnames(data),
              "Selected genes not found in input dataset!")
       )
       shiny::validate(
         need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
              "Insufficient data for disease-free survival")
       )
       cv.glm <- glmnet.os(data, input$n_folds, input$alph)[[2]]
       ggsave(file, print(plot(cv.glm, xlab = expression(paste(log, "\U03BB")))), dpi = 320)
     }
  )

  output$glm_cv_dss <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dss.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    cv.glm <- glmnet.dss(data, input$n_folds, input$alph)[[2]]
    plot(cv.glm, xlab = expression(paste(log, "\U03BB")))
  })

  output$glm_cv_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_glm", "_cv", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dss.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      cv.glm <- glmnet.dss(data, input$n_folds, input$alph)[[2]]
      ggsave(file, print(plot(cv.glm, xlab = expression(paste(log, "\U03BB")))), dpi = 320)
    }
  )

  output$glm_cv_dfi <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    cv.glm <- glmnet.dfi(data, input$n_folds, input$alph)[[2]]
    plot(cv.glm, xlab = expression(paste(log, "\U03BB")))
  })

  output$glm_cv_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_glm", "_cv", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      cv.glm <- glmnet.dfi(data, input$n_folds, input$alph)[[2]]
      ggsave(file, print(plot(cv.glm, xlab = expression(paste(log, "\U03BB")))), dpi = 320)
    }
  )

  output$glm_cv_pfi <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.pfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    cv.glm <- glmnet.pfi(data, input$n_folds, input$alph)[[2]]
    plot(cv.glm, xlab = expression(paste(log, "\U03BB")))
  })

  output$glm_cv_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_glm", "_cv", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.pfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      cv.glm <- glmnet.pfi(data, input$n_folds, input$alph)[[2]]
      ggsave(file, print(plot(cv.glm, xlab = expression(paste(log, "\U03BB")))), dpi = 320)
    }
  )

  output$glm_roc_os <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.os.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    glm.data <- glmnet.os(data, input$n_folds, input$alph)[[3]]
    shiny::validate(
      need(!is.null(glm.data),
           "No gene (coefficient) found in model!")
    )
    fit <- coxph(Surv(OS.time, OS) ~ unlist(prog_ind), data = glm.data, model = TRUE, x = TRUE, y = TRUE)
    shiny::validate (
      need(fit$iter <= 20,
           "ROC plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(fit$coefficients)),
           "ROC plot not drawn due to infinite coeff(s)!")
    )
    roc.score <- Score(list("Cox(fit$formula[[3]])" = fit), formula = Surv(OS.time, OS) ~ 1,
                       data = glm.data, plots = "ROC", metrics = "AUC")
    plotROC(roc.score)
  })

  output$glm_roc_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_glm", "_roc", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.os.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      glm.data <- glmnet.os(data, input$n_folds, input$alph)[[3]]
      shiny::validate(
        need(!is.null(glm.data),
             "No gene (coefficient) found in model!")
      )
      fit <- coxph(Surv(OS.time, OS) ~ unlist(prog_ind), data = glm.data, model = TRUE, x = TRUE, y = TRUE)
      shiny::validate (
        need(fit$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
      roc.score <- Score(list("Cox(fit$formula[[3]])" = fit), formula = Surv(OS.time, OS) ~ 1,
                         data = glm.data, plots = "ROC", metrics = "AUC")
      ggsave(file, print(plotROC(roc.score)), dpi = 320)
    }
  )
  
  output$glm_roc_dss <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dss.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    glm.data <- glmnet.dss(data, input$n_folds, input$alph)[[3]]
    shiny::validate(
      need(!is.null(glm.data),
           "No gene (coefficient) found in model!")
    )
    fit <- coxph(Surv(DSS.time, DSS) ~ unlist(prog_ind), data = glm.data, model = TRUE, x = TRUE, y = TRUE)
    shiny::validate (
      need(fit$iter <= 20,
           "ROC plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(fit$coefficients)),
           "ROC plot not drawn due to infinite coeff(s)!")
    )
    roc.score <- Score(list("Cox(fit$formula[[3]])" = fit), formula = Surv(DSS.time, DSS) ~ 1,
                       data = glm.data, plots = "ROC", metrics = "AUC")
    plotROC(roc.score)
  })


  output$glm_roc_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_glm", "_roc", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dss.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      glm.data <- glmnet.dss(data, input$n_folds, input$alph)[[3]]
      shiny::validate(
        need(!is.null(glm.data),
             "No gene (coefficient) found in model!")
      )
      fit <- coxph(Surv(DSS.time, DSS) ~ unlist(prog_ind), data = glm.data, model = TRUE, x = TRUE, y = TRUE)
      shiny::validate (
        need(fit$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
      roc.score <- Score(list("Cox(fit$formula[[3]])" = fit), formula = Surv(DSS.time, DSS) ~ 1,
                         data = glm.data, plots = "ROC", metrics = "AUC")
      ggsave(file, print(plotROC(roc.score)), dpi = 320)
    }
  )

  output$glm_roc_dfi <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    glm.data <- glmnet.dfi(data, input$n_folds, input$alph)[[3]]
    shiny::validate(
      need(!is.null(glm.data),
           "No gene (coefficient) found in model!")
    )
    fit <- coxph(Surv(DFI.time, DFI) ~ unlist(prog_ind), data = glm.data, model = TRUE, x = TRUE, y = TRUE)
    shiny::validate (
      need(fit$iter <= 20,
           "ROC plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(fit$coefficients)),
           "ROC plot not drawn due to infinite coeff(s)!")
    )
    roc.score <- Score(list("Cox(fit$formula[[3]])" = fit), formula = Surv(DFI.time, DFI) ~ 1,
                       data = glm.data, plots = "ROC", metrics = "AUC")
    plotROC(roc.score)
  })

  output$glm_roc_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_glm", "_roc", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      glm.data <- glmnet.dfi(data, input$n_folds, input$alph)[[3]]
      shiny::validate(
        need(!is.null(glm.data),
             "No gene (coefficient) found in model!")
      )
      fit <- coxph(Surv(DFI.time, DFI) ~ unlist(prog_ind), data = glm.data, model = TRUE, x = TRUE, y = TRUE)
      shiny::validate (
        need(fit$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
      roc.score <- Score(list("Cox(fit$formula[[3]])" = fit), formula = Surv(DFI.time, DFI) ~ 1,
                         data = glm.data, plots = "ROC", metrics = "AUC")
      ggsave(file, print(plotROC(roc.score)), dpi = 320)
    }
  )

  output$glm_roc_pfi <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.pfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    glm.data <- glmnet.pfi(data, input$n_folds, input$alph)[[3]]
    shiny::validate(
      need(!is.null(glm.data),
           "No gene (coefficient) found in model!")
    )
    fit <- coxph(Surv(PFI.time, PFI) ~ unlist(prog_ind), data = glm.data, model = TRUE, x = TRUE, y = TRUE)
    shiny::validate (
      need(fit$iter <= 20,
           "ROC plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(fit$coefficients)),
           "ROC plot not drawn due to infinite coeff(s)!")
    )
    roc.score <- Score(list("Cox(fit$formula[[3]])" = fit), formula = Surv(PFI.time, PFI) ~ 1,
                       data = glm.data, plots = "ROC", metrics = "AUC")
    plotROC(roc.score)
  })

  output$glm_roc_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_glm", "_roc", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.pfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      glm.data <- glmnet.pfi(data, input$n_folds, input$alph)[[3]]
      shiny::validate(
        need(!is.null(glm.data),
             "No gene (coefficient) found in model!")
      )
      fit <- coxph(Surv(PFI.time, PFI) ~ unlist(prog_ind), data = glm.data, model = TRUE, x = TRUE, y = TRUE)
      shiny::validate (
        need(fit$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
      roc.score <- Score(list("Cox(fit$formula[[3]])" = fit), formula = Surv(PFI.time, PFI) ~ 1,
                         data = glm.data, plots = "ROC", metrics = "AUC")
      ggsave(file, print(plotROC(roc.score)), dpi = 320)
   }
  )
  
########################################################################  
  output$gic_os <- renderPlot({
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.os.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      plot(bess(exp.data, surv.data, family = "cox"))
    })
  
  output$gic_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_bess", "_gic", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.os.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      ggsave(file, plot(bess(exp.data, surv.data, family = "cox")), dpi = 320)
    }
  )
  
  output$gic_dss <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dss.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    plot(bess(exp.data, surv.data, family = "cox"))  
    })
  
  output$gic_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_bess", "_gic", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dss.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      plot(bess(exp.data, surv.data, family = "cox"))  
      ggsave(file, plot(bess(exp.data, surv.data, family = "cox")), dpi = 320)
    }
  )
  
  output$gic_dfi <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    plot(bess(exp.data, surv.data, family = "cox"))   
  })
  
  output$gic_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_bess", "_gic", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      ggsave(file, plot(bess(exp.data, surv.data, family = "cox")), dpi = 320)
    }
  )
  
  output$gic_pfi <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.pfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    plot(bess(exp.data, surv.data, family = "cox"))    
  })
  
  output$gic_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_bess", "_gic", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.pfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      ggsave(file, plot(bess(exp.data, surv.data, family = "cox")), dpi = 320)
    }
  )
  
  output$pred_bess_os <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.os.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    pr <- bess.os(bess, data)[[1]]
    shiny::validate(
      need(nrow(pr) != 0,
           "No gene (coefficient) found in model!")
    )
    ggplot(pr, aes(x = genes, y = coeffs)) + geom_point(col = "darkred", size = 0.5) + xlab("genes") + ylab("coefficients")
  })

  output$pred_bess_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_bess", "_preds", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.os.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      pr <- bess.os(bess, data)[[1]]
      shiny::validate(
        need(nrow(pr) != 0,
             "No gene (coefficient) found in model!")
      )
      ggsave(file, ggplot(pr, aes(x = genes, y = coeffs)) + geom_point(col = "darkred", size = 0.5) + 
               xlab("genes") + ylab("coefficients"), dpi = 320)
    }
  )

  output$pred_bess_dss <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dss.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    pr <- bess.dss(bess, data)[[1]]
    shiny::validate(
      need(nrow(pr) != 0,
           "No gene (coefficient) found in model!")
    )  
    ggplot(pr, aes(x = genes, y = coeffs)) + geom_point(col = "darkred", size = 0.5) + xlab("genes") + ylab("coefficients")
  })

  output$pred_bess_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_bess", "_preds", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dss.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      pr <- bess.dss(bess, data)[[1]]
      shiny::validate(
        need(nrow(pr) != 0,
             "No gene (coefficient) found in model!")
      )
      ggsave(file, ggplot(pr, aes(x = genes, y = coeffs)) + geom_point(col = "darkred", size = 0.5) + 
               xlab("genes") + ylab("coefficients"), dpi = 320)
    }
  )

  output$pred_bess_dfi <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    pr <- bess.dfi(bess, data)[[1]]
    shiny::validate(
      need(nrow(pr) != 0,
           "No gene (coefficient) found in model!")
    )
    ggplot(pr, aes(x = genes, y = coeffs)) + geom_point(col = "darkred", size = 0.5) + xlab("genes") + ylab("coefficients")
  })

  output$pred_bess_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_bess", "_preds", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      pr <- bess.dfi(bess, data)[[1]]
      shiny::validate(
        need(nrow(pr) != 0,
             "No gene (coefficient) found in model!")
      )
      ggsave(file, ggplot(pr, aes(x = genes, y = coeffs)) + geom_point(col = "darkred", size = 0.5) +
               xlab("genes") + ylab("coefficients"), dpi = 320)
    }
  )

  output$pred_bess_pfi <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.pfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    pr <- bess.pfi(bess, data)[[1]]
    shiny::validate(
      need(nrow(pr) != 0,
           "No gene (coefficient) found in model!")
    )   
    ggplot(pr, aes(x = genes, y = coeffs)) + geom_point(col = "darkred", size = 0.5) + xlab("genes") + ylab("coefficients")
  })

  output$pred_bess_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_bess", "_preds", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.pfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      pr <- bess.pfi(bess, data)[[1]]
      shiny::validate(
        need(nrow(pr) != 0,
             "No gene (coefficient) found in model!")
      ) 
      ggsave(file, ggplot(pr, aes(x = genes, y = coeffs)) + geom_point(col = "darkred", size = 0.5) +
               xlab("genes") + ylab("coefficients"), dpi = 320)
    }
  )
  
  output$bess_coef_tab_os <- renderTable({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.os.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    pr <- bess.os(bess, data)[[1]]
    shiny::validate(
      need(nrow(pr) != 0,
           "No gene (coefficient) found in model!")
    )
    pr
  })
  
  output$bess_coef_tab_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_bess_coeffs", ".tsv", sep = "")},
    content = function(file) {
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.os.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    pr <- bess.os(bess, data)[[1]]
    shiny::validate(
      need(nrow(pr) != 0,
           "No gene (coefficient) found in model!")
    )
    write.table(pr, file)}
  )
  
  output$bess_coef_tab_dss <- renderTable({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dss.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    pr <- bess.dss(bess, data)[[1]]
    shiny::validate(
      need(nrow(pr) != 0,
           "No gene (coefficient) found in model!")
    )
    pr
  })
  
  output$bess_coef_tab_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_bess_coeffs", ".tsv", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dss.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      pr <- bess.dss(bess, data)[[1]]
      shiny::validate(
        need(nrow(pr) != 0,
             "No gene (coefficient) found in model!")
      )
      write.table(pr, file)}
  )
  
  output$bess_coef_tab_dfi <- renderTable({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    pr <- bess.dfi(bess, data)[[1]]
    shiny::validate(
      need(nrow(pr) != 0,
           "No gene (coefficient) found in model!")
    )
    pr
  })
  
  output$bess_coef_tab_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_bess_coeffs", ".tsv", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      pr <- bess.dfi(bess, data)[[1]]
      shiny::validate(
        need(nrow(pr) != 0,
             "No gene (coefficient) found in model!")
      )
      write.table(pr, file)}
  )
  
  output$bess_coef_tab_pfi <- renderTable({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.pfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    pr <- bess.pfi(bess, data)[[1]]
    shiny::validate(
      need(nrow(pr) != 0,
           "No gene (coefficient) found in model!")
    )
    pr
  })
  
  output$bess_coef_tab_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_bess_coeffs", ".tsv", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.pfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      pr <- bess.pfi(bess, data)[[1]]
      shiny::validate(
        need(nrow(pr) != 0,
             "No gene (coefficient) found in model!")
      )
      write.table(pr, file)}
  )

  output$bess_data_os <- renderTable({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.os.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    km.data <- bess.os(bess, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    km.data
  })

  output$bess_data_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_bessPI", "_data", ".tsv", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.os.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      km.data <- bess.os(bess, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      km.data$prog_ind <- km.data$prog_ind[[1]]
      write.table(km.data, file)}
  )

  output$bess_data_dss <- renderTable({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dss.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    km.data <- bess.dss(bess, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    km.data
  })

  output$bess_data_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_bessPI", "_data", ".tsv", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dss.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      km.data <- bess.dss(bess, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      km.data$prog_ind <- km.data$prog_ind[[1]]
      write.table(km.data, file)}
  )

  output$bess_data_dfi <- renderTable({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    km.data <- bess.dfi(bess, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    km.data
  })

  output$bess_data_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_bessPI", "_data", ".tsv", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      km.data <- bess.dfi(bess, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      km.data$prog_ind <- km.data$prog_ind[[1]]
      write.table(km.data, file)}
  )

  output$bess_data_pfi <- renderTable({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.pfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    km.data <- bess.pfi(bess, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    km.data
  })

  output$bess_data_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_bessPI", "_data", ".tsv", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.pfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      km.data <- bess.pfi(bess, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      km.data$prog_ind <- km.data$prog_ind[[1]]
      write.table(km.data, file)}
  )

  output$km_bess_os <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.os.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    km.data <- bess.os(bess, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    surv.fit.os <- survfit(Surv(OS.time, OS) ~ Risk, data = km.data)
    ggsurvplot(surv.fit.os, data = km.data, palette = c("#FF1493", "#0000CD"), pval = TRUE,
               risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")
  })

  output$km_bess_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_bessPI", "_km", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.os.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      km.data <- bess.os(bess, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      surv.fit.os <- survfit(Surv(OS.time, OS) ~ Risk, data = km.data)
      ggsave(file, print(ggsurvplot(surv.fit.os, data = km.data, palette = c("#FF1493", "#0000CD"), pval = TRUE,
                                    risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")), dpi = 320)
    }
  )

  output$km_bess_dss <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dss.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    km.data <- bess.dss(bess, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    surv.fit.dss <- survfit(Surv(DSS.time, DSS) ~ Risk, data = km.data)
    ggsurvplot(surv.fit.dss, data = km.data, palette = c("#FF1493", "#0000CD"), pval = TRUE,
               risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")
  })

  output$km_bess_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_bessPI", "_km", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dss.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      km.data <- bess.dss(bess, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      surv.fit.dss <- survfit(Surv(DSS.time, DSS) ~ Risk, data = km.data)
      ggsave(file, print(ggsurvplot(surv.fit.dss, data = km.data, palette = c("#FF1493", "#0000CD"), pval = TRUE,
                                    risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")), dpi = 320)
    }
  )

  output$km_bess_dfi <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    km.data <- bess.dfi(bess, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    surv.fit.dfi <- survfit(Surv(DFI.time, DFI) ~ Risk, data = km.data)
    ggsurvplot(surv.fit.dfi, data = km.data, palette = c("#FF1493", "#0000CD"), pval = TRUE,
               risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")
  })

  output$km_bess_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_bessPI", "_km", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      km.data <- bess.dfi(bess, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      surv.fit.dfi <- survfit(Surv(DFI.time, DFI) ~ Risk, data = km.data)
      ggsave(file, print(ggsurvplot(surv.fit.dfi, data = km.data, palette = c("#FF1493", "#0000CD"), pval = TRUE,
                                    risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")), dpi = 320)
    }
  )

  output$km_bess_pfi <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.pfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    km.data <- bess.pfi(bess, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    surv.fit.pfi <- survfit(Surv(PFI.time, PFI) ~ Risk, data = km.data)
    ggsurvplot(surv.fit.pfi, data = km.data, palette = c("#FF1493", "#0000CD"), pval = TRUE,
               risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")
  })

  output$km_bess_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_bessPI", "_km", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.pfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      km.data <- bess.pfi(bess, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      surv.fit.pfi <- survfit(Surv(PFI.time, PFI) ~ Risk, data = km.data)
      ggsave(file, print(ggsurvplot(surv.fit.pfi, data = km.data, palette = c("#FF1493", "#0000CD"), pval = TRUE,
                                    risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")), dpi = 320)
    }
  )

  output$bess_cumhz_os <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.os.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    km.data <- bess.os(bess, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    surv.fit.os <- survfit(Surv(OS.time, OS) ~ Risk, data = km.data)
    ggsurvplot(surv.fit.os, data = km.data, palette = c("#FF1493", "#0000CD"), pval = TRUE,
               fun = "cumhaz", risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")
  })

  output$bess_cumhz_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_bessPI", "_cumhz", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.os.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      km.data <- bess.os(bess, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      surv.fit.os <- survfit(Surv(OS.time, OS) ~ Risk, data = km.data)
      ggsave(file, print(ggsurvplot(surv.fit.os, data = km.data, palette = c("#FF1493", "#0000CD"), pval = TRUE,
                                    fun = "cumhaz", risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")), 
                                      dpi = 320)
    }
  )

  output$bess_cumhz_dss <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dss.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    km.data <- bess.dss(bess, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    surv.fit.dss <- survfit(Surv(DSS.time, DSS) ~ Risk, data = km.data)
    ggsurvplot(surv.fit.dss, data = km.data, palette = c("#FF1493", "#0000CD"), pval = TRUE,
               fun = "cumhaz", risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")
  })

  output$bess_cumhz_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_bessPI", "_cumhz", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dss.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      km.data <- bess.dss(bess, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      surv.fit.dss <- survfit(Surv(DSS.time, DSS) ~ Risk, data = km.data)
      ggsave(file, print(ggsurvplot(surv.fit.dss, data = km.data, palette = c("#FF1493", "#0000CD"), pval = TRUE,
                                    fun = "cumhaz", risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")), 
                                      dpi = 320)
    }
  )

  output$bess_cumhz_dfi <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    km.data <- bess.dfi(bess, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    surv.fit.dfi <- survfit(Surv(DFI.time, DFI) ~ Risk, data = km.data)
    ggsurvplot(surv.fit.dfi, data = km.data, palette = c("#FF1493", "#0000CD"), pval = TRUE,
               fun = "cumhaz", risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")
  })

  output$bess_cumhz_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_bessPI", "_cumhz", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      km.data <- bess.dfi(bess, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      surv.fit.dfi <- survfit(Surv(DFI.time, DFI) ~ Risk, data = km.data)
      ggsave(file, print(ggsurvplot(surv.fit.dfi, data = km.data, palette = c("#FF1493", "#0000CD"), pval = TRUE,
                                    fun = "cumhaz", risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")), 
                                      dpi = 320)
    }
  )

  output$bess_cumhz_pfi <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.pfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    km.data <- bess.pfi(bess, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    surv.fit.pfi <- survfit(Surv(PFI.time, PFI) ~ Risk, data = km.data)
    ggsurvplot(surv.fit.pfi, data = km.data, palette = c("#FF1493", "#0000CD"), pval = TRUE,
               fun = "cumhaz", risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")
  })

  output$bess_cumhz_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_bessPI", "_cumhz", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.pfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      km.data <- bess.pfi(bess, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      surv.fit.pfi <- survfit(Surv(PFI.time, PFI) ~ Risk, data = km.data)
      ggsave(file, print(ggsurvplot(surv.fit.pfi, data = km.data, palette = c("#FF1493", "#0000CD"), pval = TRUE,
                                    fun = "cumhaz", risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")), 
                                      dpi = 320)
    }
  )

  output$bess_roc_os <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.os.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    km.data <- bess.os(bess, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    fit <- coxph(Surv(OS.time, OS) ~ unlist(prog_ind), data = km.data, model = TRUE, x = TRUE, y = TRUE)
    shiny::validate (
      need(fit$iter <= 20,
           "ROC plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(fit$coefficients)),
           "ROC plot not drawn due to infinite coeff(s)!")
    )
    roc.score <- Score(list("Cox(fit$formula[[3]])" = fit), formula = Surv(OS.time, OS) ~ 1,
                       data = km.data, plots = "ROC", metrics = "AUC")
    plotROC(roc.score)
  })

  output$bess_roc_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_bess", "_roc", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.os.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("OS", "OS.time")])))) >= nrow(data[, c("OS", "OS.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      km.data <- bess.os(bess, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      fit <- coxph(Surv(OS.time, OS) ~ unlist(prog_ind), data = km.data, model = TRUE, x = TRUE, y = TRUE)
      shiny::validate (
        need(fit$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
      roc.score <- Score(list("Cox(fit$formula[[3]])" = fit), formula = Surv(OS.time, OS) ~ 1,
                         data = km.data, plots = "ROC", metrics = "AUC")
      ggsave(file, print(plotROC(roc.score)), dpi = 320)
    }
  )

  output$bess_roc_dss <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dss.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    km.data <- bess.dss(bess, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    fit <- coxph(Surv(DSS.time, DSS) ~ unlist(prog_ind), data = km.data, model = TRUE, x = TRUE, y = TRUE)
    shiny::validate (
      need(fit$iter <= 20,
           "ROC plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(fit$coefficients)),
           "ROC plot not drawn due to infinite coeff(s)!")
    )
    roc.score <- Score(list("Cox(fit$formula[[3]])" = fit), formula = Surv(DSS.time, DSS) ~ 1,
                       data = km.data, plots = "ROC", metrics = "AUC")
    plotROC(roc.score)
  })

  output$bess_roc_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_bess", "_roc", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dss.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DSS", "DSS.time")])))) >= nrow(data[, c("DSS", "DSS.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      km.data <- bess.dss(bess, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      fit <- coxph(Surv(DSS.time, DSS) ~ unlist(prog_ind), data = km.data, model = TRUE, x = TRUE, y = TRUE)
      shiny::validate (
        need(fit$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
      roc.score <- Score(list("Cox(fit$formula[[3]])" = fit), formula = Surv(DSS.time, DSS) ~ 1,
                         data = km.data, plots = "ROC", metrics = "AUC")
      ggsave(file, print(plotROC(roc.score)), dpi = 320)
    }
  )

  output$bess_roc_dfi <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.dfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    km.data <- bess.dfi(bess, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    fit <- coxph(Surv(DFI.time, DFI) ~ unlist(prog_ind), data = km.data, model = TRUE, x = TRUE, y = TRUE)
    shiny::validate (
      need(fit$iter <= 20,
           "ROC plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(fit$coefficients)),
           "ROC plot not drawn due to infinite coeff(s)!")
    )
    roc.score <- Score(list("Cox(fit$formula[[3]])" = fit), formula = Surv(DFI.time, DFI) ~ 1,
                       data = km.data, plots = "ROC", metrics = "AUC")
    plotROC(roc.score)
  })

  output$bess_roc_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_bess", "_roc", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.dfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("DFI", "DFI.time")])))) >= nrow(data[, c("DFI", "DFI.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      km.data <- bess.dfi(bess, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      fit <- coxph(Surv(DFI.time, DFI) ~ unlist(prog_ind), data = km.data, model = TRUE, x = TRUE, y = TRUE)
      shiny::validate (
        need(fit$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
      roc.score <- Score(list("Cox(fit$formula[[3]])" = fit), formula = Surv(DFI.time, DFI) ~ 1,
                         data = km.data, plots = "ROC", metrics = "AUC")
      ggsave(file, print(plotROC()), dpi = 320)
    }
  )

  output$bess_roc_pfi <- renderPlot({
    req(length(input$genes) >= 2)
    data <- isolate(react.cox.pfi.data())
    shiny::validate(
      need(input$genes %in% colnames(data),
           "Selected genes not found in input dataset!")
    )
    shiny::validate(
      need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
    surv.data <- as.matrix(na.omit(data)[3:2])
    bess <- bess(exp.data, surv.data, family = "cox")
    km.data <- bess.pfi(bess, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    fit <- coxph(Surv(PFI.time, PFI) ~ unlist(prog_ind), data = km.data, model = TRUE, x = TRUE, y = TRUE)
    shiny::validate (
      need(fit$iter <= 20,
           "ROC plot not drawn due to unconverging model!")
    )
    shiny::validate (
      need(all(!is.infinite(fit$coefficients)),
           "ROC plot not drawn due to infinite coeff(s)!")
    )
    roc.score <- Score(list("Cox(fit$formula[[3]])" = fit), formula = Surv(PFI.time, PFI) ~ 1,
                       data = km.data, plots = "ROC", metrics = "AUC")
    plotROC(roc.score)
  })

  output$bess_roc_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_bess", "_roc", ".png", sep = "")},
    content = function(file) {
      req(length(input$genes) >= 2)
      data <- isolate(react.cox.pfi.data())
      shiny::validate(
        need(input$genes %in% colnames(data),
             "Selected genes not found in input dataset!")
      )
      shiny::validate(
        need(!length(which(is.na(rowSums(data[, c("PFI", "PFI.time")])))) >= nrow(data[, c("PFI", "PFI.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      exp.data <- as.matrix(na.omit(data)[4:ncol(data)])
      surv.data <- as.matrix(na.omit(data)[3:2])
      bess <- bess(exp.data, surv.data, family = "cox")
      km.data <- bess.pfi(bess, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      fit <- coxph(Surv(PFI.time, PFI) ~ unlist(prog_ind), data = km.data, model = TRUE, x = TRUE, y = TRUE)
      shiny::validate (
        need(fit$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
      roc.score <- Score(list("Cox(fit$formula[[3]])" = fit), formula = Surv(PFI.time, PFI) ~ 1,
                         data = km.data, plots = "ROC", metrics = "AUC")
      ggsave(file, print(plotROC(roc.score)), dpi = 320)
    }
  )

#########################################################
  output$km_os <- renderPlot({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("OS", "OS.time")])))) >= nrow(react.surv()[, c("OS", "OS.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.os(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.os(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    }
    data <- react.cox.os.data()
    km.data <- km.os(fit, data)[[2]]
    shiny::validate(
      need(!is.null(km.data), 
           "No gene (coefficient) found in model!")
    )
    surv.fit.os <- survfit(Surv(OS.time, OS) ~ Risk, data = km.data)
    ggsurvplot(surv.fit.os, data = km.data, palette = c("#CC0066", "#00CCCC"), pval = TRUE,
               risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")
  })
  

  output$km_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_PI", "_km", ".png", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("OS", "OS.time")])))) >= nrow(react.surv()[, c("OS", "OS.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.os(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      } else {
        fit <- multi.cox.model.os(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      }
      data <- react.cox.os.data()
      km.data <- km.os(fit, data)[[2]]
      shiny::validate(
        need(!is.null(km.data), 
             "No gene (coefficient) found in model!")
      )
      surv.fit.os <- survfit(Surv(OS.time, OS) ~ Risk, data = km.data)
    ggsave(file, print(ggsurvplot(surv.fit.os, data = km.data, palette = c("#CC0066", "#00CCCC"), pval = TRUE,
                                  risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")), dpi = 320)
    }
  )
  
  output$km_dss <- renderPlot({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("DSS", "DSS.time")])))) >= nrow(react.surv()[, c("DSS", "DSS.time")])*(99/100),
           "Insufficient data for disease-specific survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.dss(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.dss(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    }
    data <- react.cox.dss.data()
    km.data <- km.dss(fit, data)[[2]]
    shiny::validate(
      need(!is.null(km.data), 
           "No gene (coefficient) found in model!")
    )
    surv.fit.dss <- survfit(Surv(DSS.time, DSS) ~ Risk, data = km.data)
    ggsurvplot(surv.fit.dss, data = km.data, palette = c("#CC0066", "#00CCCC"), pval = TRUE,
               risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")
  })
  
  output$km_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_PI", "_km", ".png", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("DSS", "DSS.time")])))) >= nrow(react.surv()[, c("DSS", "DSS.time")])*(99/100),
             "Insufficient data for disease-specific survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.dss(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      } else {
        fit <- multi.cox.model.dss(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      }
      data <- react.cox.dss.data()
      km.data <- km.dss(fit, data)[[2]]
      shiny::validate(
        need(!is.null(km.data), 
             "No gene (coefficient) found in model!")
      )
      surv.fit.dss <- survfit(Surv(DSS.time, DSS) ~ Risk, data = km.data)
      ggsave(file, print(ggsurvplot(surv.fit.dss, data = km.data, palette = c("#CC0066", "#00CCCC"), pval = TRUE,
                                    risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")), dpi = 320)
    }
  )
  
  output$km_dfi <- renderPlot({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("DFI", "DFI.time")])))) >= nrow(react.surv()[, c("DFI", "DFI.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.dfi(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.dfi(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    }
    data <- react.cox.dfi.data()
    km.data <- km.dfi(fit, data)[[2]]
    shiny::validate(
      need(!is.null(km.data), 
           "No gene (coefficient) found in model!")
    )
    surv.fit.dfi <- survfit(Surv(DFI.time, DFI) ~ Risk, data = km.data)
    ggsurvplot(surv.fit.dfi, data = km.data, palette = c("#CC0066", "#00CCCC"), pval = TRUE,
               risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")
  })
  
  output$km_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_PI", "_km", ".png", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("DFI", "DFI.time")])))) >= nrow(react.surv()[, c("DFI", "DFI.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.dfi(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      } else {
        fit <- multi.cox.model.dfi(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      }
      data <- react.cox.dfi.data()
      km.data <- km.dfi(fit, data)[[2]]
      shiny::validate(
        need(!is.null(km.data), 
             "No gene (coefficient) found in model!")
      )
      surv.fit.dfi <- survfit(Surv(DFI.time, DFI) ~ Risk, data = km.data)
      ggsave(file, print(ggsurvplot(surv.fit.dfi, data = km.data, palette = c("#CC0066", "#00CCCC"), pval = TRUE,
                                    risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")), dpi = 320)
    }
  )
  
  output$km_pfi <- renderPlot({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("PFI", "PFI.time")])))) >= nrow(react.surv()[, c("PFI", "PFI.time")])*(99/100),
           "Insufficient data for progress-free survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.pfi(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.pfi(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    }
    data <- react.cox.pfi.data()
    km.data <- km.pfi(fit, data)[[2]]
    shiny::validate(
      need(!is.null(km.data), 
           "No gene (coefficient) found in model!")
    )
    surv.fit.pfi <- survfit(Surv(PFI.time, PFI) ~ Risk, data = km.data)
    ggsurvplot(surv.fit.pfi, data = km.data, palette = c("#CC0066", "#00CCCC"), pval = TRUE,
               risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")
  })
  
  output$km_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_PI", "_km", ".png", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("PFI", "PFI.time")])))) >= nrow(react.surv()[, c("PFI", "PFI.time")])*(99/100),
             "Insufficient data for progress-free survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.pfi(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      } else {
        fit <- multi.cox.model.pfi(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      }
      data <- react.cox.pfi.data()
      km.data <- km.pfi(fit, data)[[2]]
      shiny::validate(
        need(!is.null(km.data), 
             "No gene (coefficient) found in model!")
      )
      surv.fit.pfi <- survfit(Surv(PFI.time, PFI) ~ Risk, data = km.data)
      ggsave(file, print(ggsurvplot(surv.fit.pfi, data = km.data, risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE,
                                    palette = c("#CC0066", "#00CCCC"), pval = TRUE) + xlab("days")), dpi = 320)
    }
  )
  
  output$pred_os <- renderPlot({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("OS", "OS.time")])))) >= nrow(react.surv()[, c("OS", "OS.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.os(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.os(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    }
    data <- react.cox.os.data()
    coeffs <- km.os(fit, data)[[1]]
    shiny::validate(
      need(nrow(coeffs) != 0,
           "No gene (coefficient) found in model!")
    )
    coeffs <- cbind(rownames(coeffs), as.data.frame(coeffs))
    colnames(coeffs)[1] <- c("genes")
    ggplot(coeffs, aes(x = genes, y = coef)) + geom_point(col = "darkred", size = 0.5) + xlab("genes") + ylab("CPH coefficients")
  })


  output$pred_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_PI", "_coeff", ".png", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("OS", "OS.time")])))) >= nrow(react.surv()[, c("OS", "OS.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.os(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )

      } else {
        fit <- multi.cox.model.os(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      }
      data <- react.cox.os.data()
      coeffs <- km.os(fit, data)[[1]]
      shiny::validate(
        need(nrow(coeffs) != 0,
             "No gene (coefficient) found in model!")
      )
      coeffs <- cbind(rownames(coeffs), as.data.frame(coeffs))
      colnames(coeffs)[1] <- c("genes")
      ggsave(file, ggplot(coeffs, aes(x = genes, y = coef)) + geom_point(col = "darkred", size = 0.5) + 
               xlab("genes") + ylab("CPH coefficients"), dpi = 320)
    }
  )

  output$pred_dss <- renderPlot({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("DSS", "DSS.time")])))) >= nrow(react.surv()[, c("DSS", "DSS.time")])*(99/100),
           "Insufficient data for disease-specific survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.dss(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.dss(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    }
    data <- react.cox.dss.data()
    coeffs <- km.dss(fit, data)[[1]]
    shiny::validate(
      need(nrow(coeffs) != 0,
           "No gene (coefficient) found in model!")
    )
    coeffs <- cbind(rownames(coeffs), as.data.frame(coeffs))
    colnames(coeffs)[1] <- c("genes")
    ggplot(coeffs, aes(x = genes, y = coef)) + geom_point(col = "darkred", size = 0.5) + xlab("genes") + ylab("CPH coefficients")
  })
  
  
  output$pred_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_PI", "_coeff", ".png", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("DSS", "DSS.time")])))) >= nrow(react.surv()[, c("DSS", "DSS.time")])*(99/100),
             "Insufficient data for disease-specific survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.dss(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
        
      } else {
        fit <- multi.cox.model.dss(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      }
      data <- react.cox.dss.data()
      coeffs <- km.dss(fit, data)[[1]]
      shiny::validate(
        need(nrow(coeffs) != 0,
             "No gene (coefficient) found in model!")
      )
      coeffs <- cbind(rownames(coeffs), as.data.frame(coeffs))
      colnames(coeffs)[1] <- c("genes")
      ggsave(file, ggplot(coeffs, aes(x = genes, y = coef)) + geom_point(col = "darkred", size = 0.5) + 
               xlab("genes") + ylab("CPH coefficients"), dpi = 320)
    }
  )
  
  output$pred_dfi <- renderPlot({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("DFI", "DFI.time")])))) >= nrow(react.surv()[, c("DFI", "DFI.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.dfi(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.dfi(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    }
    data <- react.cox.dfi.data()
    coeffs <- km.dfi(fit, data)[[1]]
    shiny::validate(
      need(nrow(coeffs) != 0,
           "No gene (coefficient) found in model!")
    )
    coeffs <- cbind(rownames(coeffs), as.data.frame(coeffs))
    colnames(coeffs)[1] <- c("genes")
    ggplot(coeffs, aes(x = genes, y = coef)) + geom_point(col = "darkred", size = 0.5) + xlab("genes") + ylab("CPH coefficients")
  })
  
  
  output$pred_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_PI", "_coeff", ".png", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("DFI", "DFI.time")])))) >= nrow(react.surv()[, c("DFI", "DFI.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.dfi(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
        
      } else {
        fit <- multi.cox.model.dfi(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      }
      data <- react.cox.dfi.data()
      coeffs <- km.dfi(fit, data)[[1]]
      shiny::validate(
        need(nrow(coeffs) != 0,
             "No gene (coefficient) found in model!")
      )
      coeffs <- cbind(rownames(coeffs), as.data.frame(coeffs))
      colnames(coeffs)[1] <- c("genes")
      ggsave(file, ggplot(coeffs, aes(x = genes, y = coef)) + geom_point(col = "darkred", size = 0.5) +
               xlab("genes") + ylab("CPH coefficients"), dpi = 320)
    }
  )
  
  output$pred_pfi <- renderPlot({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("PFI", "PFI.time")])))) >= nrow(react.surv()[, c("PFI", "PFI.time")])*(99/100),
           "Insufficient data for progress-free survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.pfi(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.pfi(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    }
    data <- react.cox.pfi.data()
    coeffs <- km.pfi(fit, data)[[1]]
    shiny::validate(
      need(nrow(coeffs) != 0,
           "No gene (coefficient) found in model!")
    )
    coeffs <- cbind(rownames(coeffs), as.data.frame(coeffs))
    colnames(coeffs)[1] <- c("genes")
    ggplot(coeffs, aes(x = genes, y = coef)) + geom_point(col = "darkred", size = 0.5) + xlab("genes") + ylab("CPH coefficients")
  })
  
  
  output$pred_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_PI", "_coeff", ".png", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("PFI", "PFI.time")])))) >= nrow(react.surv()[, c("PFI", "PFI.time")])*(99/100),
             "Insufficient data for progress-free survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.pfi(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
        
      } else {
        fit <- multi.cox.model.pfi(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      }
      data <- react.cox.pfi.data()
      coeffs <- km.pfi(fit, data)[[1]]
      shiny::validate(
        need(nrow(coeffs) != 0,
             "No gene (coefficient) found in model!")
      )
      coeffs <- cbind(rownames(coeffs), as.data.frame(coeffs))
      colnames(coeffs)[1] <- c("genes")
      ggsave(file, ggplot(coeffs, aes(x = genes, y = coef)) + geom_point(col = "darkred", size = 0.5) + 
               xlab("genes") + ylab("CPH coefficients"), dpi = 320)
    }
  )

  output$coef_tab_os <- renderTable({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("OS", "OS.time")])))) >= nrow(react.surv()[, c("OS", "OS.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.os(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.os(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    }
    data <- react.cox.os.data()
    coeffs <- km.os(fit, data)[[1]]
    shiny::validate(
      need(nrow(coeffs) != 0,
           "No gene (coefficient) found in model!")
    )
    coeffs <- cbind(rownames(coeffs), as.data.frame(coeffs))
    colnames(coeffs)[1] <- c("genes")
    coeffs
  })
  
  
  output$coef_tab_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_coeff_tab", ".tsv", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("OS", "OS.time")])))) >= nrow(react.surv()[, c("OS", "OS.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.os(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
        
      } else {
        fit <- multi.cox.model.os(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      }
      data <- react.cox.os.data()
      coeffs <- km.os(fit, data)[[1]]
      shiny::validate(
        need(nrow(coeffs) != 0,
             "No gene (coefficient) found in model!")
      )
      coeffs <- cbind(rownames(coeffs), as.data.frame(coeffs))
      colnames(coeffs)[1] <- c("genes")
      write.table(coeffs, file)
    }
  )
  
  output$coef_tab_dss <- renderTable({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("DSS", "DSS.time")])))) >= nrow(react.surv()[, c("DSS", "DSS.time")])*(99/100),
           "Insufficient data for disease-specific survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.dss(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.dss(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    }
    data <- react.cox.dss.data()
    coeffs <- km.dss(fit, data)[[1]]
    shiny::validate(
      need(nrow(coeffs) != 0,
           "No gene (coefficient) found in model!")
    )
    coeffs <- cbind(rownames(coeffs), as.data.frame(coeffs))
    colnames(coeffs)[1] <- c("genes")
    coeffs
  })
  
  
  output$coef_tab_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_coeff_tab", ".tsv", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("DSS", "DSS.time")])))) >= nrow(react.surv()[, c("DSS", "DSS.time")])*(99/100),
             "Insufficient data for disease-specific survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.dss(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
        
      } else {
        fit <- multi.cox.model.dss(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      }
      data <- react.cox.dss.data()
      coeffs <- km.dss(fit, data)[[1]]
      shiny::validate(
        need(nrow(coeffs) != 0,
             "No gene (coefficient) found in model!")
      )
      coeffs <- cbind(rownames(coeffs), as.data.frame(coeffs))
      colnames(coeffs)[1] <- c("genes")
      write.table(coeffs, file)
    }
  )
  
  output$coef_tab_dfi <- renderTable({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("DFI", "DFI.time")])))) >= nrow(react.surv()[, c("DFI", "DFI.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.dfi(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.dfi(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    }
    data <- react.cox.dfi.data()
    coeffs <- km.dfi(fit, data)[[1]]
    shiny::validate(
      need(nrow(coeffs) != 0,
           "No gene (coefficient) found in model!")
    )
    coeffs <- cbind(rownames(coeffs), as.data.frame(coeffs))
    colnames(coeffs)[1] <- c("genes")
    coeffs
  })
  
  
  output$coef_tab_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_coeff_tab", ".tsv", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("DFI", "DFI.time")])))) >= nrow(react.surv()[, c("DFI", "DFI.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.dfi(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
        
      } else {
        fit <- multi.cox.model.dfi(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      }
      data <- react.cox.dfi.data()
      coeffs <- km.dfi(fit, data)[[1]]
      shiny::validate(
        need(nrow(coeffs) != 0,
             "No gene (coefficient) found in model!")
      )
      coeffs <- cbind(rownames(coeffs), as.data.frame(coeffs))
      colnames(coeffs)[1] <- c("genes")
      write.table(coeffs, file)
    }
  )
  
  output$coef_tab_pfi <- renderPlot({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("PFI", "PFI.time")])))) >= nrow(react.surv()[, c("PFI", "PFI.time")])*(99/100),
           "Insufficient data for progress-free survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.pfi(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.pfi(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    }
    data <- react.cox.pfi.data()
    coeffs <- km.pfi(fit, data)[[1]]
    shiny::validate(
      need(nrow(coeffs) != 0,
           "No gene (coefficient) found in model!")
    )
    coeffs <- cbind(rownames(coeffs), as.data.frame(coeffs))
    colnames(coeffs)[1] <- c("genes")
    coeffs
  })
  
  
  output$coef_tab_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_coeff_tab", ".tsv", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("PFI", "PFI.time")])))) >= nrow(react.surv()[, c("PFI", "PFI.time")])*(99/100),
             "Insufficient data for progress-free survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.pfi(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
        
      } else {
        fit <- multi.cox.model.pfi(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      }
      data <- react.cox.pfi.data()
      coeffs <- km.pfi(fit, data)[[1]]
      shiny::validate(
        need(nrow(coeffs) != 0,
             "No gene (coefficient) found in model!")
      )
      coeffs <- cbind(rownames(coeffs), as.data.frame(coeffs))
      colnames(coeffs)[1] <- c("genes")
      write.table(coeffs, file)
    }
  )
  
  
  
  
  
  
  
  
  
  
  
  
  
  output$cumhz_os <- renderPlot({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("OS", "OS.time")])))) >= nrow(react.surv()[, c("OS", "OS.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.os(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.os(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    }
    data <- react.cox.os.data()
    km.data <- km.os(fit, data)[[2]]
    shiny::validate(
    need(!is.null(km.data),
         "No gene (coefficient) found in model!")
    )
    surv.fit.os <- survfit(Surv(OS.time, OS) ~ Risk, data = km.data)
    ggsurvplot(surv.fit.os, data = km.data, palette = c("#990099", "#00FFFF"), pval = TRUE,
             fun = "cumhaz", risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")
  })

  output$cumhz_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_PI", "_cumhz", ".png", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("OS", "OS.time")])))) >= nrow(react.surv()[, c("OS", "OS.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.os(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      } else {
        fit <- multi.cox.model.os(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      }
      data <- react.cox.os.data()
      km.data <- km.os(fit, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      surv.fit.os <- survfit(Surv(OS.time, OS) ~ Risk, data = km.data)
      ggsave(file, print(ggsurvplot(surv.fit.os, data = km.data, palette = c("#990099", "#00FFFF"), pval = TRUE,
                                    fun = "cumhaz", risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")),
                                      dpi = 320)
    }
  )
  
  output$cumhz_dss <- renderPlot({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("DSS", "DSS.time")])))) >= nrow(react.surv()[, c("DSS", "DSS.time")])*(99/100),
           "Insufficient data for disease-specific survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.dss(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.dss(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    }
    data <- react.cox.dss.data()
    km.data <- km.dss(fit, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    surv.fit.dss <- survfit(Surv(DSS.time, DSS) ~ Risk, data = km.data)
    ggsurvplot(surv.fit.dss, data = km.data, palette = c("#990099", "#00FFFF"), pval = TRUE,
               fun = "cumhaz", risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")
  })
  
  output$cumhz_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_PI", "_cumhz", ".png", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("DSS", "DSS.time")])))) >= nrow(react.surv()[, c("DSS", "DSS.time")])*(99/100),
             "Insufficient data for disease-specific survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.dss(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      } else {
        fit <- multi.cox.model.dss(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      }
      data <- react.cox.dss.data()
      km.data <- km.dss(fit, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      surv.fit.dss <- survfit(Surv(DSS.time, DSS) ~ Risk, data = km.data)
      ggsave(file, print(ggsurvplot(surv.fit.dss, data = km.data, palette = c("#990099", "#00FFFF"), pval = TRUE,
                                    fun = "cumhaz", risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")),
                                      dpi = 320)
    }
  )
  
  output$cumhz_dfi <- renderPlot({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("DFI", "DFI.time")])))) >= nrow(react.surv()[, c("DFI", "DFI.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.dfi(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.dfi(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    }
    data <- react.cox.dfi.data()
    km.data <- km.dfi(fit, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    surv.fit.dfi <- survfit(Surv(DFI.time, DFI) ~ Risk, data = km.data)
    ggsurvplot(surv.fit.dfi, data = km.data, palette = c("#990099", "#00FFFF"), pval = TRUE,
               fun = "cumhaz", risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")
  })
  
  output$cumhz_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_PI", "_cumhz", ".png", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("DFI", "DFI.time")])))) >= nrow(react.surv()[, c("DFI", "DFI.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.dfi(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      } else {
        fit <- multi.cox.model.dfi(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      }
      data <- react.cox.dfi.data()
      km.data <- km.dfi(fit, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      surv.fit.dfi <- survfit(Surv(DFI.time, DFI) ~ Risk, data = km.data)
      ggsave(file, print(ggsurvplot(surv.fit.dfi, data = km.data, palette = c("#990099", "#00FFFF"), pval = TRUE,
                                    fun = "cumhaz", risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")),
                                      dpi = 320)
    }
  )
  
  output$cumhz_pfi <- renderPlot({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("PFI", "PFI.time")])))) >= nrow(react.surv()[, c("PFI", "PFI.time")])*(99/100),
           "Insufficient data for progress-free survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.pfi(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.pfi(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "K-M plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "K-M plot not drawn due to infinite coeff(s)!")
      )
    }
    data <- react.cox.pfi.data()
    km.data <- km.pfi(fit, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    surv.fit.pfi <- survfit(Surv(PFI.time, PFI) ~ Risk, data = km.data)
    ggsurvplot(surv.fit.pfi, data = km.data, palette = c("#990099", "#00FFFF"), pval = TRUE,
               fun = "cumhaz", risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")
  })
  
  output$cumhz_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_PI", "_cumhz", ".png", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("PFI", "PFI.time")])))) >= nrow(react.surv()[, c("PFI", "PFI.time")])*(99/100),
             "Insufficient data for progress-free survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.pfi(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      } else {
        fit <- multi.cox.model.pfi(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      }
      data <- react.cox.pfi.data()
      km.data <- km.pfi(fit, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      surv.fit.pfi <- survfit(Surv(PFI.time, PFI) ~ Risk, data = km.data)
      ggsave(file, print(ggsurvplot(surv.fit.pfi, data = km.data, palette = c("#990099", "#00FFFF"), pval = TRUE,
                                    fun = "cumhaz", risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE) + xlab("days")),
                                      dpi = 320)
    }
  )
  
  output$km_tab_os <- renderTable({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("OS", "OS.time")])))) >= nrow(react.surv()[, c("OS", "OS.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.os(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "Unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "Infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.os(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "Unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "Infinite coeff(s)!")
      )
    }
    data <- react.cox.os.data()
    km.data <- km.os(fit, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    km.data
  })
  
  output$km_tab_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_PI", "_data", ".tsv", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("OS", "OS.time")])))) >= nrow(react.surv()[, c("OS", "OS.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.os(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "Unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "Infinite coeff(s)!")
        )
      } else {
        fit <- multi.cox.model.os(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "Unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "Infinite coeff(s)!")
        )
      }
      data <- react.cox.os.data()
      km.data <- km.os(fit, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      km.data$prog_ind <- km.data$prog_ind[[1]]
      write.table(km.data, file)}
  )
  
  output$km_tab_dss <- renderTable({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("DSS", "DSS.time")])))) >= nrow(react.surv()[, c("DSS", "DSS.time")])*(99/100),
           "Insufficient data for disease-specific survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.dss(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "Unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "Infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.dss(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "Unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "Infinite coeff(s)!")
      )
    }
    data <- react.cox.dss.data()
    km.data <- km.dss(fit, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    km.data
  })
  
  output$km_tab_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_PI", "_data", ".tsv", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("DSS", "DSS.time")])))) >= nrow(react.surv()[, c("DSS", "DSS.time")])*(99/100),
             "Insufficient data for disease-specific survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.dss(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "Unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "Infinite coeff(s)!")
        )
      } else {
        fit <- multi.cox.model.dss(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "Unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "Infinite coeff(s)!")
        )
      }
      data <- react.cox.dss.data()
      km.data <- km.dss(fit, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      km.data$prog_ind <- km.data$prog_ind[[1]]
      write.table(km.data, file)}
  )
  
  output$km_tab_dfi <- renderTable({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("DFI", "DFI.time")])))) >= nrow(react.surv()[, c("DFI", "DFI.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.dfi(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "Unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "Infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.dfi(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "Unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "Infinite coeff(s)!")
      )
    }
    data <- react.cox.dfi.data()
    km.data <- km.dfi(fit, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    km.data
  })
  
  output$km_tab_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_PI", "_data", ".tsv", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("DFI", "DFI.time")])))) >= nrow(react.surv()[, c("DFI", "DFI.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.dfi(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "Unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "Infinite coeff(s)!")
        )
      } else {
        fit <- multi.cox.model.dfi(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "Unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "Infinite coeff(s)!")
        )
      }
      data <- react.cox.dfi.data()
      km.data <- km.dfi(fit, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      km.data$prog_ind <- km.data$prog_ind[[1]]
      write.table(km.data, file)}
  )
  
  output$km_tab_pfi <- renderTable({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("PFI", "PFI.time")])))) >= nrow(react.surv()[, c("PFI", "PFI.time")])*(99/100),
           "Insufficient data for progress-free survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.pfi(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "Unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "Infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.pfi(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "Unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "Infinite coeff(s)!")
      )
    }
    data <- react.cox.pfi.data()
    km.data <- km.pfi(fit, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    km.data
  })
  
  output$km_tab_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_PI", "_data", ".tsv", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("PFI", "PFI.time")])))) >= nrow(react.surv()[, c("PFI", "PFI.time")])*(99/100),
             "Insufficient data for progress-free survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.pfi(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "Unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "Infinite coeff(s)!")
        )
      } else {
        fit <- multi.cox.model.pfi(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "Unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "Infinite coeff(s)!")
        )
      }
      data <- react.cox.pfi.data()
      km.data <- km.pfi(fit, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      km.data$prog_ind <- km.data$prog_ind[[1]]
      write.table(km.data, file)}
  )

  output$risk_roc_os <- renderPlot({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("OS", "OS.time")])))) >= nrow(react.surv()[, c("OS", "OS.time")])*(99/100),
           "Insufficient data for overall survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.os(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.os(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
    }
    data <- react.cox.os.data()
    km.data <- km.os(fit, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    km.fit <- coxph(Surv(OS.time, OS) ~ unlist(prog_ind), data = km.data, model = TRUE, x = TRUE, y = TRUE)
    roc.score <- Score(list("Cox(fit$formula[[3]])" = km.fit), formula = Surv(OS.time, OS) ~ 1,
                       data = km.data, plots = "ROC", metrics = "AUC")
    plotROC(roc.score)
  })


  output$risk_roc_os_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_OS", "_PI", "_roc", ".png", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("OS", "OS.time")])))) >= nrow(react.surv()[, c("OS", "OS.time")])*(99/100),
             "Insufficient data for overall survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.os(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      } else {
        fit <- multi.cox.model.os(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      }
      data <- react.cox.os.data()
      km.data <- km.os(fit, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      km.fit <- coxph(Surv(OS.time, OS) ~ unlist(prog_ind), data = km.data, model = TRUE, x = TRUE, y = TRUE)
      roc.score <- Score(list("Cox(fit$formula[[3]])" = km.fit), formula = Surv(OS.time, OS) ~ 1,
                         data = km.data, plots = "ROC", metrics = "AUC")
      ggsave(file, print(plotROC(roc.score)), dpi = 320)
    }
  )
  
  output$risk_roc_dss <- renderPlot({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("DSS", "DSS.time")])))) >= nrow(react.surv()[, c("DSS", "DSS.time")])*(99/100),
           "Insufficient data for disease-specific survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.dss(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.dss(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
    }
    data <- react.cox.dss.data()
    km.data <- km.dss(fit, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    km.fit <- coxph(Surv(DSS.time, DSS) ~ unlist(prog_ind), data = km.data, model = TRUE, x = TRUE, y = TRUE)
    roc.score <- Score(list("Cox(fit$formula[[3]])" = km.fit), formula = Surv(DSS.time, DSS) ~ 1,
                       data = km.data, plots = "ROC", metrics = "AUC")
    plotROC(roc.score)
  })
  
  
  output$risk_roc_dss_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DSS", "_PI", "_roc", ".png", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("DSS", "DSS.time")])))) >= nrow(react.surv()[, c("DSS", "DSS.time")])*(99/100),
             "Insufficient data for disease-specific survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.dss(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      } else {
        fit <- multi.cox.model.dss(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      }
      data <- react.cox.dss.data()
      km.data <- km.dss(fit, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      km.fit <- coxph(Surv(DSS.time, DSS) ~ unlist(prog_ind), data = km.data, model = TRUE, x = TRUE, y = TRUE)
      roc.score <- Score(list("Cox(fit$formula[[3]])" = km.fit), formula = Surv(DSS.time, DSS) ~ 1,
                         data = km.data, plots = "ROC", metrics = "AUC")
      ggsave(file, print(plotROC(roc.score)), dpi = 320)
    }
  )
  
  output$risk_roc_dfi <- renderPlot({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("DFI", "DFI.time")])))) >= nrow(react.surv()[, c("DFI", "DFI.time")])*(99/100),
           "Insufficient data for disease-free survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.dfi(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.dfi(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
    }
    data <- react.cox.dfi.data()
    km.data <- km.dfi(fit, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    km.fit <- coxph(Surv(DFI.time, DFI) ~ unlist(prog_ind), data = km.data, model = TRUE, x = TRUE, y = TRUE)
    roc.score <- Score(list("Cox(fit$formula[[3]])" = km.fit), formula = Surv(DFI.time, DFI) ~ 1,
                       data = km.data, plots = "ROC", metrics = "AUC")
    plotROC(roc.score)
  })
  
  
  output$risk_roc_dfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_DFI", "_PI", "_roc", ".png", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("DFI", "DFI.time")])))) >= nrow(react.surv()[, c("DFI", "DFI.time")])*(99/100),
             "Insufficient data for disease-free survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.dfi(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      } else {
        fit <- multi.cox.model.dfi(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      }
      data <- react.cox.dfi.data()
      km.data <- km.dfi(fit, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      km.fit <- coxph(Surv(DFI.time, DFI) ~ unlist(prog_ind), data = km.data, model = TRUE, x = TRUE, y = TRUE)
      roc.score <- Score(list("Cox(fit$formula[[3]])" = km.fit), formula = Surv(DFI.time, DFI) ~ 1,
                         data = km.data, plots = "ROC", metrics = "AUC")
      ggsave(file, print(plotROC(roc.score)), dpi = 320)
    }
  )
  
  output$risk_roc_pfi <- renderPlot({
    req(input$genes)
    shiny::validate(
      need(!length(which(is.na(rowSums(react.surv()[, c("PFI", "PFI.time")])))) >= nrow(react.surv()[, c("PFI", "PFI.time")])*(99/100),
           "Insufficient data for progress-free survival")
    )
    if(!is.null(input$fold_genes)) {
      fit <- multi.cox.model.pfi(react.surv(), react.exp.fold())
      shiny::validate (
        need(fit$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
    } else {
      fit <- multi.cox.model.pfi(react.surv(), react.exp())
      shiny::validate (
        need(fit$iter <= 20,
             "ROC plot not drawn due to unconverging model!")
      )
      shiny::validate (
        need(all(!is.infinite(fit$coefficients)),
             "ROC plot not drawn due to infinite coeff(s)!")
      )
    }
    data <- react.cox.pfi.data()
    km.data <- km.pfi(fit, data)[[2]]
    shiny::validate(
      need(!is.null(km.data),
           "No gene (coefficient) found in model!")
    )
    km.fit <- coxph(Surv(PFI.time, PFI) ~ unlist(prog_ind), data = km.data, model = TRUE, x = TRUE, y = TRUE)
    roc.score <- Score(list("Cox(fit$formula[[3]])" = km.fit), formula = Surv(PFI.time, PFI) ~ 1,
                       data = km.data, plots = "ROC", metrics = "AUC")
    plotROC(roc.score)
  })
  
  
  output$risk_roc_pfi_down <- downloadHandler(
    filename = function() {paste(react.dataset()$`cancer type abbreviation`[1], "_PFI", "_PI", "_roc", ".png", sep = "")},
    content = function(file) {
      req(input$genes)
      shiny::validate(
        need(!length(which(is.na(rowSums(react.surv()[, c("PFI", "PFI.time")])))) >= nrow(react.surv()[, c("PFI", "PFI.time")])*(99/100),
             "Insufficient data for progress-free survival")
      )
      if(!is.null(input$fold_genes)) {
        fit <- multi.cox.model.pfi(react.surv(), react.exp.fold())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      } else {
        fit <- multi.cox.model.pfi(react.surv(), react.exp())
        shiny::validate (
          need(fit$iter <= 20,
               "K-M plot not drawn due to unconverging model!")
        )
        shiny::validate (
          need(all(!is.infinite(fit$coefficients)),
               "K-M plot not drawn due to infinite coeff(s)!")
        )
      }
      data <- react.cox.pfi.data()
      km.data <- km.pfi(fit, data)[[2]]
      shiny::validate(
        need(!is.null(km.data),
             "No gene (coefficient) found in model!")
      )
      km.fit <- coxph(Surv(PFI.time, PFI) ~ unlist(prog_ind), data = km.data, model = TRUE, x = TRUE, y = TRUE)
      roc.score <- Score(list("Cox(fit$formula[[3]])" = km.fit), formula = Surv(PFI.time, PFI) ~ 1,
                         data = km.data, plots = "ROC", metrics = "AUC")
      ggsave(file, print(plotROC(roc.score)), dpi = 320)
    }
  )
}

