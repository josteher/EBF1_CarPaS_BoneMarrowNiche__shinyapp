### shiny Server

library(shiny)
library(DT)
library(RColorBrewer)
library(FateID)
library(data.table)
library(RaceID)
library(ggpubr)
library(plotly)


source("./RaceID_functions_modified.R")

# Load data for
#           Niche WT only, 
sc_carpas_wt    <- readRDS("./shiny__sc_object__CarPas_female_WT_only_mintotal_2000_minnumber_1_cln_17_probthr_1e-04_FGenes_6_Malat1_CGenes_2_Jun.Rds")
cdiff_carpas_wt <- readRDS("./shiny__cdiff__CarPas_female_WT_only_mintotal_2000_minnumber_1_cln_17_probthr_1e-04_FGenes_6_Malat1_CGenes_2_Jun.Rds")
sc_m_carpas_wt  <- readRDS("./shiny__sc_mice__CarPas_female_WT_only_mintotal_2000_minnumber_1_cln_17_probthr_1e-04_FGenes_6_Malat1_CGenes_2_Jun.Rds")

#           Niche WT + KO, 
sc_carpas_wtko    <- readRDS("./shiny__sc_object_rot__CarPas_female_new_RaceID_mintotal_2000_minnumber_1_cln_21_probthr_1e-04_FGenes_7_Malat1_CGenes_2_Jun.Rds")
cdiff_carpas_wtko <- readRDS("./shiny__cdiff__CarPas_female_new_RaceID_mintotal_2000_minnumber_1_cln_21_probthr_1e-04_FGenes_7_Malat1_CGenes_2_Jun.Rds")
sc_m_carpas_wtko  <- readRDS("./shiny__sc_mice__CarPas_female_new_RaceID_mintotal_2000_minnumber_1_cln_21_probthr_1e-04_FGenes_7_Malat1_CGenes_2_Jun.Rds")



#           HSC,  
sc_hsc_wtko    <- readRDS("./shiny__sc_object__HSC_only_mintotal_4000_cln_5_metric_pearson_probthr_1e-04_FGenes_0__CGenes_9_Kcnq1ot1.Rds")
cdiff_hsc_wtko <- readRDS("./shiny__cdiff__HSC_only_mintotal_4000_cln_5_metric_pearson_probthr_1e-04_FGenes_0__CGenes_9_Kcnq1ot1.Rds")
sc_m_hsc_wtko  <- readRDS("./shiny__sc_mice__HSC_only_mintotal_4000_cln_5_metric_pearson_probthr_1e-04_FGenes_0__CGenes_9_Kcnq1ot1.Rds")



#           HSC + LSK
sc_hsclsk_wtko    <- readRDS("./shiny__sc_object__HSC_LSK_Niche_Ebf1ko_mintotal_3000_cln_13_probthr_1e-05_FGenes_1_Kcnq1ot1_CGenes_2_Mki67.Rds")
cdiff_hsclsk_wtko <- readRDS("./shiny__cdiff__HSC_LSK_Niche_Ebf1ko_mintotal_3000_cln_13_probthr_1e-05_FGenes_1_Kcnq1ot1_CGenes_2_Mki67.Rds")
sc_m_hsclsk_wtko  <- readRDS("./shiny__sc_mice__HSC_LSK_Niche_Ebf1ko_mintotal_3000_cln_13_probthr_1e-05_FGenes_1_Kcnq1ot1_CGenes_2_Mki67.Rds")



# Create lists for all objects
# Create lists for all objects
sc_objects <- list("CarPas_WtKo"=sc_carpas_wtko,    "CarPas_Wt"=sc_carpas_wt,    "CarPas_Ko"=NULL, "HSC"=sc_hsc_wtko ,    "HSC_LSK"=sc_hsclsk_wtko)
cdiffs     <- list("CarPas_WtKo"=cdiff_carpas_wtko, "CarPas_Wt"=cdiff_carpas_wt, "CarPas_Ko"=NULL, "HSC"=cdiff_hsc_wtko , "HSC_LSK"=cdiff_hsclsk_wtko)
sc_mice_ob <- list("CarPas_WtKo"=sc_m_carpas_wtko,  "CarPas_Wt"=sc_m_carpas_wt,  "CarPas_Ko"=NULL, "HSC"=sc_m_hsc_wtko ,  "HSC_LSK"=sc_m_hsclsk_wtko)




# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  
  ###
  ### CAR PaS from WT only
  ###
  
  output$plotmap_cp1 <- renderPlot({
    plotmap_m(sc_objects[["CarPas_Wt"]], cex=1.5)
  })
  
  output$plotsymbolsmap_cp1 <- renderPlot({
    plotsymbolsmap_m(sc_mice_ob[["CarPas_Wt"]], types = substring(colnames(sc_mice_ob[["CarPas_Wt"]]@ndata), first = 1, last=4), 
                     cex=1.5, 
                     leg.pos = "topleft", 
                     leg.cex = 1, 
                     samples_col = c(c("blue",colorRampPalette(c("dodgerblue", "palegreen"))(3)) ))
  })
  

  ########################################################
  ### cdiff output
  ###
  ### Plot a particular cluster
  output$plotclust_cp1 <- renderPlot({
    plotmap_m(sc_objects[["CarPas_Wt"]], my_part = input$cdiff_cl_cp1, cex = 1.5)
  })

  # Make cdiff list and user input reactive
  r_cdiff_cp1  <- reactive(cdiffs[["CarPas_Wt"]])
  cd_p_val_cp1 <- reactive(input$cdiff_pval_cp1)
  cd_p_adj_cp1 <- reactive(input$cdiff_padj_cp1)


  # Show two data tables of cdiff output for one cluster
  output$upgenes_cp1   <- renderDataTable( apply( r_cdiff_cp1()[[input$cdiff_cl_cp1]][r_cdiff_cp1()[[input$cdiff_cl_cp1]]$pv <= cd_p_val_cp1() & r_cdiff_cp1()[[input$cdiff_cl_cp1]]$padj <= cd_p_adj_cp1(),],
                                                  2,
                                                  function(x){round(x,digits = 8)} ) )
  output$downgenes_cp1 <- renderDataTable( apply( r_cdiff_cp1()[[input$cdiff_cl_cp1]][r_cdiff_cp1()[[input$cdiff_cl_cp1]]$pv <= cd_p_val_cp1() & r_cdiff_cp1()[[input$cdiff_cl_cp1]]$padj <= cd_p_adj_cp1(),],
                                                  2,
                                                  function(x){round(x,digits = 8)} ) )


  ########################################################
  # plotexpmap
  ###
  ### Plotexpmap Gene expression in tSNE Map
  output$plotexptsne_cp1       <- renderPlot(plotexpmap_m(sc_objects[["CarPas_Wt"]],input$ptexp_gene_cp1, logsc=F))
  output$plotexptsne_logsc_cp1 <- renderPlot(plotexpmap_m(sc_objects[["CarPas_Wt"]],input$ptexp_gene_cp1, logsc=T, n = paste(input$ptexp_gene_cp1,"log2")))




  ########################################################
  ### Compare two clusters
  output$plotclust_cl1_cp1 <- renderPlot(plotmap_m(sc_objects[["CarPas_Wt"]], my_part = input$cluster1_cp1, cex = 1.5))
  output$plotclust_cl2_cp1 <- renderPlot(plotmap_m(sc_objects[["CarPas_Wt"]], my_part = input$cluster2_cp1, cex = 1.5))


  # Two clusters to compare
  cl1_cp1 <- reactive(names(sc_objects[["CarPas_Wt"]]@cpart[sc_objects[["CarPas_Wt"]]@cpart %in% input$cluster1_cp1]))
  cl2_cp1 <- reactive(names(sc_objects[["CarPas_Wt"]]@cpart[sc_objects[["CarPas_Wt"]]@cpart %in% input$cluster2_cp1]))

  # RaceID function to compare the two clusters

  diffexp_cp1 <- reactive( diffexpnb(getfdata(sc_objects[["CarPas_Wt"]],n=c(cl1_cp1(),cl2_cp1() )), A=cl1_cp1(), B=cl2_cp1() ) )
  p_val_cp1   <- reactive( input$c2c_pval_cp1)
  p_adj_cp1   <- reactive( input$padj_cp1)

  # Show two data tables of differentially expressed genes
  output$diffexpnb_2cl_1_cp1   <- renderDataTable(
    apply( diffexp_cp1()$res[diffexp_cp1()$res$pval <= p_val_cp1() & diffexp_cp1()$res$padj <= p_adj_cp1(),],
           2,
           function(x){round(x,digits = 8)}))

  output$diffexpnb_2cl_2_cp1   <- renderDataTable(
    apply( diffexp_cp1()$res[diffexp_cp1()$res$pval <= p_val_cp1() & diffexp_cp1()$res$padj <= p_adj_cp1(),],
           2,
           function(x){round(x,digits = 8)}))


  # Get input for top genes
  top_g_ma_cp1 <- reactive(input$top_genes_cp1)

  output$maplot_cp1   <- renderPlot({
    maplot(diffexp_cp1(),
           top = top_g_ma_cp1()
    )
  })


  output$maplot_cp1_2   <- renderPlotly({
    print(
      ggplotly(
        maplot(diffexp_cp1(), top = top_g_ma_cp1()) +
          aes(log2(mean), lfc, label=name,alpha=0.6) +
          geom_text(aes(label=name), data = . %>% filter(padj < 0.05), size=3, alpha=0.8) +
          xlab("Log2 mean expression") +
          ylab("Log2 fold change") +
          scale_y_continuous(breaks=seq(-20,20,0.5)) +
          scale_x_continuous(breaks=seq(-20,20,0.5))
      )
    )
  })












  ###
  ### CAR PaS from WT and KO
  ###
  output$plotmap_cp2 <- renderPlot({
    plotmap_m(sc_objects[["CarPas_WtKo"]], cex=1.5)
  })

  output$plotsymbolsmap_cp2 <- renderPlot({
    plotsymbolsmap_m(sc_mice_ob[["CarPas_WtKo"]], types = substring(colnames(sc_mice_ob[["CarPas_WtKo"]]@ndata), first = 1, last=4),
                     cex=1.5,
                     leg.pos = "topleft",
                     leg.cex = 1,
                     samples_col = c(c("brown",colorRampPalette(c("red", "gold"))(3)), c("blue",colorRampPalette(c("dodgerblue", "palegreen"))(3)) ))
  })


  ########################################################
  ### cdiff output
  ###
  ### Plot a particular cluster
  output$plotclust_cp2 <- renderPlot({
    plotmap_m(sc_objects[["CarPas_WtKo"]], my_part = input$cdiff_cl_cp2, cex = 1.5)
  })

  # Make cdiff list and user input reactive
  r_cdiff  <- reactive(cdiffs[["CarPas_WtKo"]])
  cd_p_val <- reactive(input$cdiff_pval_cp2)
  cd_p_adj <- reactive(input$cdiff_padj_cp2)


  # Show two data tables of cdiff output for one cluster
  output$upgenes_cp2   <- renderDataTable( apply( r_cdiff()[[input$cdiff_cl_cp2]][r_cdiff()[[input$cdiff_cl_cp2]]$pv <= cd_p_val() & r_cdiff()[[input$cdiff_cl_cp2]]$padj <= cd_p_adj(),],
                                              2,
                                              function(x){round(x,digits = 8)} ) )
  output$downgenes_cp2 <- renderDataTable( apply( r_cdiff()[[input$cdiff_cl_cp2]][r_cdiff()[[input$cdiff_cl_cp2]]$pv <= cd_p_val() & r_cdiff()[[input$cdiff_cl_cp2]]$padj <= cd_p_adj(),],
                                              2,
                                              function(x){round(x,digits = 8)} ) )


  ########################################################
  # plotexpmap
  ###
  ### Plotexpmap Gene expression in tSNE Map
  output$plotexptsne_cp2       <- renderPlot(plotexpmap_m(sc_objects[["CarPas_WtKo"]],input$ptexp_gene_cp2, logsc=F))
  output$plotexptsne_logsc_cp2 <- renderPlot(plotexpmap_m(sc_objects[["CarPas_WtKo"]],input$ptexp_gene_cp2, logsc=T, n = paste(input$ptexp_gene_cp2,"log2")))




  ########################################################
  ### Compare two clusters
  output$plotclust_cl1_cp2 <- renderPlot(plotmap_m(sc_objects[["CarPas_WtKo"]], my_part = input$cluster1_cp2, cex = 1.5))
  output$plotclust_cl2_cp2 <- renderPlot(plotmap_m(sc_objects[["CarPas_WtKo"]], my_part = input$cluster2_cp2, cex = 1.5))


  # Two clusters to compare
  cl1 <- reactive(names(sc_objects[["CarPas_WtKo"]]@cpart[sc_objects[["CarPas_WtKo"]]@cpart %in% input$cluster1_cp2]))
  cl2 <- reactive(names(sc_objects[["CarPas_WtKo"]]@cpart[sc_objects[["CarPas_WtKo"]]@cpart %in% input$cluster2_cp2]))

  # RaceID function to compare the two clusters

  diffexp <- reactive( diffexpnb(getfdata(sc_objects[["CarPas_WtKo"]],n=c(cl1(),cl2() )), A=cl1(), B=cl2() ) )
  p_val   <- reactive( input$c2c_pval_cp2)
  p_adj   <- reactive( input$padj_cp2)

  # Show two data tables of differentially expressed genes
  output$diffexpnb_2cl_1_cp2   <- renderDataTable(
    apply( diffexp()$res[diffexp()$res$pval <= p_val() & diffexp()$res$padj <= p_adj(),],
           2,
           function(x){round(x,digits = 8)}))

  output$diffexpnb_2cl_2_cp2   <- renderDataTable(
    apply( diffexp()$res[diffexp()$res$pval <= p_val() & diffexp()$res$padj <= p_adj(),],
           2,
           function(x){round(x,digits = 8)}))


  # Get input for top genes
    top_g_ma <- reactive(input$top_genes_cp2)

    output$maplot_cp2   <- renderPlot({
      maplot(diffexp(),
             top = top_g_ma()
      )
    })


    output$maplot_cp2_2   <- renderPlotly({
      print(
        ggplotly(
          maplot(diffexp(), top = top_g_ma()) +
            aes(log2(mean), lfc, label=name,alpha=0.6) +
            geom_text(aes(label=name), data = . %>% filter(padj < 0.05), size=3, alpha=0.8) +
            xlab("Log2 mean expression") +
            ylab("Log2 fold change") +
            scale_y_continuous(breaks=seq(-20,20,0.5)) +
            scale_x_continuous(breaks=seq(-20,20,0.5))
        )
      )
    })











    ###
    ### HSC from WT and KO
    ###


    output$plotmap_h2 <- renderPlot({
      plotmap_m(sc_objects[["HSC"]], cex=1.5)
    })

    output$plotsymbolsmap_h2 <- renderPlot({
      plotsymbolsmap_m(sc_mice_ob[["HSC"]], types = substring(colnames(sc_mice_ob[["HSC"]]@ndata), first = 1, last=4),
                       cex=1.5,
                       leg.pos = "bottomleft",
                       leg.cex = 1,
                       samples_col = c(c("brown",colorRampPalette(c("red", "gold"))(5)), c("blue",colorRampPalette(c("dodgerblue", "palegreen"))(4)) ))
    })


    ########################################################
    ### cdiff output
    ###
    ### Plot a particular cluster
    output$plotclust_h2 <- renderPlot({
      plotmap_m(sc_objects[["HSC"]], my_part = input$cdiff_cl_h2, cex = 1.5)
    })

    # Make cdiff list and user input reactive
    r_cdiff_h2  <- reactive(cdiffs[["HSC"]])
    cd_p_val_h2 <- reactive(input$cdiff_pval_h2)
    cd_p_adj_h2 <- reactive(input$cdiff_padj_h2)


    # Show two data tables of cdiff output for one cluster
    output$upgenes_h2   <- renderDataTable( apply( r_cdiff_h2()[[input$cdiff_cl_h2]][r_cdiff_h2()[[input$cdiff_cl_h2]]$pv <= cd_p_val_h2() & r_cdiff_h2()[[input$cdiff_cl_h2]]$padj <= cd_p_adj_h2(),],
                                                   2,
                                                   function(x){round(x,digits = 8)} ) )
    output$downgenes_h2 <- renderDataTable( apply( r_cdiff_h2()[[input$cdiff_cl_h2]][r_cdiff_h2()[[input$cdiff_cl_h2]]$pv <= cd_p_val_h2() & r_cdiff_h2()[[input$cdiff_cl_h2]]$padj <= cd_p_adj_h2(),],
                                                   2,
                                                   function(x){round(x,digits = 8)} ) )


    ########################################################
    # plotexpmap
    ###
    ### Plotexpmap Gene expression in tSNE Map
    output$plotexptsne_h2       <- renderPlot(plotexpmap_m(sc_objects[["HSC"]],input$ptexp_gene_h2, logsc=F))
    output$plotexptsne_logsc_h2 <- renderPlot(plotexpmap_m(sc_objects[["HSC"]],input$ptexp_gene_h2, logsc=T, n = paste(input$ptexp_gene_h2,"log2")))




    ########################################################
    ### Compare two clusters
    output$plotclust_cl1_h2 <- renderPlot(plotmap_m(sc_objects[["HSC"]], my_part = input$cluster1_h2, cex = 1.5))
    output$plotclust_cl2_h2 <- renderPlot(plotmap_m(sc_objects[["HSC"]], my_part = input$cluster2_h2, cex = 1.5))


    # Two clusters to compare
    cl1_h2 <- reactive(names(sc_objects[["HSC"]]@cpart[sc_objects[["HSC"]]@cpart %in% input$cluster1_h2]))
    cl2_h2 <- reactive(names(sc_objects[["HSC"]]@cpart[sc_objects[["HSC"]]@cpart %in% input$cluster2_h2]))

    # RaceID function to compare the two clusters

    diffexp_h2 <- reactive( diffexpnb(getfdata(sc_objects[["HSC"]],n=c(cl1_h2(),cl2_h2() )), A=cl1_h2(), B=cl2_h2() ) )
    p_val_h2   <- reactive( input$c2c_pval_h2)
    p_adj_h2   <- reactive( input$padj_h2)

    # Show two data tables of differentially expressed genes
    output$diffexpnb_2cl_1_h2   <- renderDataTable(
      apply( diffexp_h2()$res[diffexp_h2()$res$pval <= p_val_h2() & diffexp_h2()$res$padj <= p_adj_h2(),],
             2,
             function(x){round(x,digits = 8)}))

    output$diffexpnb_2cl_2_h2   <- renderDataTable(
      apply( diffexp_h2()$res[diffexp_h2()$res$pval <= p_val_h2() & diffexp_h2()$res$padj <= p_adj_h2(),],
             2,
             function(x){round(x,digits = 8)}))


    # Get input for top genes
    top_g_ma_h2 <- reactive(input$top_genes_h2)

    output$maplot_h2   <- renderPlot({
      maplot(diffexp_h2(),
             top = top_g_ma_h2()
      )
    })


    output$maplot_h2_2   <- renderPlotly({
      print(
        ggplotly(
          maplot(diffexp_h2(), top = top_g_ma_h2()) +
            aes(log2(mean), lfc, label=name,alpha=0.6) +
            geom_text(aes(label=name), data = . %>% filter(padj < 0.05), size=3, alpha=0.8) +
            xlab("Log2 mean expression") +
            ylab("Log2 fold change") +
            scale_y_continuous(breaks=seq(-20,20,0.5)) +
            scale_x_continuous(breaks=seq(-20,20,0.5))
        )
      )
    })














    ###
    ### HSC & LSK from WT and KO
    ###


    output$plotmap_h1 <- renderPlot({
      plotmap_m(sc_objects[["HSC_LSK"]], cex=1.5)
    })

    output$plotsymbolsmap_h1 <- renderPlot({
      plotsymbolsmap_m(sc_mice_ob[["HSC_LSK"]], types = substring(colnames(sc_mice_ob[["HSC_LSK"]]@ndata), first = 1, last=4),
                       cex=1.5,
                       leg.pos = "bottomleft",
                       leg.cex = 1,
                       samples_col = c(c("brown",colorRampPalette(c("red", "gold"))(5)), c("blue",colorRampPalette(c("dodgerblue", "palegreen"))(4)) ))
    })


    ########################################################
    ### cdiff output
    ###
    ### Plot a particular cluster
    output$plotclust_h1 <- renderPlot({
      plotmap_m(sc_objects[["HSC_LSK"]], my_part = input$cdiff_cl_h1, cex = 1.5)
    })

    # Make cdiff list and user input reactive
    r_cdiff_h1  <- reactive(cdiffs[["HSC_LSK"]])
    cd_p_val_h1 <- reactive(input$cdiff_pval_h1)
    cd_p_adj_h1 <- reactive(input$cdiff_padj_h1)


    # Show two data tables of cdiff output for one cluster
    output$upgenes_h1   <- renderDataTable( apply( r_cdiff_h1()[[input$cdiff_cl_h1]][r_cdiff_h1()[[input$cdiff_cl_h1]]$pv <= cd_p_val_h1() & r_cdiff_h1()[[input$cdiff_cl_h1]]$padj <= cd_p_adj_h1(),],
                                                    2,
                                                    function(x){round(x,digits = 8)} ) )
    output$downgenes_h1 <- renderDataTable( apply( r_cdiff_h1()[[input$cdiff_cl_h1]][r_cdiff_h1()[[input$cdiff_cl_h1]]$pv <= cd_p_val_h1() & r_cdiff_h1()[[input$cdiff_cl_h1]]$padj <= cd_p_adj_h1(),],
                                                    2,
                                                    function(x){round(x,digits = 8)} ) )


    ########################################################
    # plotexpmap
    ###
    ### Plotexpmap Gene expression in tSNE Map
    output$plotexptsne_h1       <- renderPlot(plotexpmap_m(sc_objects[["HSC_LSK"]],input$ptexp_gene_h1, logsc=F))
    output$plotexptsne_logsc_h1 <- renderPlot(plotexpmap_m(sc_objects[["HSC_LSK"]],input$ptexp_gene_h1, logsc=T, n = paste(input$ptexp_gene_h1,"log2")))




    ########################################################
    ### Compare two clusters
    output$plotclust_cl1_h1 <- renderPlot(plotmap_m(sc_objects[["HSC_LSK"]], my_part = input$cluster1_h1, cex = 1.5))
    output$plotclust_cl2_h1 <- renderPlot(plotmap_m(sc_objects[["HSC_LSK"]], my_part = input$cluster2_h1, cex = 1.5))


    # Two clusters to compare
    cl1_h1 <- reactive(names(sc_objects[["HSC_LSK"]]@cpart[sc_objects[["HSC_LSK"]]@cpart %in% input$cluster1_h1]))
    cl2_h1 <- reactive(names(sc_objects[["HSC_LSK"]]@cpart[sc_objects[["HSC_LSK"]]@cpart %in% input$cluster2_h1]))

    # RaceID function to compare the two clusters

    diffexp_h1 <- reactive( diffexpnb(getfdata(sc_objects[["HSC_LSK"]],n=c(cl1_h1(),cl2_h1() )), A=cl1_h1(), B=cl2_h1() ) )
    p_val_h1   <- reactive( input$c2c_pval_h1)
    p_adj_h1   <- reactive( input$padj_h1)

    # Show two data tables of differentially expressed genes
    output$diffexpnb_2cl_1_h1   <- renderDataTable(
      apply( diffexp_h1()$res[diffexp_h1()$res$pval <= p_val_h1() & diffexp_h1()$res$padj <= p_adj_h1(),],
             2,
             function(x){round(x,digits = 8)}))

    output$diffexpnb_2cl_2_h1   <- renderDataTable(
      apply( diffexp_h1()$res[diffexp_h1()$res$pval <= p_val_h1() & diffexp_h1()$res$padj <= p_adj_h1(),],
             2,
             function(x){round(x,digits = 8)}))


    # Get input for top genes
    top_g_ma_h1 <- reactive(input$top_genes_h1)

    output$maplot_h1   <- renderPlot({
      maplot(diffexp_h1(),
             top = top_g_ma_h1()
      )
    })


    output$maplot_h1_2   <- renderPlotly({
      print(
        ggplotly(
          maplot(diffexp_h1(), top = top_g_ma_h1()) +
            aes(log2(mean), lfc, label=name,alpha=0.6) +
            geom_text(aes(label=name), data = . %>% filter(padj < 0.05), size=3, alpha=0.8) +
            xlab("Log2 mean expression") +
            ylab("Log2 fold change") +
            scale_y_continuous(breaks=seq(-20,20,0.5)) +
            scale_x_continuous(breaks=seq(-20,20,0.5))
        )
      )
    })
    
    
    
    

  
})
