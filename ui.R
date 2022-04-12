# shiny Server

library(shiny)
library(plotly)
genenames <- readRDS("./genenames.Rds")
shinyUI(navbarPage("Role of Ebf1 in the hematopoietic bone marrow niche", position="fixed-top",
                   
                   
                   ###################################
                   ###################################
                   ###################################
                   # CAR & PaS cells (Wt)
                   
                  
                   tabPanel("CAR & PaS cells (Wt)",
                            br(),
                            br(),
                            br(),
                            a(href="https://imprint.ie-freiburg.mpg.de/imprint/imprint.php?ref=12", "Imprint"),
                            br(),
                            a(href="https://imprint.ie-freiburg.mpg.de/privacy/privacy.php?ref=12", "Privacy Notice"),
                            br(),
                            
                            p("Shiny app contact:", style = "font-size:18px", align="left",
                              a(href  = "mailto:herman@ie-freiburg.mpg.de", "Josip S. Herman",
                                style = "font-size:18px")),
                            br(),
                            
                            fluidRow(
                              column(6,
                                     h1("EBF1-deficient bone marrow stroma elicits persistent changes in HSC potential"),
                                     a(href  = "https://doi.org/10.1038/s41590-020-0595-7", "https://doi.org/10.1038/s41590-020-0595-7",
                                       style = "font-size:18px"),
                                     p(
                                       a(href  = "mailto:derecka@ie-freiburg.mpg.de", "Marta Derecka",
                                         style = "font-size:20px"),
                                       style = "font-size:20px",
                                       ", Josip Stefan Herman, Pierre Cauchy, Senthilkumar Ramamoorthy, Ekaterina Lupar, Dominic Grün, ",
                                       a(href  = "mailto:grosschedl@ie-freiburg.mpg.de", "Rudolf Grosschedl",
                                         style = "font-size:20px")
                                     )
                              )),
                            
                            h2("Identified clusters after k-medoids clustering & used sorting gates"),
                            h5("(Loading of data takes app. 30s)"),
                            
                            
                            fluidRow(
                              column(6,
                                     plotOutput("plotmap_cp1",       
                                                width = "100%",
                                                height = "600px")),
                              column(6,
                                     plotOutput("plotsymbolsmap_cp1",
                                                width = "100%",
                                                height = "600px"))
                            ),
                            
                            
                            br(),
                            h2("Compare one cluster against all the other clusters"),
                            
                            fluidRow(
                              column(2,
                                     numericInput("cdiff_cl_cp1",
                                                  value=1,
                                                  label = "Choose a cluster", 
                                                  min = 1, 
                                                  max = 1000  )),
                              column(5,
                                     plotOutput("plotclust_cp1",
                                                width = "100%",
                                                height = "500px"), 
                                     offset = 2)
                              
                            ),
                            
                            h5("mean.ncl = mean expression of gene not in cluster"),
                            h5("mean.cl  = mean expression of gene in cluster"),
                            h5("fc  = fold change compared to all other clusters"),
                            h5("pv = p-value"),
                            h5("padj = adjusted p-value"),
                            br(),
                            
                            fluidRow(
                              column(2, 
                                     numericInput("cdiff_pval_cp1", 
                                                  value = 0.05, 
                                                  "Set p-value threshhold"),
                                     numericInput("cdiff_padj_cp1",
                                                  value = 1, 
                                                  "Set adjusted p-value threshhold"))
                            ),
                            fluidRow(
                              column(6,
                                     DT::dataTableOutput("upgenes_cp1")),
                              column(6,
                                     DT::dataTableOutput("downgenes_cp1"))
                            ),
                            
                            br(),
                            br(),
                            br(),
                            
                            h2("Plot gene expression in tSNE-Map"),
                            fluidRow(
                              column(2,
                                     selectInput("ptexp_gene_cp1", 
                                                 label = "Enter a genesymbol",
                                                 genenames,
                                                 selected = "Cxcl12")),
                              column(5,
                                     plotOutput("plotexptsne_cp1",       
                                                width = "100%",
                                                height = "500px")),
                              column(5,
                                     plotOutput("plotexptsne_logsc_cp1", 
                                                width = "100%",
                                                height = "500px"))
                            ),
                            
                            
                            br(),
                            br(),
                            br(),
                            
                            
                            h2("Compare gene expression of two clusters"),
                            
                            fluidRow(column(2, 
                                            numericInput("cluster1_cp1", 
                                                         value = 1, 
                                                         "Choose cluster 1", 
                                                         min = 1, 
                                                         max = 44 ),
                                            numericInput("cluster2_cp1", 
                                                         value = 2, 
                                                         "Choose cluster 2", 
                                                         min = 1, 
                                                         max = 44 ))
                            ),
                            
                            fluidRow(
                              column(6,
                                     plotOutput("plotclust_cl1_cp1",
                                                width = "100%",
                                                height = "600px")),
                              column(6,
                                     plotOutput("plotclust_cl2_cp1",
                                                width = "100%",
                                                height = "600px"))
                            ),
                            
                            
                            fluidRow(
                              column(2, 
                                     numericInput("c2c_pval_cp1", 
                                                  value = 0.05, 
                                                  "Set p-value threshhold"),
                                     numericInput("padj_cp1",
                                                  value = 1, 
                                                  "Set adjusted p-value threshhold"))
                            ),
                            fluidRow(
                              column(6,
                                     DT::dataTableOutput("diffexpnb_2cl_1_cp1")),
                              column(6,
                                     DT::dataTableOutput("diffexpnb_2cl_2_cp1"))
                            ),
                            
                            
                            br(),
                            br(),
                            br(),
                            
                            fluidRow(
                              column(2,
                                     numericInput("top_genes_cp1", 
                                                  value = 100, 
                                                  "Choose number of top genes to display"))
                            ),
                            
                            # fluidRow(
                            #   column(2,
                            #          selectInput("maplot_iactiv_cp2", 
                            #                        label = "Choose MA plot mode",
                            #                        c("interactive", "non-interactive"),
                            #                        selected = "non-interactive"))
                            # ),
                            
                            br(),
                            h2("MA Plot"),
                            fluidRow(
                              column(12, 
                                     plotOutput("maplot_cp1", 
                                                width = "100%",
                                                height = "1000px"))
                            ),
                            
                            br(),
                            h2("Interactive MA Plot"),
                            fluidRow(
                              column(12, 
                                     plotlyOutput("maplot_cp1_2", 
                                                  width = "100%",
                                                  height = "1000px"))
                            ),
                            br(),
                            br(),
                            br()
                   ),
                   
                   
                   
                   
                   
                   
                   ###################################
                   ###################################
                   ###################################
                   # CAR & PaS cells (Wt & Ko)
                   
                   
                   
                   tabPanel("CAR & PaS cells (Wt & Ko)",
                            br(),
                            br(),
                            br(),
                            a(href="https://imprint.ie-freiburg.mpg.de/imprint/imprint.php?ref=12", "Imprint"),
                            br(),
                            a(href="https://imprint.ie-freiburg.mpg.de/privacy/privacy.php?ref=12", "Privacy Notice"),
                            br(),
                            
                            p("Shiny app contact:", style = "font-size:18px", align="left",
                              a(href  = "mailto:herman@ie-freiburg.mpg.de", "Josip S. Herman",
                                style = "font-size:18px")),
                            br(),
                            
                            fluidRow(
                              column(6,
                                     h1("EBF1-deficient bone marrow stroma elicits persistent changes in HSC potential"),
                                     a(href  = "https://doi.org/10.1038/s41590-020-0595-7", "https://doi.org/10.1038/s41590-020-0595-7",
                                       style = "font-size:18px"),
                                     p(
                                       a(href  = "mailto:derecka@ie-freiburg.mpg.de", "Marta Derecka",
                                         style = "font-size:20px"),
                                       style = "font-size:20px",
                                       ", Josip Stefan Herman, Pierre Cauchy, Senthilkumar Ramamoorthy, Ekaterina Lupar, Dominic Grün, ",
                                       a(href  = "mailto:grosschedl@ie-freiburg.mpg.de", "Rudolf Grosschedl",
                                         style = "font-size:20px")
                                     )
                              )),
                            
                            h2("Identified clusters after k-medoids clustering & used sorting gates"),
                            h5("(Loading of data takes app. 30s)"),
                            
                            fluidRow(
                              column(6,
                                     plotOutput("plotmap_cp2",       
                                                width = "100%",
                                                height = "600px")),
                              column(6,
                                     plotOutput("plotsymbolsmap_cp2",
                                                width = "100%",
                                                height = "600px"))
                              ),


                            br(),
                            h2("Compare one cluster against all the other clusters"),

                            fluidRow(
                              column(2,
                                     numericInput("cdiff_cl_cp2",
                                                  value=1,
                                                  label = "Choose a cluster", 
                                                  min = 1, 
                                                  max = 1000  )),
                              column(5,
                                     plotOutput("plotclust_cp2",
                                                  width = "100%",
                                                  height = "500px"), 
                                     offset = 2)
                              
                            ),

                            h5("mean.ncl = mean expression of gene not in cluster"),
                            h5("mean.cl  = mean expression of gene in cluster"),
                            h5("fc  = fold change compared to all other clusters"),
                            h5("pv = p-value"),
                            h5("padj = adjusted p-value"),
                            br(),
                            
                            fluidRow(
                              column(2, 
                                     numericInput("cdiff_pval_cp2", 
                                                  value = 0.05, 
                                                  "Set p-value threshhold"),
                                     numericInput("cdiff_padj_cp2",
                                                  value = 1, 
                                                  "Set adjusted p-value threshhold"))
                            ),
                            fluidRow(
                              column(6,
                                     DT::dataTableOutput("upgenes_cp2")),
                              column(6,
                                     DT::dataTableOutput("downgenes_cp2"))
                            ),

                            br(),
                            br(),
                            br(),
                            
                            h2("Plot gene expression in tSNE-Map"),
                            fluidRow(
                              column(2,
                                     selectInput("ptexp_gene_cp2", 
                                                 label = "Enter a genesymbol",
                                                 genenames,
                                                 selected = "Cxcl12")),
                              column(5,
                                     plotOutput("plotexptsne_cp2",       
                                                  width = "100%",
                                                  height = "500px")),
                              column(5,
                                     plotOutput("plotexptsne_logsc_cp2", 
                                                width = "100%",
                                                height = "500px"))
                            ),
                            
                            
                            br(),
                            br(),
                            br(),
                            

                            h2("Compare gene expression of two clusters"),

                            fluidRow(column(2, 
                                            numericInput("cluster1_cp2", 
                                                         value = 1, 
                                                         "Choose cluster 1", 
                                                         min = 1, 
                                                         max = 44 ),
                                            numericInput("cluster2_cp2", 
                                                         value = 2, 
                                                         "Choose cluster 2", 
                                                         min = 1, 
                                                         max = 44 ))
                            ),
                            
                            fluidRow(
                              column(6,
                                     plotOutput("plotclust_cl1_cp2",
                                                  width = "100%",
                                                height = "600px")),
                              column(6,
                                     plotOutput("plotclust_cl2_cp2",
                                                width = "100%",
                                                height = "600px"))
                            ),


                            fluidRow(
                              column(2, 
                                     numericInput("c2c_pval_cp2", 
                                                  value = 0.05, 
                                                  "Set p-value threshhold"),
                                     numericInput("padj_cp2",
                                                  value = 1, 
                                                  "Set adjusted p-value threshhold"))
                              ),
                            fluidRow(
                              column(6,
                                     DT::dataTableOutput("diffexpnb_2cl_1_cp2")),
                              column(6,
                                     DT::dataTableOutput("diffexpnb_2cl_2_cp2"))
                            ),
                            
                            
                            br(),
                            br(),
                            br(),
                            
                            fluidRow(
                              column(2,
                                     numericInput("top_genes_cp2", 
                                                  value = 100, 
                                                  "Choose number of top genes to display"))
                            ),
                            
                            # fluidRow(
                            #   column(2,
                            #          selectInput("maplot_iactiv_cp2", 
                            #                        label = "Choose MA plot mode",
                            #                        c("interactive", "non-interactive"),
                            #                        selected = "non-interactive"))
                            # ),
                            
                            br(),
                            h2("MA Plot"),
                            fluidRow(
                              column(12, 
                                     plotOutput("maplot_cp2", 
                                                width = "100%",
                                                height = "1000px"))
                            ),
                            
                            br(),
                            h2("Interactive MA Plot"),
                            fluidRow(
                              column(12, 
                                     plotlyOutput("maplot_cp2_2", 
                                                width = "100%",
                                                height = "1000px"))
                            ),
                            br(),
                            br(),
                            br()#,
                            
                            
                            
                            
                            # Within cluster comparison of Wt vs Ko
                            # br(),
                            # h2("Compare Wt vs. Ko within one cluster"),
                            # 
                            # fluidRow(
                            #   column(2,
                            #          numericInput("wth_cl",
                            #                       value=1,
                            #                       label = "Choose a cluster", 
                            #                       min = 1, 
                            #                       max = 1000  )),
                            #   column(5,
                            #          plotOutput("plotclust_wth_cp2",
                            #                     width = "100%",
                            #                     height = "500px"), 
                            #          offset = 2)
                            #   
                            # ),
                            # 
                            # h5("mean.ncl = mean expression of gene not in cluster"),
                            # h5("mean.cl  = mean expression of gene in cluster"),
                            # h5("fc  = fold change compared to all other clusters"),
                            # h5("pv = p-value"),
                            # h5("padj = adjusted p-value"),
                            # br(),
                            # 
                            # fluidRow(
                            #   column(2, 
                            #          numericInput("wth_pval_cp2", 
                            #                       value = 0.05, 
                            #                       "Set p-value threshhold"),
                            #          numericInput("wth_padj_cp2",
                            #                       value = 1, 
                            #                       "Set adjusted p-value threshhold"))
                            # ),
                            # fluidRow(
                            #   column(6,
                            #          DT::dataTableOutput("wth_upgenes_cp2")),
                            #   column(6,
                            #          DT::dataTableOutput("wth_downgenes_cp2"))
                            # )
                            
                          ),
                   
                   
                   
                   
                   
                   ###################################
                   ###################################
                   ###################################
                   # HSC & LSK (Wt & Ko)
                   
                   
                   tabPanel("HSC (Wt & Ko)",
                            br(),
                            br(),
                            br(),
                            a(href="https://imprint.ie-freiburg.mpg.de/imprint/imprint.php?ref=12", "Imprint"),
                            br(),
                            a(href="https://imprint.ie-freiburg.mpg.de/privacy/privacy.php?ref=12", "Privacy Notice"),
                            br(),
                            
                            p("Shiny app contact:", style = "font-size:18px", align="left",
                              a(href  = "mailto:herman@ie-freiburg.mpg.de", "Josip S. Herman",
                                style = "font-size:18px")),
                            br(),
                            
                            fluidRow(
                              column(6,
                                     h1("EBF1-deficient bone marrow stroma elicits persistent changes in HSC potential"),
                                     a(href  = "https://doi.org/10.1038/s41590-020-0595-7", "https://doi.org/10.1038/s41590-020-0595-7",
                                       style = "font-size:18px"),
                                     p(
                                       a(href  = "mailto:derecka@ie-freiburg.mpg.de", "Marta Derecka",
                                         style = "font-size:20px"),
                                       style = "font-size:20px",
                                       ", Josip Stefan Herman, Pierre Cauchy, Senthilkumar Ramamoorthy, Ekaterina Lupar, Dominic Grün, ",
                                       a(href  = "mailto:grosschedl@ie-freiburg.mpg.de", "Rudolf Grosschedl",
                                         style = "font-size:20px")
                                     )
                              )),
                            
                            h2("Identified clusters after k-medoids clustering & used sorting gates"),
                            h5("(Loading of data takes app. 30s)"),
                            
                            
                            fluidRow(
                              column(6,
                                     plotOutput("plotmap_h2",       
                                                width = "100%",
                                                height = "600px")),
                              column(6,
                                     plotOutput("plotsymbolsmap_h2",
                                                width = "100%",
                                                height = "600px"))
                            ),
                            
                            
                            br(),
                            h2("Compare one cluster against all the other clusters"),
                            
                            fluidRow(
                              column(2,
                                     numericInput("cdiff_cl_h2",
                                                  value=1,
                                                  label = "Choose a cluster", 
                                                  min = 1, 
                                                  max = 1000  )),
                              column(5,
                                     plotOutput("plotclust_h2",
                                                width = "100%",
                                                height = "500px"), 
                                     offset = 2)
                              
                            ),
                            
                            h5("mean.ncl = mean expression of gene not in cluster"),
                            h5("mean.cl  = mean expression of gene in cluster"),
                            h5("fc  = fold change compared to all other clusters"),
                            h5("pv = p-value"),
                            h5("padj = adjusted p-value"),
                            br(),
                            
                            fluidRow(
                              column(2, 
                                     numericInput("cdiff_pval_h2", 
                                                  value = 0.05, 
                                                  "Set p-value threshhold"),
                                     numericInput("cdiff_padj_h2",
                                                  value = 1, 
                                                  "Set adjusted p-value threshhold"))
                            ),
                            fluidRow(
                              column(6,
                                     DT::dataTableOutput("upgenes_h2")),
                              column(6,
                                     DT::dataTableOutput("downgenes_h2"))
                            ),
                            
                            br(),
                            br(),
                            br(),
                            
                            h2("Plot gene expression in tSNE-Map"),
                            fluidRow(
                              column(2,
                                     selectInput("ptexp_gene_h2", 
                                                 label = "Enter a genesymbol",
                                                 genenames,
                                                 selected = "Cxcl12")),
                              column(5,
                                     plotOutput("plotexptsne_h2",       
                                                width = "100%",
                                                height = "500px")),
                              column(5,
                                     plotOutput("plotexptsne_logsc_h2", 
                                                width = "100%",
                                                height = "500px"))
                            ),
                            
                            
                            br(),
                            br(),
                            br(),
                            
                            
                            h2("Compare gene expression of two clusters"),
                            
                            fluidRow(column(2, 
                                            numericInput("cluster1_h2", 
                                                         value = 1, 
                                                         "Choose cluster 1", 
                                                         min = 1, 
                                                         max = 44 ),
                                            numericInput("cluster2_h2", 
                                                         value = 2, 
                                                         "Choose cluster 2", 
                                                         min = 1, 
                                                         max = 44 ))
                            ),
                            
                            fluidRow(
                              column(6,
                                     plotOutput("plotclust_cl1_h2",
                                                width = "100%",
                                                height = "600px")),
                              column(6,
                                     plotOutput("plotclust_cl2_h2",
                                                width = "100%",
                                                height = "600px"))
                            ),
                            
                            
                            fluidRow(
                              column(2, 
                                     numericInput("c2c_pval_h2", 
                                                  value = 0.05, 
                                                  "Set p-value threshhold"),
                                     numericInput("padj_h2",
                                                  value = 1, 
                                                  "Set adjusted p-value threshhold"))
                            ),
                            fluidRow(
                              column(6,
                                     DT::dataTableOutput("diffexpnb_2cl_1_h2")),
                              column(6,
                                     DT::dataTableOutput("diffexpnb_2cl_2_h2"))
                            ),
                            
                            
                            br(),
                            br(),
                            br(),
                            
                            fluidRow(
                              column(2,
                                     numericInput("top_genes_h2", 
                                                  value = 100, 
                                                  "Choose number of top genes to display"))
                            ),
                            
                            # fluidRow(
                            #   column(2,
                            #          selectInput("maplot_iactiv_cp2", 
                            #                        label = "Choose MA plot mode",
                            #                        c("interactive", "non-interactive"),
                            #                        selected = "non-interactive"))
                            # ),
                            
                            br(),
                            h2("MA Plot"),
                            fluidRow(
                              column(12, 
                                     plotOutput("maplot_h2", 
                                                width = "100%",
                                                height = "1000px"))
                            ),
                            
                            br(),
                            h2("Interactive MA Plot"),
                            fluidRow(
                              column(12, 
                                     plotlyOutput("maplot_h2_2", 
                                                  width = "100%",
                                                  height = "1000px"))
                            ),
                            br(),
                            br(),
                            br()
                            
                            
                            
                   ),
                   
                   
                   
                   
                   
                   
                   
                   ###################################
                   ###################################
                   ###################################
                   # HSC (Wt & Ko) 
                   
                   tabPanel("HSC & LSK (Wt & Ko)",
                            br(),
                            br(),
                            br(),
                            a(href="https://imprint.ie-freiburg.mpg.de/imprint/imprint.php?ref=12", "Imprint"),
                            br(),
                            a(href="https://imprint.ie-freiburg.mpg.de/privacy/privacy.php?ref=12", "Privacy Notice"),
                            br(),
                            
                            p("Shiny app contact:", style = "font-size:18px", align="left",
                              a(href  = "mailto:herman@ie-freiburg.mpg.de", "Josip S. Herman",
                                style = "font-size:18px")),
                            br(),
                            
                            fluidRow(
                              column(6,
                                     h1("EBF1-deficient bone marrow stroma elicits persistent changes in HSC potential"),
                                     a(href  = "https://doi.org/10.1038/s41590-020-0595-7", "https://doi.org/10.1038/s41590-020-0595-7",
                                       style = "font-size:18px"),
                                     p(
                                       a(href  = "mailto:derecka@ie-freiburg.mpg.de", "Marta Derecka",
                                         style = "font-size:20px"),
                                       style = "font-size:20px",
                                       ", Josip Stefan Herman, Pierre Cauchy, Senthilkumar Ramamoorthy, Ekaterina Lupar, Dominic Grün, ",
                                       a(href  = "mailto:grosschedl@ie-freiburg.mpg.de", "Rudolf Grosschedl",
                                         style = "font-size:20px")
                                     )
                              )),
                            
                            h2("Identified clusters after k-medoids clustering & used sorting gates"),
                            h5("(Loading of data takes app. 30s)"),
                            
                            
                            fluidRow(
                              column(6,
                                     plotOutput("plotmap_h1",       
                                                width = "100%",
                                                height = "600px")),
                              column(6,
                                     plotOutput("plotsymbolsmap_h1",
                                                width = "100%",
                                                height = "600px"))
                            ),
                            
                            
                            br(),
                            h2("Compare one cluster against all the other clusters"),
                            
                            fluidRow(
                              column(2,
                                     numericInput("cdiff_cl_h1",
                                                  value=1,
                                                  label = "Choose a cluster", 
                                                  min = 1, 
                                                  max = 1000  )),
                              column(5,
                                     plotOutput("plotclust_h1",
                                                width = "100%",
                                                height = "500px"), 
                                     offset = 2)
                              
                            ),
                            
                            h5("mean.ncl = mean expression of gene not in cluster"),
                            h5("mean.cl  = mean expression of gene in cluster"),
                            h5("fc  = fold change compared to all other clusters"),
                            h5("pv = p-value"),
                            h5("padj = adjusted p-value"),
                            br(),
                            
                            fluidRow(
                              column(2, 
                                     numericInput("cdiff_pval_h1", 
                                                  value = 0.05, 
                                                  "Set p-value threshhold"),
                                     numericInput("cdiff_padj_h1",
                                                  value = 1, 
                                                  "Set adjusted p-value threshhold"))
                            ),
                            fluidRow(
                              column(6,
                                     DT::dataTableOutput("upgenes_h1")),
                              column(6,
                                     DT::dataTableOutput("downgenes_h1"))
                            ),
                            
                            br(),
                            br(),
                            br(),
                            
                            h2("Plot gene expression in tSNE-Map"),
                            fluidRow(
                              column(2,
                                     selectInput("ptexp_gene_h1", 
                                                 label = "Enter a genesymbol",
                                                 genenames,
                                                 selected = "Cxcl12")),
                              column(5,
                                     plotOutput("plotexptsne_h1",       
                                                width = "100%",
                                                height = "500px")),
                              column(5,
                                     plotOutput("plotexptsne_logsc_h1", 
                                                width = "100%",
                                                height = "500px"))
                            ),
                            
                            
                            br(),
                            br(),
                            br(),
                            
                            
                            h2("Compare gene expression of two clusters"),
                            
                            fluidRow(column(2, 
                                            numericInput("cluster1_h1", 
                                                         value = 1, 
                                                         "Choose cluster 1", 
                                                         min = 1, 
                                                         max = 44 ),
                                            numericInput("cluster2_h1", 
                                                         value = 2, 
                                                         "Choose cluster 2", 
                                                         min = 1, 
                                                         max = 44 ))
                            ),
                            
                            fluidRow(
                              column(6,
                                     plotOutput("plotclust_cl1_h1",
                                                width = "100%",
                                                height = "600px")),
                              column(6,
                                     plotOutput("plotclust_cl2_h1",
                                                width = "100%",
                                                height = "600px"))
                            ),
                            
                            
                            fluidRow(
                              column(2, 
                                     numericInput("c2c_pval_h1", 
                                                  value = 0.05, 
                                                  "Set p-value threshhold"),
                                     numericInput("padj_h1",
                                                  value = 1, 
                                                  "Set adjusted p-value threshhold"))
                            ),
                            fluidRow(
                              column(6,
                                     DT::dataTableOutput("diffexpnb_2cl_1_h1")),
                              column(6,
                                     DT::dataTableOutput("diffexpnb_2cl_2_h1"))
                            ),
                            
                            
                            br(),
                            br(),
                            br(),
                            
                            fluidRow(
                              column(2,
                                     numericInput("top_genes_h1", 
                                                  value = 100, 
                                                  "Choose number of top genes to display"))
                            ),
                            
                            # fluidRow(
                            #   column(2,
                            #          selectInput("maplot_iactiv_cp2", 
                            #                        label = "Choose MA plot mode",
                            #                        c("interactive", "non-interactive"),
                            #                        selected = "non-interactive"))
                            # ),
                            
                            br(),
                            h2("MA Plot"),
                            fluidRow(
                              column(12, 
                                     plotOutput("maplot_h1", 
                                                width = "100%",
                                                height = "1000px"))
                            ),
                            
                            br(),
                            h2("Interactive MA Plot"),
                            fluidRow(
                              column(12, 
                                     plotlyOutput("maplot_h1_2", 
                                                  width = "100%",
                                                  height = "1000px"))
                            ),
                            br(),
                            br(),
                            br()
                   )
                   
                   
))
