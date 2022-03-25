#Author : Ayushi Pathak
#Date : 19 March 2022

library(shiny)
#install.packages("shinydashboard")
library(shinydashboard)

ui<- dashboardPage(
  skin='black',
  title='SIMULATION-Gene',
  
  #HEADER----------------------------------------------------------------
  
  dashboardHeader(
    title = 'SIMULATION-Gene',
    dropdownMenu(
      type='notifications',
      headerText = strong('HELP'),
      icon = icon('question'),
      badgeStatus = NULL,
      notificationItem(
        text = 'Creates data frames based on models',
         
       
        icon = icon('spinner')
      )
    ),
    
    tags$li(
      tags$a(
        strong('ABOUT Simulation-G'),
        height=40,
        href='https://github.com/AyushiPathak',
        title = "",
        target = "_blank"
      ),
      class='dropdown'
    )
  ),
  
  #SIDEBAR----------------------------------------------------------------
  
  dashboardSidebar(
    sidebarMenu(
      menuSubItem( "MARKOV MODEL" , tabName = 'mm', href = NULL, newtab = TRUE,
                   icon = shiny::icon("arrow-right"), selected = NULL),
      menuSubItem('HARDY WEIBERG MODEL',tabName = 'hw',href = NULL, newtab = TRUE,
                  icon = shiny::icon("arrow-right"), selected = NULL),
      menuSubItem('MIGRATION MODELS', tabName = 'migrations', href = NULL, newtab = TRUE,
                  icon = shiny::icon("arrow-right"), selected = NULL)
      
    )
  ),
  
  #BODY----------------------------------------------------------------
  
  dashboardBody(
    tabItems(
      tabItem('mm',
              fluidRow(
                tabBox(width = 12,
                  tabPanel(title = 'GRAPH',
                           numericInput(inputId='base',
                                        label='Choose the number of baspairs',
                                        value = "0", 
                                        width = NULL,
                                        min = NA,
                                        max = NA, 
                                        step = NA ),
                           numericInput(inputId='propbability_a',
                                        label='Choose the probability of A',
                                        value = "0", 
                                        width = NULL,
                                        min = NA,
                                        max = NA, 
                                        step = NA),
                           numericInput(inputId='propbability_t',
                                        label='Choose the probability of T',
                                        value = "0", 
                                        width = NULL,
                                        min = NA,
                                        max = NA, 
                                        step = NA),
                           numericInput(inputId='propbability_g',
                                        label='Choose the probability of G',
                                        value = "0", 
                                        width = NULL,
                                        min = NA,
                                        max = NA, 
                                        step = NA),
                           numericInput(inputId='propbability_c',
                                        label='Choose the probability of C',
                                        value = "0", 
                                        width = NULL,
                                        min = NA,
                                        max = NA, 
                                        step = NA),
                           plotOutput('bargraph'),
                           tableOutput('bargraph_1'),
                           tableOutput('bargraph_2'),
                           tableOutput('bargraph_3'),
                           tableOutput('bargraph_4')
                           # downloadButton('downloadData','Download Data'),
                           # downloadButton('downloadPlot','Download Plot')
                           
                           
                           )
                  
                )
              )
              
        
      ),
      tabItem('hw',
              tabBox(width = 12,
                tabPanel(title='Graph',
                         numericInput(inputId='difference',
                                      label='Enter the sequence gap value',
                                      value = ".01", 
                                      width = NULL,
                                      min = NA,
                                      max = NA, 
                                      step = NA ),
                         numericInput(inputId='initial_value',
                                      label='Enter the initial value in range 0-1',
                                      value = "0", 
                                      width = NULL,
                                      min = NA,
                                      max = NA, 
                                      step = NA ),
                         numericInput(inputId='final_value',
                                      label='Enter the final value (0<initail<final<1',
                                      value = "0", 
                                      width = NULL,
                                      min = NA,
                                      max = NA, 
                                      step = NA ),
                         plotOutput('hw_out_plot'),
                         tableOutput('hw_out_plot_1')
                         )
              )
        
      ),
      tabItem('migrations',
              fluidRow(
                tabBox(width = 12,
                  tabPanel(title = 'Island Model',
                           numericInput(inputId='migration_rate',
                                        label='Enter the migration rate',
                                        value = ".01", 
                                        width = NULL,
                                        min = NA,
                                        max = NA, 
                                        step = NA ),
                           numericInput(inputId='number_of_generations',
                                        label='Choose the Number of Generations',
                                        value = "0", 
                                        width = NULL,
                                        min = NA,
                                        max = NA, 
                                        step = NA ),
                           numericInput(inputId='avg_allele_freq',
                                        label='Choose the Average Allele Frequence (P bar)',
                                        value = "0", 
                                        width = NULL,
                                        min = NA,
                                        max = NA, 
                                        step = NA ),
                           numericInput(inputId='initial_pop',
                                        label='Choose the Initial Population (pX[1])',
                                        value = "0", 
                                        width = NULL,
                                        min = NA,
                                        max = NA, 
                                        step = NA ),
                           plotOutput('out_migrations_tim'),
                           tableOutput('out_migrations_tim_1')
                           
                  ),
                  tabPanel(title = 'Island Mainland Model',
                           numericInput(inputId='migration_rate1',
                                        label='Enter the migration rate',
                                        value = ".01", 
                                        width = NULL,
                                        min = NA,
                                        max = NA, 
                                        step = NA ),
                           numericInput(inputId='number_of_generations1',
                                        label='Choose the Number of Generations',
                                        value = "0", 
                                        width = NULL,
                                        min = NA,
                                        max = NA, 
                                        step = NA ),
                           plotOutput('out_migrations_imm'),
                           tableOutput('out_migrations_imm_1')
                  
                ),
                tabPanel(title = 'The General Model',
                         numericInput(inputId='number_of_generations2',
                                      label='Choose the Number of Generations',
                                      value = "0", 
                                      width = NULL,
                                      min = NA,
                                      max = NA, 
                                      step = NA ),
                         numericInput(inputId='initial_pop1',
                                      label='Choose the Initial Population (pX[1])',
                                      value = "0.1", 
                                      width = NULL,
                                      min = NA,
                                      max = NA, 
                                      step = NA ),
                         numericInput(inputId='initial_pop_y',
                                      label='Choose the Initial Population (pY[1])',
                                      value = "0.5", 
                                      width = NULL,
                                      min = NA,
                                      max = NA, 
                                      step = NA),
                         numericInput(inputId='initial_pop_z',
                                      label='Choose the Initial Population (pZ[1])',
                                      value = "0.9", 
                                      width = NULL,
                                      min = NA,
                                      max = NA, 
                                      step = NA),
                         plotOutput('out_migrations_tgm'),
                         tableOutput('out_migrations_tgm_1')
                         )
              )
            )
    )
    
    
  )
  
)
)  

server<-function(input,output){
  
  #MARKOV MODEL OF DNA MUTATION---------------------------------------------------------------------------------------------------------
  output$bargraph <- renderPlot({
    nucleotides <- c('A','T','G','C')
    prob_A <- input$propbability_a
    prob_T <- input$propbability_t
    prob_G <- input$propbability_g
    prob_C <- input$propbability_c
    order00 <- c(prob_A,prob_T,prob_G,prob_C)
    print(order00)
    names(order00)<-nucleotides
    order00seq<- sample(nucleotides,input$base,rep=T,prob=order00)
    require('seqinr')
    order00Freq <- count(order00seq,1,alphabet=nucleotides,freq=TRUE)
    order00FreqDiNt <- count(order00seq,2,alphabet=nucleotides,freq=TRUE)
    
    df_0<- data.frame(order00Freq,order00FreqDiNt)
    output$bargraph_1 <- renderTable(order00seq)
    output$bargraph_2 <- renderTable(df_0)
    
    #storing probability distribution in matrix
    
    #create a matrix
    order01_mat<- matrix(NA,nr=4,nc=4)
    
    #put column and row names 
    
    colnames(order01_mat) <- rownames(order01_mat) <- nucleotides
    #probability distribution per base
    #the rows are given according to the placement of Nt 
    order01_mat[1,] <- c(0.6,0.1,0.1,0.2)   #for A
    order01_mat[2,] <- c(0.4,0.05,0.05,0.5) #for T
    order01_mat[3,] <- c(0.05,0.2,0.7,0.05) #for G
    order01_mat[4,] <- c(0.1,0.5,0.3,0.1)   #for C
    
    #checkpoint03
    #print(order01_mat)
    
    #initial prob distribution to see porb of base to start the sequence 
    prob_initial <- order00
    
    #append the nucleotides
    names(prob_initial) <- nucleotides
    
    #checkpoint04
    #print(prob_initial)
    order01seq_generator <- function(seqlen, nucleotides, prob_initial, order01_mat){
      # seqlen = length of sequence
      # nucleotides = A,T,G,c  
      # prob_initial = initial probability distribution
      # order01seq = first order Markov chain probability distribution 
      
      #vactor to store values
      output_seq <- rep(NA,seqlen)
      #assign first base 
      output_seq[1] <- sample(nucleotides,1,prob= prob_initial)
      #creating the rest of the sequence  
      for (i in 2:length(output_seq)) {
        previousNt <- output_seq[i-1]
        currentProb <- order01_mat[previousNt,]
        output_seq[i] <- sample(nucleotides,1,prob = currentProb)
      }
      cat('|| order01seq_generator | DNA Sequence Computed ||')
      return(output_seq)
      
    }
    
    
    order01seq <- order01seq_generator(input$base,nucleotides,prob_initial,order01_mat)
    
    #checkpoint05
    #print(order01seq)
    
    #counting and calculating frequency\
    order01Freq <- count(order01seq, 1, alphabet = nucleotides, freq = TRUE)
    
    #frequency of dinucleotide
    order01FreqDiNt <- count(order01seq, 2, alphabet = nucleotides, freq = TRUE)
    
    
    layout (matrix(1:4,nr=2,nc=2))
    barplot(order00Freq,
            col=c('#145DA0','#0C2D48','#2E8BC0','#B1D4E0'),
            main='Compositional bais of each diNt Zero Order Markov Chain', 
            xlab = 'Nt', 
            ylab= 'Nt proportion')
    barplot(order00FreqDiNt,
            col=c('#050A30','#145DA0','#0C2D48','#2E8BC0','#B1D4E0'),
            main='Bias in Each Dinucleotides',
            xlab='Nt',
            ylab='Nt proportion')
    barplot(order01Freq, 
            col=c('#145DA0','#0C2D48','#2E8BC0','#B1D4E0'),
            main='Compositional bais of each diNt First Order Markov Chain', 
            xlab = 'Nt', 
            ylab= 'Nt proportion')
    barplot(order01FreqDiNt,col=c('#050A30','#145DA0','#0C2D48','#2E8BC0','#B1D4E0'),
            main='Bias in Each Dinucleotides', 
            xlab = 'Nt', 
            ylab= 'Nt proportion')
    
    df_1 <- data.frame(order01Freq,order01FreqDiNt)
    # print(df_1)
    
    output$bargraph_3 <- renderTable(order01seq)
    output$bargraph_4 <- renderTable(df_1)
    # output$bargraph_3 <- renderTable(order01FreqDiNt)
    
    
  })
  
  #DOWNLOAD BUTTON MARKOV MODEL---------------------------------------------------------------------------------------------
  
  # output$downloadData <- downloadHandler(
  #   filename=function() {
  #     paste('bargraph_4','csv',sep = '.')
  #   },
  #   content = function(file) {
  #     write.csv( output$bargraph_4(),file)
  #   }
  # )
  # 
  # output$downloadPlot <- downloadHandler(
  #   filename = function(){
  #     paste('bargraph','csv',sep = '.')
  #   },
  #   content = function(file){
  #     png(file)
  #     
  #     dev.off()
  #   }
  # )
  
  #HARDY WEINBERG MODLE ---------------------------------------------------------------------------------------------------
  output$hw_out_plot<- renderPlot({
    
    p <- seq(input$initial_value,input$final_value,input$difference)
    
    #range for q
    q <- 1-p
    
    #give alleles values 
    a1a1 <- p^2
    a1a2 <- 2 * (p * q)
    a2a2 <- q^2
    # print(c(a1a1, a1a2, a2a2))
    
    #create data frame
    hw_data_1<- data.frame(p,q,a1a1, a1a2, a2a2)
    hw_data <- data.frame(a1a1, a1a2, a2a2)
    
    hw_plot <- matplot(hw_data,type = 'b',pch=1)
    
    output$hw_out_plot_1<- renderTable(hw_data_1)
  })
  
  #MIGRATION MODLE (1) ISLAND MODEL---------------------------------------------------------------------------------------------------------
  output$out_migrations_tim <- renderPlot({
    m <- input$migration_rate
    pX <- rep(NA,T)
    pbar <- input$avg_allele_freq
    pX[1] <- input$initial_pop
    T <- input$number_of_generations
    
    for( t in 2:T)
      pX[t] <- pbar + (pX[1]-pbar)*(1-m)^t
    df <- data.frame( Generation = 1:T, Frequency = pX)
    # print(df)
    plot(y=df$Frequency,
         x=df$Generation,
         col=c('#145DA0','#0C2D48','#2E8BC0','#B1D4E0'),
         main='Island Model', 
         xlab = 'Generation', 
         ylab= 'Allele Frequency')
    output$out_migrations_tim_1 <- renderTable(df)
  })


  #MIGRATION MODLE (2) ISLAND MAINLAND MODEL---------------------------------------------------------------------------------------------------------
  output$out_migrations_imm <- renderPlot({
    migration_rates <- c(input$migration_rate1)
    results <- data.frame(m=rep(migration_rates,each=input$number_of_generations1), 
                          Generation=rep(1:input$number_of_generations1,times=1),
                          p=NA)
    
    for( m in migration_rates) {
      px <- 0
      py <- 1
      results$p[ results$m==m ] <- py
      for( t in 2:input$number_of_generations1){
        p.0 <- results$p[ results$m==m & results$Generation == (t-1) ]
        p.1 <- (1-m)*p.0 + px*m
        results$p[ results$m==m & results$Generation == t ] <- p.1
      }
    }
    results$m <- factor(results$m)
    print(results)
    
    plot(y=results$p, 
         x=results$Generation,
         col=c('#145DA0','#0C2D48','#2E8BC0','#B1D4E0'),
         main='Island Mainland Model', 
         xlab = 'Generation', 
         ylab= 'Island Frequency')
    output$out_migrations_imm_1<- renderTable(results)
  })
  
  #MIGRATION MODLE (3) THE GENERAL MODEL---------------------------------------------------------------------------------------------------------
  output$out_migrations_tgm <- renderPlot({
    T <- input$number_of_generations2
    pX <- rep(NA,T)
    pY <- rep(NA,T)
    pZ <- rep(NA,T)
    pX[1] <- input$initial_pop1
    pY[1] <- input$initial_pop_y
    pZ[1] <- input$initial_pop_z
    mXY <- 0.04
    mXZ <- 0.02
    mYZ <- 0.08
    for( gen in 2:T){
      pX[gen] <- mXY*pY[gen-1] + mXZ*pZ[gen-1] + ( 1 - (mXY+mXZ))*pX[gen-1]
      pY[gen] <- mXY*pX[gen-1] + mYZ*pZ[gen-1] + ( 1 - (mXY+mYZ))*pY[gen-1]
      pZ[gen] <- mXZ*pX[gen-1] + mYZ*pY[gen-1] + ( 1 - (mYZ+mXZ))*pZ[gen-1]
      
    }
    df <- data.frame( Generation=rep(1:T,times=3))
    df$Frequency=c(pX,pY,pZ)
    df$Population <- rep(c("X","Y","Z"), each=T)
    plot(x=df$Generation, y=df$Frequency,col=c('#145DA0','#0C2D48','#2E8BC0','#B1D4E0'),
         main='The General Model', 
         xlab = 'Generation', 
         ylab= 'Population Allele Frequency')
    output$out_migrations_tgm_1<- renderTable(df)
    
  })

}

shinyApp(ui, server)
  