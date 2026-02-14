
library(shiny)
library(dplyr)
library(ggplot2)
library(cowplot)

Name<-c("AUT","BLR","BEL","BGR","CAN","CHL",
        "HRV","CZE","DNK","EST","FIN",
        "DEUTNP","HUN","ISL","IRL",
        "ITA","JPN","LTU","NLD","NOR","POL","PRT",
        "KOR","RUS","SVK","SVN","ESP","SWE","CHE",
        "TWN","GBR_NP","GBRTENW","GBR_NIR","GBR_SCO",
        "UKR","USA")

Name2<-c("Austria","Belarus","Belgium","Bulgaria",
         "Canada","Chile","Croatia","Czechia",
         "Denmark","Estonia","Finland",
         "Germany",
         "Hungary","Iceland","Ireland","Italy",
         "Japan","Lithuania","The Netherlands","Norway",
         "Poland","Portugal","Republic of Korea",
         "Russia","Slovakia","Slovenia","Spain",
         "Sweden","Switzerland","Taiwan","United Kingdom",
         "England and Wales","Northern Ireland","Scotland",
         "Ukraine","USA") 

country_lookup <- setNames(Name, Name2)


## data needed


D1<-read.table("TFRComponents.csv",sep=",",header=TRUE, skip=0)

D2<-read.table("TimeTrends.csv",sep=",",header=TRUE, skip=0)

D3<-read.table("AgeComponents.csv",sep=",",header=TRUE, skip=0)

D4<-read.table("AgeDecomposition.csv",sep=",",header=TRUE, skip=0)


# ---- UI ----
ui <- fluidPage(
  
  titlePanel("TFR Components Explorer"),
  
  sidebarLayout(
    sidebarPanel(
      
      selectInput("Country",
                  "Select Country:",
                  choices = country_lookup,
                  selected = country_lookup[1]),
      
     
      radioButtons("figure",
                   "Select Output:",
                   choices = c("Figure 1: Time Trends in TFR Components",
                               "Figure 2: TFR Decomposition",
                               "Figure 3: Age-Components of TFR",
                               "Figure 4: Age Decomposition"),
                   selected = "Figure 1: Time Trends in TFR Components")
    ),
    
    mainPanel(
      plotOutput("main_plot", height = "500px")
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {
  
  # Reactive filtered data
  data_D1 <- reactive({
    D1 %>% filter(Country == input$Country)
      })
  
  data_D2 <- reactive({
    D2 %>% filter(Country == input$Country)
  })
  
  data_D3 <- reactive({
    D3 %>% filter(Country == input$Country)
  })
  
  
  data_D4 <- reactive({
    D4 %>% filter(Country == input$Country)
  })
  
   output$main_plot <- renderPlot({
    
     country_full <- names(country_lookup)[country_lookup == input$Country]
    
    # -------------------------------
    # FIGURE 1: Components of TFR
    # -------------------------------
    if (input$figure == "Figure 1: Time Trends in TFR Components") {
     
     df <- data_D1()
     
     d <- range(df$Year)
     
     df_long <- tidyr::pivot_longer(
       df,
       cols = c(TFR, TFR1),
       names_to = "Measure",
       values_to = "Value"
     )
     
     
     ggplot(df_long, aes(x = PC, y = Value, color = Measure)) +
       geom_line(linewidth = 1.2) +
       geom_point(size = 2) +
       labs(
         x = "Proportion Childless",
         y = "TFR and TFR+1",
         color = "",
         title = paste0(
           "Time Trends in TFR and its Components, ", country_full,", ",d[1],"–",d[2],":  
 TFR of Mothers (TFR+1) and Proportion Childless (PC)") 
       ) +
       scale_color_manual(
         values = c("TFR" = "#1b9e77", 
                    "TFR1" = "#d95f02"),
         labels = c("TFR" = "TFR",
                    "TFR1" = "TFR+1")
       )+
       theme_minimal(base_size = 15)+ 
       theme(
         axis.text  = element_text(size = 14),
         axis.title = element_text(size = 16),
         legend.title = element_text(size = 14),
         legend.text  = element_text(size = 13),
         plot.title = element_text(
           size = 18,        # bigger
           face = "bold",    # bold
           hjust = 0.5       # center it
         ))
     
    }
    
     
     # -------------------------------
     # FIGURE 2: TFR decomposition
     # -------------------------------
     else if (input$figure == "Figure 2: TFR Decomposition") {
       
       d <- range(data_D3()$Year, na.rm = TRUE)
       
       data <- data_D2()
       
       ggplot(data, aes(x = Period, y = values, fill = factor(Eq3))) +
         geom_col() +
         geom_hline(yintercept = 0, linewidth = 0.4) +
         labs(x="Period (years)",
              y="Contribution to TFR change",fill="Components",
              title = paste0(
              "Decomposition of TFR Change, ", country_full,", ",d[1],"–",d[2],":  
 TFR of Mothers (TFR+1) and Proportion Childless (PC)") 
              )+
         scale_fill_discrete(
           breaks = c("1", "2"),
           labels = c("TFR+1", "PC")
         ) +
         theme_minimal() + 
         theme(
           axis.text  = element_text(size = 14),
           axis.title = element_text(size = 16),
           legend.title = element_text(size = 14),
           legend.text  = element_text(size = 13),
           plot.title = element_text(
             size = 18,        # bigger
             face = "bold",    # bold
             hjust = 0.5       # center it
           ))
       
       
     }
     
     
     
    # -------------------------------
    # FIGURE 3: Age-Components of TFR 
    # -------------------------------
    else if (input$figure == "Figure 3: Age-Components of TFR") {
      
      data <- data_D3()
      
      data$Year <- factor(data$Year)
      
      d <- unique(data$Year)
      
      p1 <-ggplot(data, aes(x = Age, y = AFR,
                            color = Year)) +
        geom_line(linewidth = 1) +
        scale_color_brewer(type = "qual", palette = 2) +
        labs(
          x = "Age",
          y = "AFR"
        ) +
        theme_minimal(base_size = 14) +
        guides(
          color = guide_legend(order = 1)
        )+
        theme(legend.position = "none")
      
      
      p2 <- ggplot(data, aes(x = Age, y = PC,
                             color = Year)) +
        geom_line(linewidth = 1) +
        scale_color_brewer(type = "qual", palette = 2) +
        scale_y_continuous(limits = c(0,1)) +
        labs(x = "Age",
             y = "Proportion Childless") +
        theme_minimal(base_size = 14) +
        guides(
          color = guide_legend(order = 1),
          linetype = guide_legend(order = 2)
        )+
        theme(legend.position = "none")
      
      
      p3 <-ggplot(data, aes(x = Age, y = AFR1,
                            color = Year)) +
        geom_line(linewidth = 1) +
        scale_color_brewer(type = "qual", palette = 2) +
        labs(
          x = "Age",
          y = "AFR+1"
        ) +
        theme_minimal(base_size = 14) +
        guides(
          color = guide_legend(order = 1),
          linetype = guide_legend(order = 2)
        )+
        theme(legend.position = "bottom")
      
      
      
      
      ############
      plegend <- ggplot(data, aes(x = Age, y = AFR,
                                  color = Year)) +
        geom_line(linewidth = 1) +
        scale_color_brewer(type = "qual", palette = 2) + 
        theme_minimal(base_size = 14) +
        theme(legend.position = "bottom")
      
      
      legend <- get_legend(plegend)
      
      p1 <- p1 + theme(legend.position = "none")
      p2 <- p2 + theme(legend.position = "none")
      p3 <- p3 + theme(legend.position = "none")
      
      combined <- plot_grid(
        p1, p2, p3,
        ncol = 1,
        align = "v",
        axis = "lr"
      )
      
      
      final_plot_base <- plot_grid(
        combined,
        legend,
        ncol = 1,
        rel_heights = c(1, 0.1)
      )
      
      final_plot <- ggdraw() +
        draw_label(paste0(
            "The Accumulated Fertility Rate (AFR), Proportion Childless (PC), 
and AFR for Mothers (AFR+1), ", country_full,", ",d[1], "–", d[2]),
          x = 0.5,
          y = 0.999,        # push title to very top
          hjust = 0.5,
          vjust = 1,
          fontface = "bold",
          size = 16
        ) +
        draw_plot(final_plot_base, y = 0, height = 0.93)
      
      final_plot
    }
    
      
      
      
      
      
      
     
    # -------------------------------
    # FIGURE 4: Age Decomposition
    # -------------------------------
    else if (input$figure == "Figure 4: Age Decomposition") {
      
      d <- range(data_D3()$Year, na.rm = TRUE)
      
      data <- data_D4() 
      
        ggplot(data, aes(x = age, y = value, fill = component)) +
        geom_col(width = 0.85) +
        geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.6) +
        scale_x_continuous(limits = c(12, 50)) +
        labs(x = "Age", 
             y = "Age-contribution", 
             fill = "Components",
             title = paste0(
               "Age- and Component-Decomposition of the Change in TFR, ",d[1], "–", d[2], ", 
", country_full,": Proportion Childless (PC) and Fertility of Mothers (f+1)" ))+
        theme_minimal() +
        theme(
          axis.text  = element_text(size = 14),
          axis.title = element_text(size = 16),
          legend.title = element_text(size = 14),
          legend.text  = element_text(size = 13),
          legend.position = "bottom",
          plot.title = element_text(
            size = 18,        # bigger
            face = "bold",    # bold
            hjust = 0.5       # center it
          )
        )
  }
    
})}

shinyApp(ui, server)