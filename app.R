library(shiny)
library(tidyverse)
library(shinyjs)
library(shinyWidgets)
library(brms)
library(tidybayes)
library(cowplot)

load("model_data.RData")

ui <- bootstrapPage(
    useShinyjs(),
    withMathJax(),
    shinyUI(
        navbarPage("Estimation of Patient-Specific Efficacy from RCTs", 
                   id = "tabs",
                   tabPanel(title = "Home",
                            fluidPage(
                                tags$style(HTML(".irs-bar {width: 100%; height: 5px; background: royalblue; border-top: 1px solid royalblue; border-bottom: 1px solid royalblue;}")),
                                tags$style(HTML(".irs-bar-edge {background: royalblue; border: 1px solid royalblue; height: 5px; border-radius: 15px 15px 15px 15px;}")),
                                tags$style(HTML(".irs-line {border: 1px solid black; height: 5px;}")),
                                tags$style(HTML(".irs-grid-text {font-family: 'arial'; color: royalblue}")),
                                tags$style(HTML(".irs-max {font-family: 'arial'; color: black;}")),
                                tags$style(HTML(".irs-min {font-family: 'arial'; color: black;}")),
                                tags$style(HTML(".irs-single {color:white; background:royalblue;}")),
                                fluidRow(
                                    column(width = 12,
                                           h4("About this Application:"),
                                           uiOutput("link_main"),
                                           hr()
                                    )),
                                fluidRow(
                                    sidebarPanel(width = 4,
                                        prettyRadioButtons("killip_in",
                                                           "Killip Class",
                                                           choices = list("I" = "I", 
                                                                          "II" = "II",
                                                                          "III" = "III",
                                                                          "IV" = "IV"),
                                                           selected = "I",
                                                           inline = TRUE,
                                                           status = "primary"),
                                        prettyRadioButtons("pmi_in", 
                                                           label = "Previous MI",
                                                           choices = list("Yes" = "yes", "No" = "no"), 
                                                           inline = TRUE,
                                                           selected = "no",
                                                           status = "primary"),
                                        prettyRadioButtons("miloc_in",
                                                           "MI Location",
                                                           choices = list("Anterior" = "Anterior", 
                                                                          "Inferior" = "Inferior",
                                                                          "Other" = "Other"),
                                                           selected = "Anterior",
                                                           inline = TRUE,
                                                           status = "primary"),
                                        sliderInput(inputId = "age_in",
                                                    label = "Age (years)",
                                                    min = 18,
                                                    max = 110,
                                                    value = 50,
                                                    step = 1,
                                                    ticks = FALSE),
                                        sliderInput(inputId = "pulse_in",
                                                    label = "Pulse (bpm)",
                                                    min = 20,
                                                    max = 250,
                                                    value = 100,
                                                    step = 1,
                                                    ticks = FALSE),
                                        sliderInput(inputId = "sbp_in",
                                                    label = "Systolic Blood Pressure (mmHg)",
                                                    min = 50,
                                                    max = 250,
                                                    value = 100,
                                                    step = 1,
                                                    ticks = FALSE),
                                        hr(),
                                        actionButton("submit", "Submit")
                                    ),
                                    mainPanel(
                                        plotOutput("plot")
                                    )
                                ))))))

server <- function(input, output, session) {
    
    plot_vals <- eventReactive(input$submit, {
        
        
        withProgress(message = "Generating Plots:", value = 0, {
        
        age_s = (input$age_in - mean(gusto$age)) / sd(gusto$age)
        pulse_s = (input$pulse_in - mean(gusto$pulse)) / sd(gusto$pulse)
        sbp_s = (input$sbp_in - mean(gusto$sysbp)) / sd(gusto$sysbp)
        Killip = input$killip_in
        miloc = input$miloc_in
        pmi = input$pmi_in
        
        incProgress(0.25, detail = "sampling (25%)")
        
        draw_sk <- fitted_draws(m_int, newdata = tibble(age_s = age_s,
                                                        pulse_s = pulse_s,
                                                        sbp_s = sbp_s,
                                                        Killip = Killip,
                                                        pmi = pmi,
                                                        miloc = miloc,
                                                        tx = "SK"),
                                n = 4000)
        
        incProgress(0.25, detail = "sampling (50%)")
        
        draw_tpa <- fitted_draws(m_int, newdata = tibble(age_s = age_s,
                                                         pulse_s = pulse_s,
                                                         sbp_s = sbp_s,
                                                         Killip = Killip,
                                                         pmi = pmi,
                                                         miloc = miloc,
                                                         tx = "tPA"),
                                n = 4000)
       
        incProgress(0.45, detail = "sampling (95%)")
        
        })
        
        tibble(sk_risk = draw_sk$.value, 
               tpa_risk = draw_tpa$.value,
               arr = sk_risk - tpa_risk,
               rr = tpa_risk / sk_risk,
               or = exp(qlogis(tpa_risk) - qlogis(sk_risk)))
        
    })
    
    

    # Dynamic Plot
    output$plot <- renderPlot({
        
        plot_grid(
            plot_grid(
                plot_vals() %>%
                    ggplot(aes(x = sk_risk, y = 1)) +
                    geom_halfeyeh(fill = "darkred", .width = 0.89) + 
                    labs(
                        x = "P(death | tx = SK)",
                        y = "Density",
                        subtitle = "Absolute Risk (SK)"
                    ) + 
                    theme_classic() +
                    theme(text = element_text(family = "Gill Sans MT"),
                          axis.text.y = element_blank(),
                          axis.ticks.y = element_blank()),
                
                plot_vals() %>%
                    ggplot(aes(x = tpa_risk, y = 1)) +
                    geom_halfeyeh(fill = "darkred", .width = 0.89) + 
                    labs(
                        x = "P(death | tx = tPA)",
                        y = "Density",
                        subtitle = "Absolute Risk (tPA)"
                    ) + 
                    theme_classic() +
                    theme(text = element_text(family = "Gill Sans MT"),
                          axis.text.y = element_blank(),
                          axis.ticks.y = element_blank()),
                ncol = 2, align = "hv", axis = "b"
            ),
            
            plot_grid(
                plot_vals() %>%
                    ggplot(aes(x = arr, y = 1)) +
                    geom_halfeyeh(fill = "blue", .width = 0.89) + 
                    geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
                    labs(
                        x = "Absolute Risk Reduction",
                        y = "Density",
                        subtitle = "Absolute Risk Reduction"
                    ) + 
                    theme_classic() +
                    theme(text = element_text(family = "Gill Sans MT"),
                          axis.text.y = element_blank(),
                          axis.ticks.y = element_blank()),
                
                plot_vals() %>%
                    ggplot(aes(x = rr, y = 1)) +
                    geom_halfeyeh(fill = "darkgreen", .width = 0.89) + 
                    geom_vline(xintercept = 1, linetype = "dashed", color = "black") + 
                    labs(
                        x = "Risk Ratio",
                        y = "Density",
                        subtitle = "Risk Ratio"
                    ) + 
                    theme_classic() +
                    theme(text = element_text(family = "Gill Sans MT"),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank()),
                
                plot_vals() %>%
                    ggplot(aes(x = or, y = 1)) +
                    geom_halfeyeh(fill = "purple", .width = 0.89) +
                    geom_vline(xintercept = 1, linetype = "dashed", color = "black") + 
                    labs(
                        x = "Odds Ratio",
                        y = "Density",
                        subtitle = "Odds Ratio"
                    ) + 
                    theme_classic() +
                    theme(text = element_text(family = "Gill Sans MT"),
                          axis.text.y = element_blank(),
                          axis.ticks.y = element_blank()),
                
                ncol = 3, align = "hv", axis = "b"),
            ncol = 1, align = "hv"
        )
    })
    
    # Link for paper
    url_post<- a("in this project summary.", 
                 href="https://benyandrew.netlify.com/blog/bayesian_rct/")
    url_gusto <- a("GUSTO trial",
                   href = "https://www.nejm.org/doi/pdf/10.1056/NEJM199309023291001")
    url_email <- a("benjamin.andrew@duke.edu", 
                   href="mailto:benjamin.andrew@duke.edu")
    output$link_main <- renderUI({
        tagList("This is an interactive tool to generate and display patient-specific efficacy estimates based on Bayesian models of RCT data. Briefly, we use Bayesian modeling to generate posterior probability distributions of covariate adjusted outcomes for patients under various treatments. For this example, data from the ", url_gusto, " was used. You can find a full description of the methods used, and a walkthrough of the modeling process ", url_post,
                "To begin, use the sliders and radio buttons to change patient characteristics. When you are done, hit the submit button to update the plots. Note: please be patient after submitting your changes as the sampling process can take 45-60 seconds. This is very much a work in progress and any comments or suggestions are always welcome: ", url_email) 
    })
    
    
}
# Run the application 
shinyApp(ui = ui, server = server)
