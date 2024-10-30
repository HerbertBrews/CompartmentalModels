# SIR Simulation Plot ----
SIRSim <- function(n, initial_infected = 1, initial_recovered = 0, beta = NA, gamma = NA, R0 = NA, infectious_period = NA, period_of_time = 365, time_interval = 1) {
  # Load required libraries
  if (!requireNamespace("deSolve", quietly = TRUE)) {
    install.packages("deSolve")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
  if (!requireNamespace("showtext", quietly = TRUE)) {
    install.packages("showtext")
  }
  library(deSolve)
  library(ggplot2)
  library(showtext)
  
  # Function for calculating SIR parameters
  calculate_SIR_parameters <- function(
    n, initial_infected, initial_recovered, beta, gamma, R0, infectious_period, period_of_time, time_interval
  ) {
    # Calculate initial susceptible population
    initial_S <- n - initial_infected - initial_recovered
    
    # Calculate gamma if infectious period is known
    if (is.na(gamma) && !is.na(infectious_period)) {
      gamma <- 1 / infectious_period
    }
    
    # Calculate gamma if beta and R0 are known
    if (is.na(gamma) && !is.na(beta) && !is.na(R0)) {
      gamma <- (beta * initial_S) / R0
    }
    
    # Calculate beta if R0 and gamma are known
    if (is.na(beta) && !is.na(R0) && !is.na(gamma)) {
      beta <- (R0 * gamma) / initial_S
    }
    
    # Calculate R0 if beta and gamma are known
    if (is.na(R0) && !is.na(beta) && !is.na(gamma)) {
      R0 <- (beta * initial_S) / gamma
    }
    
    # Calculate beta if R0 and infectious period are known
    if (is.na(beta) && !is.na(R0) && !is.na(infectious_period)) {
      gamma <- 1 / infectious_period
      beta <- (R0 * gamma) / initial_S
    }
    
    # Check if all required parameters are now known
    if (is.na(beta) || is.na(gamma) || is.na(R0)) {
      stop("Not enough information to calculate all parameters. Please provide at least two of beta, gamma, R0, or infectious_period.")
    }
    
    # Return the calculated parameters and initial states
    list(
      beta = beta,
      gamma = gamma,
      R0 = R0,
      initial_S = initial_S,
      initial_infected = initial_infected,
      initial_recovered = initial_recovered,
      period_of_time = period_of_time,
      time_interval = time_interval
    )
  }
  
  # Function for SIR model equations
  SIR_func <- function(time, state, parms) {
    with(as.list(c(state, parms)), {
      dS <- -b * I * S
      dI <- b * I * S - g * I
      dR <- g * I
      return(list(c(dS, dI, dR)))
    })
  }
  
  # Calculate SIR parameters
  SIR_values <- calculate_SIR_parameters(n, initial_infected, initial_recovered, beta, gamma, R0, infectious_period, period_of_time, time_interval)
  
  # Initial states
  init_states <- c(S = SIR_values$initial_S, I = SIR_values$initial_infected, R = SIR_values$initial_recovered)
  
  # Parameters
  parameters <- c(b = SIR_values$beta, g = SIR_values$gamma)
  
  # Time sequence
  tps <- seq(0, SIR_values$period_of_time, by = SIR_values$time_interval)
  
  # Solve the differential equations
  out_SIR <- ode(y = init_states, times = tps, func = SIR_func, parms = parameters)
  
  SIR_df <- as.data.frame(out_SIR)
  
  # Plotting function
  SIR_Plot <- function(SIR_df, SIR_values) {
    # Enable showtext and add Google font "Roboto"
    showtext_auto()
    font_add_google("Roboto")
    
    # Find the time at which the number of infected individuals (I) is maximized
    time_at_max_infected <- SIR_df$time[which.max(SIR_df$I)]
    
    ggplot(data = SIR_df, aes(x = time)) +
      geom_line(aes(y = S, colour = "Susceptible"), size = 1) +
      geom_line(aes(y = I, colour = "Infected"), size = 1) +
      geom_line(aes(y = R, colour = "Recovered"), size = 1) +
      # Add dashed vertical line at max infected time
      geom_vline(xintercept = time_at_max_infected, linetype = "dashed", colour = "grey30") +
      labs(
        title = "SIR Model Simulation",
        caption = sprintf(
          "R0 = %.2f, beta = %.4f, gamma = %.4f",
          SIR_values$R0, SIR_values$beta, SIR_values$gamma
        ),
        subtitle = paste("Max Infected at t =", round(time_at_max_infected, 2)),
        x = "Time",
        y = "Individuals (N)",
        colour = "Compartment"
      ) +
      scale_colour_manual(
        values = c("Susceptible" = "red", "Infected" = "blue", "Recovered" = "green")
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(colour = "black", size = 17.5, face = "bold", family = "Roboto"),
        plot.subtitle = element_text(colour = "grey30", size = 12.5, face = "italic", family = "Roboto"),
        axis.title = element_text(family = "Roboto", face = "bold"),
        axis.text = element_text(size = 12.5, family = "Roboto", colour = "black"),
        legend.position = "right",
        legend.title = element_text(face = "bold", family = "Roboto"),
        legend.text = element_text(family = "Roboto")
      )
  }
  
  # Plot the results
  PLOT <- SIR_Plot(SIR_df, SIR_values)
  print(PLOT)
  
  # Return the values
  return(SIR_values)
}

# Example use:
SIRSim(1000, beta = 0.0002, gamma = 0.04, period_of_time = 100)