###### The following codes are used to generate Figure 3 in Wang et al. "Towards a unified understanding: the linkage of MaxEnt, ETRHEQ, and SFE Models in estimating evapotranspiration"
###### The codes was developed by Yi Wang.


# Load necessary library
library(deSolve) 
library(ggplot2) 
library(gridExtra)
library(grid)

# Initialize lists for storing variables globally
LE_list <<- list() # latent heat (LE, W/m2)
RH_list <<- list() # relative humidity (RH)
RHs_list <<- list() # surface relative humidity (RHs)
delta_RH_list <<- list() # the difference between RHs and RH
var_RH_list <<-list() # the vertical variance of RH
Ts_list <<- list() # surface temperature (K)
qs_list <<- list() # surface specific humidity (kg/kg)
D_list <<- list() # Dissipation
T_list <<-list() # air temperature in the box (K)
q_list <<-list() # air specific humidity in the box
LE_SFE_list <<-list() # estimated LE based on the SFE model (McColl et al. 2019)

# Define parameters and constants
params <- list(ga = 1/50, gs = 1/900, Tg = 8 + 273.15, kg = 0.3, 
               Fsd = 700, Fld = 300, alpha = 0.5, h = 1000, emi = 0.98, 
               SB = 5.67037442e-8, epsilon = 0.622, cp = 1005, 
               lambda = 2.502e6, gamma <- 1005/2.502e6, g = 9.81, rho = 1.293, 
               dg = 1)

# ga is the aerodynamic conductance (m/s)
# gs is the surface conductance (m/s)
# Tg is the ground heat flux reference temperature (K)
# kg is the ground heat flux conductivity (J·m−1·s−1·K−1)
# Fsd is the downwelling shortwave radiation (W/m2)
# Fld is the downwelling longwave radiation (W/m2)
# alpha is the surface albedo 
# h is the boundary layer height (m)
# emi is the surface emissivity
# SB is the Stefan–Boltzmann constant
# epsilon is the dimensionless ratio of the gas constant for dry air to water vapor
# cp is the specific heat of air at constant pressure (J·kg−1·K−1)
# lambda is the latent heat of vaporization (J·kg−1)
# gamma = cp/lambda
# g is the gravitational acceleration (m·s−2)
# rho is the air density (kg/m3)
# dg is the the depth to a deep soil layer with constant reference temperature (m)


# Define a function to solve for Ts
solve_for_Ts <- function(T, q, params) {
  fn <- function(Ts) {
    Rn <- (1 - params$alpha) * params$Fsd + params$Fld - params$emi * params$SB * Ts^4 #under Eq. 8 in McColl et al. 2019
    G <- params$kg * (Ts - params$Tg) / params$dg #under Eq. 8 in McColl et al. 2019
    H <- params$rho * params$cp * params$ga * (Ts - T) #Eq. 9 in McColl et al. 2019
    
    # Calculate q*_{Ts}
    estar_Ts <- 611.2 * exp(17.67 * (Ts - 273.15) / (Ts - 29.65)) #Eq. 11 in Salvucci and Gentine 2013
    qstar_Ts <- params$epsilon * estar_Ts / (101325 - (1 - params$epsilon) * estar_Ts) #Eq.9 in Salvucci and Gentine 2013
    LE <- params$rho * params$lambda * params$ga * params$gs / (params$ga + params$gs) * (qstar_Ts - q) #Eq. 10 in McColl et al. 2019
    
    Rn - G - H - LE 
  }
  
  # Use a numerical solver to find Ts
  Ts_solution <- uniroot(fn, lower = 273, upper = 350)$root  
  return(Ts_solution)
}

# Define a function to solve for Ts
solve_for_qs <- function(Ts, q, params) {
  estar_Ts <- 611.2 * exp(17.67 * (Ts - 273.15) / (Ts - 29.65)) #Eq. 11 in Salvucci and Gentine 2013
  qstar_Ts <- params$epsilon * estar_Ts / (101325 - (1 - params$epsilon) * estar_Ts) #Eq.9 in Salvucci and Gentine 2013
  
  # Calculate LE using the first LE equation
  LE <- params$rho * params$lambda * params$ga * params$gs / (params$ga + params$gs) * (qstar_Ts - q) #Eq. 10 in McColl et al. 2019
  
  
  qs <- qstar_Ts - LE / (params$rho * params$lambda * params$gs) #Eq. 10 in McColl et al. 2019
  
  return(qs)
}


# The temperature and moisture budgets for the box model
model <- function(t, state, params) {
  with(as.list(c(state, params)), {
    # Solve for Ts and qs
    Ts <- solve_for_Ts(T, q, params) 
    qs <- solve_for_qs(Ts, q, params) 
    
    # Access the global lists
    globalLE_list <- get("LE_list", envir = .GlobalEnv)
    globalRH_list <- get("RH_list", envir = .GlobalEnv)
    globalRHs_list <- get("RHs_list", envir = .GlobalEnv)
    globaldelta_RH_list <- get("delta_RH_list", envir = .GlobalEnv)
    globalvar_RH_list <- get("var_RH_list", envir = .GlobalEnv)
    globalTs_list <- get("Ts_list", envir = .GlobalEnv)
    globalqs_list <- get("qs_list", envir = .GlobalEnv)
    globalD_list <- get("D_list", envir = .GlobalEnv)
    globalT_list <- get("T_list", envir = .GlobalEnv)
    globalq_list <- get("q_list", envir = .GlobalEnv)
    globalLE_SFE_list <- get("LE_SFE_list", envir = .GlobalEnv)
    
    
    # Define additional parameters
    gr <- h / (10 * 86400)  # Relaxation conductance (h divided by 10 days in seconds)
    Tr <- 10 + 273.15       # Atmospheric cooling relaxation temperature in Kelvin
    qr <- 0.0011574074074074073    # Atmospheric drying relaxation specific humidity in kg/kg
    
    # Calculate Rn (Net Radiation) Ref: McColl et al. 2019
    Rn <- (1 - alpha) * Fsd + Fld - emi * SB * Ts^4 
    
    # Calculate G (Ground Heat Flux) Ref: McColl et al. 2019
    G <- kg * (Ts - Tg) / dg
    
    # Primary calculation of H (Sensible Heat Flux) Eq. 9 in McColl et al. 2019
    H <- rho * cp * ga * (Ts - T)
    
    
    # calculating surface and atmospheric saturated specific humidity
    estar_T <- 611.2 * exp(17.67 * (T - 273.15) / (T - 29.65))
    qstar_T <- params$epsilon * estar_T / (101325 - (1 - params$epsilon) * estar_T)
    estar_Ts <- 611.2 * exp(17.67 * (Ts - 273.15) / (Ts - 29.65))
    qstar_Ts <- params$epsilon * estar_Ts / (101325 - (1 - params$epsilon) * estar_Ts)
    
    
    # Calculate air RH (Relative Humidity)
    RH <- q / qstar_T
    
    # Check for condensation
    if (RH >= 1) {
      # Calculate condensation rate
      condensation_rate <- max(0, q - qstar_T)
      
      # Adjust LE for condensation
      LE <- LE - condensation_rate * params$lambda
      
      # Update dTdt and dqdt for condensation
      dTdt <- dTdt + condensation_rate * params$lambda / (rho * cp)
      dqdt <- dqdt - condensation_rate / (rho * h)
    }
    
    
    # Calculate surface RH
    RHs <- qs / qstar_Ts
    
    # Calculate the vertical difference in RH
    delta_RH <- RH - RHs
    
    # Calculate the vertical variance of RH
    var_RH <- 1/4*(RHs - RH)^2
    
    
    # Calculate the slope of the relation between saturation vapor pressure and temperature
    delta <- (qstar_Ts - qstar_T) / (Ts - T)
    
    
    # Calculate LE (Latent Heat Flux) from the energy balance equation
    LE <- Rn - G - H
    
    # Calculate LE based on the SFE model (Eq. 5 in McColl et al. 2019) 
    LE_SFE <-(RH*delta*lambda/cp/((RH*delta*lambda/cp)+1))*(Rn-G)
    
    # Thermal inertia for G
    Is <- 800
    
    # Thermal inertia for H
    Ia <- rho * cp * sqrt(ga)
    
    # Thermal inertia for LE
    Ie <- delta/gamma*RHs*Ia
    
    # Calculte the dissipation
    D <- 2 * ((G)^2) / Is + 2 * (H^2) / Ia + 2 * (LE^2) / Ie
    
    
    # Store LE and RH in global lists
    globalLE_list[[as.character(t)]] <- LE
    globalRH_list[[as.character(t)]] <- RH
    globalRHs_list[[as.character(t)]] <- RHs
    globaldelta_RH_list[[as.character(t)]] <- delta_RH
    globalvar_RH_list[[as.character(t)]] <- var_RH
    globalTs_list[[as.character(t)]] <- Ts
    globalqs_list[[as.character(t)]] <- qs
    globalD_list[[as.character(t)]] <- D
    globalT_list[[as.character(t)]] <- T
    globalq_list[[as.character(t)]] <- q
    globalLE_SFE_list[[as.character(t)]] <- LE_SFE
    
    assign("LE_list", globalLE_list, envir = .GlobalEnv)
    assign("RH_list", globalRH_list, envir = .GlobalEnv)
    assign("RHs_list", globalRHs_list, envir = .GlobalEnv)
    assign("delta_RH_list", globaldelta_RH_list, envir = .GlobalEnv)
    assign("var_RH_list", globalvar_RH_list, envir = .GlobalEnv)
    assign("Ts_list", globalTs_list, envir = .GlobalEnv)
    assign("qs_list", globalqs_list, envir = .GlobalEnv)
    assign("D_list", globalD_list, envir = .GlobalEnv)
    assign("T_list", globalT_list, envir = .GlobalEnv)
    assign("q_list", globalq_list, envir = .GlobalEnv)
    assign("LE_SFE_list", globalLE_SFE_list, envir = .GlobalEnv)
    
    
    # Force RH to equal RHs if RH > RHs (the steady state)
    if (RH > RHs) {
      RH <- RHs
      
      # Set the model to steady state
      dTdt <- 0
      dqdt <- 0
      
    }else{
      # Differential equations for T and q following McColl et al. (2019)
      dTdt <- H / (rho * cp * h) + gr / h * (Tr - T)
      dqdt <- LE / (rho * lambda * h) + gr / h * (qr - q)
      
      # Differential equations for T and q ignoring the atmospheric cooling and drying terms (optional)
      #dTdt <- H / (rho * cp * h)
      #dqdt <- LE / (rho * lambda * h)
    }
    
    
    return(list(c(dTdt, dqdt)))
    
  })
  
  
}


# Initial conditions, now including Ts and qs
initial_state <- c(T = 300, q = 0.001) # Initial Ts and qs need to be set

# Time vector (25 days)
time <- seq(0, 25 * 86400, by = 3600) 

# Initialize global variables for LE and RH
LE <- numeric(length = length(time))
RH <- numeric(length = length(time))
RHs <- numeric(length = length(time))
delta_RH <- numeric(length = length(time))
var_RH <- numeric(length = length(time))
Ts <- numeric(length = length(time))
qs <- numeric(length = length(time))
D <- numeric(length = length(time))
T <- numeric(length = length(time))
q <- numeric(length = length(time))
LE_SFE <- numeric(length = length(time))

names(LE) <- time
names(RH) <- time
names(RHs) <- time
names(delta_RH) <- time
names(var_RH) <- time
names(Ts) <- time
names(qs) <- time
names(D) <- time
names(T) <- time
names(q) <- time
names(LE_SFE) <- time



# Solve the model which will provide T and q at each timestep. Output only stores T and q while the lists stored every variable of interest.
output <- ode(y = initial_state, times = time, func = model, parms = params)

# Processing the output (optional, depending on your analysis needs)
output_df <- as.data.frame(output)


# Combine the lists into a data frame
results_df <- data.frame(
  Time = as.numeric(names(LE_list)),  # assuming the names of the lists are the time values
  LE = unlist(LE_list),
  RH = unlist(RH_list),
  RHs = unlist(RHs_list),
  delta_RH = unlist(delta_RH_list),
  var_RH = unlist(var_RH_list),
  Ts = unlist(Ts_list),
  qs = unlist(qs_list),
  D = unlist(D_list),
  T = unlist(T_list),
  q = unlist(q_list),
  LE_SFE = unlist(LE_SFE_list)
)

results_df$time <- results_df[, "Time"] / 86400  # Convert time from seconds to days


# Plot the LE and the SFE estimated LE (LE-SFE)

plot_LE <- ggplot(data = results_df, aes(x = time)) +
  geom_ribbon(aes(ymin = LE - 10, ymax = LE + 10), fill = "grey", alpha = 0.5) +
  geom_line(aes(y = LE, color = "LE")) +
  geom_line(aes(y = LE_SFE, color = "SFE-LE")) +
  labs(title = 'a', x = 'Days', y = 'LE ('*W~m^-2*')') +
  scale_color_manual(name = "Legend", values = c("LE" = "blue", "SFE-LE" = "red")) +
  theme_bw()+theme(axis.title.y = element_text(face="plain", size=11))+
  theme(
    legend.background = element_rect(linetype="solid",color ="1"),
    legend.text = element_text(size = 7),
    legend.title = element_blank(),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid=element_blank(),
    panel.border = element_rect(colour = "black", size=1.6),
    axis.ticks = element_line(size = 1.6),
    legend.position = c(0.8,0.2),
    legend.box = "vertical")

plot_LE

# Plot T and Ts over Time
plot_T_Ts <- ggplot(data = results_df, aes(x = time)) +
  geom_line(aes(y = Ts, linetype = "Ts"), color = 'blue') +
  geom_line(aes(y = T, linetype = "T"), color = 'blue') +
  labs(title = 'b', x = 'Days', y = 'Temperature (K)') +
  scale_linetype_manual(name = "Line Type", values = c( "Ts" = "dashed","T" = "solid")) +
  theme_bw()+theme(axis.title.y = element_text(face="plain", size=11))+
  theme(
    legend.background = element_rect(linetype="solid",color ="1"),
    legend.text = element_text(size = 7),
    legend.title = element_blank(),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid=element_blank(),
    panel.border = element_rect(colour = "black", size=1.6),
    axis.ticks = element_line(size = 1.6),
    legend.position = c(0.8,0.2),
    legend.box = "vertical")

plot_T_Ts

# Plot RH change over time
plot_RH <- ggplot(data = results_df, aes(x = time)) +
  geom_line(aes(y = RH), color = 'blue') +
  labs(title = 'c', x = 'Days', y = 'Relative humidity') +
  ylim(0,1)+
  theme_bw()+theme(axis.title.y = element_text(face="plain", size=11))+
  theme(
    legend.background = element_rect(linetype="solid",color ="1"),
    legend.text = element_text(size = 7),
    legend.title = element_blank(),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid=element_blank(),
    panel.border = element_rect(colour = "black", size=1.6),
    axis.ticks = element_line(size = 1.6),
    legend.position = "none")

plot_RH

# Plot RH and RHs over Time
plot_RH_RHs <- ggplot(data = results_df, aes(x = time)) +
  geom_line(aes(y = RHs, linetype = "RHs"), color = 'blue') +
  geom_line(aes(y = RH, linetype = "RH"), color = 'blue') +
  labs(title = 'b', x = 'Days', y = 'Relative humidity') +
  scale_linetype_manual(name = "Line Type", values = c( "RHs" = "dashed","RH" = "solid")) +
  theme_bw()+theme(axis.title.y = element_text(face="plain", size=11))+
  theme(
    legend.background = element_rect(linetype="solid",color ="1"),
    legend.text = element_text(size = 7),
    legend.title = element_blank(),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid=element_blank(),
    panel.border = element_rect(colour = "black", size=1.6),
    axis.ticks = element_line(size = 1.6),
    legend.position = c(0.8,0.2),
    legend.box = "vertical")

plot_RH_RHs

#plot vertical difference in RH and dissipation with time

scale_factor <- max(results_df$D) / 0.05

plot_delta_RH_D <- ggplot(data = results_df, aes(x = time)) +
  geom_line(aes(y = D, color = "Dissipation")) +
  geom_line(aes(y = delta_RH * scale_factor, color = "RH-RHs")) +
  scale_y_continuous(
    name = "Dissipation",
    limits = c(-0.05 * scale_factor, max(results_df$D)),  # Adjust limits for primary axis
    sec.axis = sec_axis(~ . / scale_factor, name = "RH-RHs")  # Scale secondary axis
  ) +
  scale_color_manual(name = "Legend", values = c("Dissipation" = "red", "RH-RHs" = "blue")) +
  labs(title = 'd',x = 'Days') +
  theme_bw() +
  theme(
    axis.title.y = element_text(face="plain", size=11),
    legend.background = element_rect(linetype="solid", color ="1"),
    legend.text = element_text(size = 7),
    legend.title = element_blank(),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", size=1.6),
    axis.ticks = element_line(size = 1.6),
    legend.position = c(0.8,0.2),
    legend.box = "vertical"
  )
plot_delta_RH_D


#plot vertical variance of RH and dissipation over time
plot_var_RH_D <- ggplot(data = results_df, aes(x = time)) +
  geom_line(aes(y = D, color = "Dissipation")) +
  geom_line(aes(y = (var_RH - min(results_df$var_RH)) / (max(results_df$var_RH) - min(results_df$var_RH)) * max(results_df$D), color = "Var(RH)")) +
  scale_y_continuous(
    name = "Dissipation",
    sec.axis = sec_axis(~ . / max(results_df$D) * (max(results_df$var_RH) - min(results_df$var_RH)) + min(results_df$var_RH), name = "Var(RH)")
  ) +
  scale_color_manual(name = "Legend", values = c("Dissipation" = "red", "Var(RH)" = "green")) +
  labs(title = 'e', x = 'Days') + 
  theme_bw() +
  theme(
    axis.title.y = element_text(face="plain", size=11),
    legend.background = element_rect(linetype="solid", color ="1"),
    legend.text = element_text(size = 7),
    legend.title = element_blank(),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", size=1.6),
    axis.ticks = element_line(size = 1.6),
    legend.position = c(0.8,0.8),
    legend.box = "vertical"
  )
plot_var_RH_D


#Output the graphs
model_simulation <- grid.arrange(plot_LE,plot_T_Ts,plot_RH,plot_delta_RH_D,plot_var_RH_D ,nrow=3)
ggsave("Figure 3 .png", plot=model_simulation, dpi =800, width =10 , height = 10)


