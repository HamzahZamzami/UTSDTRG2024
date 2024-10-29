### 1-D Response, Linear Soil, Rigid Rock, Damped Soil
library(egg)
library(ggplot2)
library(ggpubr)

#INPUT SOIL DATA
H <- 360       # soil thickness in m (Total covering til surface)
Vs <- 500     # shear wave velocity in m/s
damping_ratio <- 2 / 100 # damping ratio in %

# Load the data
data <- read.csv("Input Motion UTS 2024.csv")  # for the code, set the column name into "time" and "acceleration" first
time <- data$time
acceleration <- data$acceleration

#Plot ground motion data (bedrock)
plot(time, acceleration, type = "l", 
     xlab = "Time (s)", ylab = "Acceleration (g)", 
     main = "Ground Motion Data (Bedrock)")

## Fourier Transform of Acceleration Data
# Perform Fourier Transform
n <- length(acceleration)
dt <- time[2] - time[1]  # Sampling interval
frequencies <- seq(0, 1 / (2 * dt), length.out = n / 2 + 1)

# Compute FFT and extract positive frequencies only
fft_result <- fft(acceleration)
amplitude <- Mod(fft_result[1:(n / 2 + 1)]) / n  # Normalize amplitude for FFT scaling (normalize complex number)

# Plot the frequency response (bedrock)
plot(frequencies, amplitude, type = "l", 
     xlab = "Frequency (Hz)", ylab = "Amplitude", 
     main = "Fourier Spectrum of Ground Motion (Bedrock)")

## Apply a Transfer Function
# F2(w) for 1D response, Linear soil, Rigid Rock, Damped soil
transfer_function <- 1 / sqrt(cos(2*pi*frequencies*H / Vs)^2 + 
                                abs(damping_ratio*(2*pi*frequencies*H / Vs))^2)
# Plot transfer_function
plot(frequencies, transfer_function, type = "l", xlab = "Frequency (Hz)", ylab = "Transfer Function", main = "Transfer Function Plot")

## SURFACE AMPLIFICATION
# Amplify Fourier-transformed data
surface_amplitude <- amplitude * transfer_function

# Plot amplified response (surface)
plot(frequencies, surface_amplitude, type = "l", xlab = "Frequency (Hz)", ylab = "Amplified Amplitude", main = "Surface Ground Motion Spectrum")

## INVERSE Fourier Transform to Get Surface Motion in Time Domain
# Combine amplified result to form a symmetric complex sequence
fft_amplified <- fft_result
fft_amplified[1:(n / 2 + 1)] <- fft_amplified[1:(n / 2 + 1)] * transfer_function
fft_amplified[(n / 2 + 2):n] <- Conj(fft_amplified[(n / 2):1])

# Inverse FFT to get time-domain response
surface_motion <- Re(fft(fft_amplified, inverse = TRUE)) / n

# Plot the original and surface motion in time domain
plot(time, surface_motion, type = "l", col = "blue", xlab = "Time (s)", ylab = "Acceleration", main = "Original vs. Surface Motion")
lines(time, acceleration, col = "red")
legend("topright", legend = c("Surface", "Original/Bedrock"), col = c("blue", "red"), lty = 1)


### PLOT SUMMARY
# 0. Plot the ground motion data
p0 <- ggplot(data.frame(time, acceleration), aes(x = time, y = acceleration)) +
  geom_line() +
  theme_bw()+
  labs(x = "Time (s)", y = "Acceleration (g)", title = "Ground Motion Data (Bedrock)")

# 1. Plot the frequency response (bedrock)
p1 <- ggplot(data.frame(frequencies, amplitude), aes(x = frequencies, y = amplitude)) +
  geom_line() +
  theme_bw() +
  labs(x = "Frequency (Hz)", y = "Amplitude", title = "Fourier Spectrum of Ground Motion (Bedrock)")

# 2. Plot transfer function
p2 <- ggplot(data.frame(frequencies, transfer_function), aes(x = frequencies, y = transfer_function)) +
  geom_line() +
  theme_bw() +
  labs(x = "Frequency (Hz)", y = "Transfer Function", title = "Transfer Function Plot")

# 3. Plot amplified response (surface)
p3 <- ggplot(data.frame(frequencies, surface_amplitude), aes(x = frequencies, y = surface_amplitude)) +
  geom_line() +
  theme_bw() +
  labs(x = "Frequency (Hz)", y = "Amplified Amplitude", title = "Surface Ground Motion Spectrum")

# 4. Plot the original and surface motion in time domain with an internal legend
p4 <- ggplot(data.frame(time, surface_motion, acceleration), aes(x = time)) +
  geom_line(aes(y = surface_motion, color = "Surface")) +
  geom_line(aes(y = acceleration, color = "Original/Bedrock")) +
  scale_color_manual(name = "Legend", values = c("Surface" = "blue", "Original/Bedrock" = "red")) +
  theme_bw() +
  labs(x = "Time (s)", y = "Acceleration", title = "Original vs. Surface Motion")

# Arrange all plots with ggarrange
ggarrange(p0, p1, p2, p3, p4, nrow = 5)
#save plots into pdf / jpg file
ggsave("All plot.pdf", # save in pdf format
       plot = last_plot(),
       width = 15, 
       height = 25,
       units = "in")
ggsave("p1.jpg", # save in pdf format
       plot = p1,
       width = 15, 
       height = 5,
       units = "in")
ggsave("p2.jpg", # save in pdf format
       plot = p2,
       width = 15, 
       height = 5,
       units = "in")
ggsave("p3.jpg", # save in pdf format
       plot = p3,
       width = 15, 
       height = 5,
       units = "in")
ggsave("p4.jpg", # save in pdf format
       plot = p4,
       width = 15, 
       height = 5,
       units = "in")
ggsave("p0.jpg", # save in pdf format
       plot = p0,
       width = 15, 
       height = 5,
       units = "in")
