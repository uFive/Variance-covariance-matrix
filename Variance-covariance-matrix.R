# Load necessary packages
library(MASS)  # Inverse matrix calculation

# 1. Input data
data <- data.frame(
  Treatment = c("A", "B", "C"),  
  yi = c(-0.1431, -0.1087, -0.0625),  # Calculate ln(Response Ratio) (corresponding to RRA, RRB, RRC)
  vi = c(0.0014, 0.0016, 0.0024)   # Variance (corresponding to the variance of each RR)
)

# 2. Shared control group information
SD_C <- 1.224744871   # Control group standard deviation
X_C <- 16.5   # Control group mean
N_C <- 6       # Control group sample size

# 3. Calculate the variance-covariance matrix V
cov_AB <- (SD_C^2) / (N_C * X_C^2)  # Covariance between A and B (shared control)
cov_AC <- (SD_C^2) / (N_C * X_C^2)  # Covariance between A and C (shared control)
cov_BC <- (SD_C^2) / (N_C * X_C^2)  # Covariance between B and C (shared control)

V <- matrix(c(data$vi[1], cov_AC, cov_AB,
              cov_AC, data$vi[2], cov_BC,
              cov_AB, cov_BC, data$vi[3]), 
            nrow=3, byrow=TRUE)

# 4. Calculate the weighted GLS aggregated response ratio (RR_agg)
X <- matrix(1, nrow=3, ncol=1)  # X is a column vector of ones
E <- matrix(data$yi, nrow=3, ncol=1)  # E is the effect size column vector

# Calculate GLS estimated aggregated RR
V_inv <- ginv(V)  # Calculate the inverse of V (using generalized inverse for rank-deficient matrix)
RR_agg <- solve(t(X) %*% V_inv %*% X) %*% (t(X) %*% V_inv %*% E)

# Calculate aggregated variance
var_RR_agg <- solve(t(X) %*% V_inv %*% X)

# Calculate correlation (Eq. 9)
r_RRA_RRB <- cov_AB / (sqrt(data$vi[1]) * sqrt(data$vi[2]))

# 5. Output results
cat("Aggregated ln(Response Ratio):", RR_agg, "\n")
cat("Aggregated variance:", var_RR_agg, "\n")
cat("Correlation r(RRA, RRB):", r_RRA_RRB, "\n")
cat("Aggregated Response Ratio (inverse transformation):", exp(RR_agg), "\n")
