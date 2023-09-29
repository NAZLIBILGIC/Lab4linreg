linreg <- setRefClass(
  "linreg",
  fields = list(
    formula = "formula",
    data = "data.frame",
    coefficients = "matrix",
    residuals = "array",
    predicted = "matrix",
    degrees_of_freedom = "numeric",
    residual_variance = "numeric",
    variance_of_coefficients = "matrix",
    t_values = "matrix",
    p_values = "matrix"
  ),
  methods = list(
    initialize = function(formula, data) {
      .self$formula <- formula
      .self$data <- data
      .self$fit_model()
    },
    fit_model <- function() {
      X <- model.matrix(.self$formula, data = .self$data)
      y <- .self$data[[all.vars(.self$formula)[1]]]
      n <- nrow(X)
      p <- ncol(X)
      
      beta_hat <- solve(t(X) %*% X) %*% t(X) %*% y
      y_hat <- X %*% beta_hat
      residuals <- y - y_hat
      df <- n - p
      sigma_sq_hat <- sum(residuals^2) / df
      Var_beta_hat <- sigma_sq_hat * solve(t(X) %*% X)
      
      # Calculate t-values and p-values
      t_values <- beta_hat / sqrt(diag(Var_beta_hat))
      p_values <- 2 * (1 - pt(abs(t_values), df))
      
      # Store the computed statistics
      .self$coefficients <- beta_hat
      .self$residuals <- residuals
      .self$predicted <- y_hat
      .self$degrees_of_freedom <- df
      .self$residual_variance <- sigma_sq_hat
      .self$variance_of_coefficients <- Var_beta_hat
      .self$t_values <- t_values
      .self$p_values <- p_values
    },
    resid <- function() {
      residuals_vector <- .self$residuals
      return(residuals_vector)
    },
    pred <- function() {
      return(.self$predicted)
    },
    coef <- function() {
      coef_vector <- as.vector(.self$coefficients)
      coef_names <- colnames(.self$coefficients)
      return(setNames(coef_vector, coef_names))
    },
    summary <- function() {
      cat("Regression Summary:\n")
      cat("Residuals:\n")
      print(head(.self$residuals))
      cat("Degrees of Freedom: ", .self$degrees_of_freedom, "\n")
      cat("Residual Variance: ", .self$residual_variance, "\n")
      cat("Variance of Coefficients:\n")
      print(.self$variance_of_coefficients)
      cat("T-values:\n")
      print(.self$t_values)
      cat("P-values:\n")
      print(.self$p_values)
    },
    show <- function(){
      cat("Coefficients:\n")
      print(.self$coef())
    },
    plot <- function() {
      library(ggplot2)
      
      df <- data.frame(Residuals = .self$resid(), Fitted = .self$pred())
      
      p1 <- ggplot(df, aes(x = Fitted, y = Residuals)) +
        geom_point() +
        stat_summary(
          fun = mean,
          fun.args = list(trim = 0.25),
          colour = "red",
          geom = "line"
        ) +
        labs(title = "Residuals vs Fitted", x = "Fitted values", y = "Residuals") +
        theme(axis.title.y = element_text(vjust = 0.5, size = 10)
        )
      
      std_residuals <- .self$residuals / sd(.self$residuals)
      sqrt_std_residuals <- sqrt(abs(std_residuals))
      
      df_std <- data.frame(Fitted = .self$pred(), Sqrt_Std_Residuals = sqrt_std_residuals)
      p2 <- ggplot(df_std, aes(x = Fitted, y = Sqrt_Std_Residuals)) +
        geom_point() +
        stat_summary(
          fun = mean,
          fun.args = list(trim = 0.25),
          colour = "red",
          geom = "line"
        ) +
        labs(title = "Scale-Location", x = "Fitted values", y = expression(sqrt("|Standardized residuals|"))) +
        theme(axis.title.y = element_text(vjust = 0.5, size = 10)) +
        theme(axis.title.x = element_text(vjust = 0.5, size = 10)) +
        theme(plot.title = element_text(size = 10, hjust = 0.5)) +
        theme_bw()
      
      print(p1)
      print(p2)
    }
  )
)


data(iris)

mod_object <- linreg(Petal.Length ~ Species, data = iris)

print(mod_object)
summary(mod_object)

#mod_object$fit_model()

#mod_object$summary()

#mod_object$plot()
