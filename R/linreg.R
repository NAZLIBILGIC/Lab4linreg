library(ggplot2)

# Define the linreg class
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
    fit_model = function() {
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
    resid = function() {
      residuals_vector <- .self$residuals
      return(residuals_vector)
    },
    pred = function() {
      return(.self$predicted)
    },
    coef = function() {
      coef_vector <- as.vector(.self$coefficients)
      coef_names <- colnames(.self$coefficients)
      return(setNames(coef_vector, coef_names))
    },
    summary = function() {
      cat("Regression Summary:\n")

      cat("Coefficients:\n")
      cat("             Estimate  Std. Error  t value  Pr(>|t|)\n")
      cat(sprintf("(Intercept) %9.4f %9.4f %7.2f %9.2e ***\n",
                  .self$coefficients[1],
                  sqrt(.self$variance_of_coefficients[1, 1]),
                  .self$t_values[1],
                  .self$p_values[1]))
      for (i in 2:length(.self$coefficients)) {
        cat(sprintf("Species%s  %9.4f %9.4f %7.2f %9.2e\n",
                    colnames(.self$coefficients)[i],
                    .self$coefficients[i],
                    sqrt(.self$variance_of_coefficients[i, i]),
                    .self$t_values[i],
                    .self$p_values[i]))
      }

      cat("Residual standard error: ")
      cat(sprintf("%.2f on %d degrees of freedom\n",
                  sqrt(.self$residual_variance), .self$degrees_of_freedom))

      cat("Residual Variance: ")
      cat(sprintf("%.7f\n", .self$residual_variance))

      cat("Variance of Coefficients:\n")
      print(.self$variance_of_coefficients)

      cat("\nT-values:\n")
      print(.self$t_values)

      cat("\nP-values:\n")
      print(.self$p_values)
    },


    show = function() {
      cat("linreg")
      cat("\nFormula: ", as.character(.self$formula))
      cat("\nData: ", deparse(substitute(.self$data)))
      cat("\nCoefficients:\n")
      print(.self$coef())
    },

    plot = function() {
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
        theme(
          axis.title.y = element_text(vjust = 0.5, size = 10),
          axis.title.x = element_text(vjust = 0.5, size = 10),
          plot.title = element_text(size = 10, hjust = 0.5),
          theme_bw()
        )

      # Print the plots
      print(p1)
      print(p2)
    }
      )
  )

  data(iris)
  mod_object <- linreg(Petal.Length~Species, data = iris)

  mod_object$summary()


