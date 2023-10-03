#'A function for the multiple regression model
#'@import ggplot2
#'@import methods
#'@description
#'The function receives 2 arguments formula and data and uses linear algebra and returns of class linreg as an RC class.
#'A reference class called "linreg" is defined using "setRefClass". It is used to perform a linear regression analysis on a data set providing the summary of the results.The class and the methods are designed to make the linear regression model more user friendly and organised.
#'
#'The fields list is created to define the fields that an object of the linreg class will have.
#'
#'The methods list is created to define the methods (functions) that can be applied to objects of the linreg class.These methods assist to perform various operations on the linear regression model.
#'
#'Initialize: Constructor method for creating a linreg object. It takes a formula and data as arguments, stores them, and then calls fit_model to perform the regression.
#'The fit model uses the provided formula and data, calculating coefficients, residuals, and various statistics.
#'The resid returns the residuals as a vector.
#'The pred returns the predicted values.
#'The coef returns the coefficients with their names.
#'The summary prints a summary of the linear regression model, including coefficients, standard errors, t-values, p-values, residual standard error, and variance of coefficients.
#'The show displays basic information about the linreg object, including the formula, data, and coefficients.
#'The plot generates and displays two diagnostic plots for the linear regression model: "Residuals vs Fitted" and "Scale-Location".
#'
#'@param formula The formula argument should take a formula object. It includes the formula given from the user so that the function will do the calculations for the regression.
#'@param data It must be a data frame and it includes the data for the regression.
#'@param coefficients It must be a matrix. A matrix to store the regression coefficients.
#'@param residuals It must be an array. It is used to store the residuals.
#'@param predicted It must be a matrix. It is used to store the predicted values.
#'@param degrees_of_freedom The type should be numeric. It is used to store the degrees of freedom.
#'@param residual_variance A numeric value for the residual variance.
#'@param variance_of_coefficients It must be a matrix. It is used to store the variance of the coefficients.
#'@param t_values It must be a matrix to store the t values.
#'@param p_values It must be a matrix to store the p values.
#'
#'@return The function returns an object of class linreg as an RC class.
#'
#'@export
# Define the linreg class
linreg <- setRefClass(
  "linreg",
  fields = list(
    formula = "formula",
    data = "data.frame",
    formula_str = "character",
    data_str = "character",
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
      .self$formula_str <- deparse(substitute(formula))
      .self$data_str <- deparse(substitute(data))
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

      #print.default(as.vector(t_values))
      p_values <- 2 * (1 - pt(abs(t_values), df))

      # Handle cases of zero standard errors
      zero_se_indices <- which(diag(Var_beta_hat) == 0)
      t_values[zero_se_indices] <- NA
      p_values[zero_se_indices] <- NA


      formatted_p_values <- format(p_values, scientific = TRUE)


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
      return(setNames(as.vector(.self$coefficients), rownames(.self$coefficients)))
    },
    summary = function() {
      cat("Regression Summary:\n")

      # Create a data frame for the coefficients
      coefficients_df <- data.frame(
        Variable = rownames(p_values),
        Estimate = as.vector(coefficients),
        Std.Error = as.vector(sqrt(diag(variance_of_coefficients))),
        t.value = as.vector(t_values),
        p.value = as.vector(p_values),
        Significance = '***',
        stringsAsFactors = FALSE
      )

      print.data.frame(coefficients_df, row.names = FALSE)

      cat("Residual standard error:", sprintf("%.2f on %d degrees of freedom\n", sqrt(.self$residual_variance), .self$degrees_of_freedom))
      # returns variable without printing anything
      invisible(coefficients_df)
    },

    print = function() {
      cat("Call:\n")
      cat(sprintf("linreg(formula = %s, data = %s)\n", .self$formula_str, .self$data_str))

      cat("Coefficients:\n")
      print.default(.self$coef())
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

# data(iris)
# mod_object <- linreg(Petal.Length~Species, data = iris)
# mod_object$print()
# mod_object$summary()


