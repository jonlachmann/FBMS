
#' Impute Missing Values in the Data
#'
#' Imputes missing values in the data using median imputation based on the data  set.
#'
#' @param object A fitted model object with an "imputed" attribute indicating columns to impute.
#' @param x A matrix or data frame
#' @return A matrix with imputed values and additional columns for missingness indicators.
#' @export
#' @examples
#' \donttest{
#' set.seed(123)
#' x <- matrix(rnorm(60), 10, 6)
#' colnames(x) <- paste0("X", 1:6)
#' x[1:2, 1] <- NA  # Introduce missing values
#' model <- list(imputed = c(1))
#' attr(model, "imputed") <- c(1)
#' x_imputed <- impute_x(model, x)
#' dim(x_imputed)  # 10 rows, 7 columns (6 original + 1 missingness indicator)
#' any(is.na(x_imputed))  # FALSE, no missing values
#' }
impute_x <- function (object, x) {
  if (!is.null(attr(object, which = "imputed"))) {
    df <- data.frame(x)
    na.matr <- data.frame(1 * (is.na(df)))
    cm <- colMeans(na.matr)
    na.matr <- na.matr[, attr(object, which = "imputed")]
    names(na.matr) <- paste0("mis_", names(na.matr))
    for (i in which(cm != 0)){
      med <- median(df[[i]], na.rm = TRUE)
      if(is.na(med))
        stop("No data for imputation in test set, provide x_train in predict!")
      df[[i]][is.na(df[[i]])] <- med
    }
    return(as.matrix(data.frame(df,na.matr)))
  }
  return(as.matrix(x))
}

#' Impute Missing Values in Test Data Using Training Data
#'
#' Imputes missing values in the test data using median imputation based on the training set.
#'
#' @param object A fitted model object with an "imputed" attribute indicating columns to impute.
#' @param x_test A matrix or data frame containing the test data.
#' @param x_train A matrix or data frame containing the training data.
#' @return A matrix with imputed values and additional columns for missingness indicators.
#' @export
#' @examples
#' \donttest{
#' set.seed(123)
#' x_test <- matrix(rnorm(60), 10, 6)
#' colnames(x_test) <- paste0("X", 1:6)
#' x_test[1:2, 1] <- NA  # Introduce missing values
#' x_train <- matrix(rnorm(300), 50, 6)
#' colnames(x_train) <- paste0("X", 1:6)
#' model <- list(imputed = c(1))
#' attr(model, "imputed") <- c(1)
#' x_imputed <- impute_x_pred(model, x_test, x_train)
#' dim(x_imputed)  # 10 rows, 7 columns (6 original + 1 missingness indicator)
#' any(is.na(x_imputed))  # FALSE, no missing values
#' }
impute_x_pred <- function (object, x_test, x_train) {
  if (!is.null(attr(object, which = "imputed"))) {
    df <- data.frame(x_test)
    x_train <- data.frame(x_train)
    na.matr <- data.frame(1 * (is.na(df)))
    cm <- colMeans(na.matr)
    na.matr <- na.matr[, attr(object, which = "imputed")]
    names(na.matr) <- paste0("mis_", names(na.matr))
    for (i in which(cm != 0)){
      med <- median(x_train[[i]], na.rm = TRUE)
      if(is.na(med))
      {
        warning("One or more missing in test columns do not have any data in x_train, test set will be used for imputations!")
        med <-  median(df[[i]], na.rm = TRUE)
      }
      df[[i]][is.na(df[[i]])] <- med
    }
    return(as.matrix(data.frame(df,na.matr)))
  }
  return(as.matrix(x_test))
}
