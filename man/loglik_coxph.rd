% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/source_permfit_survival_20230610.r
\name{loglik_coxph}
\alias{loglik_coxph}
\title{Title Calculate the Cox-PH Log-Likelihood Based on Event Status, Event Time, and Risk Prediction}
\usage{
loglik_coxph(Status, Times, f_hat_y)
}
\arguments{
\item{Status}{Patients' event statuses.}

\item{Times}{Patients' event times.}

\item{f_hat_y}{The product of X and Beta (X \%*\% Beta).}
}
\value{
Returns the log-likelihood of a Cox-PH model.
}
\description{
Title Calculate the Cox-PH Log-Likelihood Based on Event Status, Event Time, and Risk Prediction
}
