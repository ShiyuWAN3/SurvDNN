\name{mod_permfit}
\alias{mod_permfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 mod_permfit to do ...
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mod_permfit(method, model.type, object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{
%%     ~~Describe \code{x} here~~
  Different methods of machine learning
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Methods implemeted in mod_permfit:

"random_forest": Random Survival Forest;

"survival_aft": accelerated failure model;

"survival_cox": Cox-PH regression;

"DeepSurv": deep learning based DeepSurv;

"DeepHit": deep learning based DeepHit;

"Xgboost": Xgboost with Cox-based loss;

"Survival_SVM": survival support vector machine;

"lasso": Lasso-Cox;

"ensemble_dnnet": proposed SurvDNN
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Shiyu Wan
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
