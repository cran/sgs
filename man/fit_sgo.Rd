% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_sgo.R
\name{fit_sgo}
\alias{fit_sgo}
\title{Fit an SGO model.}
\usage{
fit_sgo(
  X,
  y,
  groups,
  type = "linear",
  lambda = "path",
  path_length = 20,
  min_frac = 0.05,
  alpha = 0.95,
  max_iter = 5000,
  backtracking = 0.7,
  max_iter_backtracking = 100,
  tol = 1e-05,
  standardise = "l2",
  intercept = TRUE,
  screen = TRUE,
  verbose = FALSE,
  w_weights = NULL,
  v_weights = NULL,
  warm_start = NULL
)
}
\arguments{
\item{X}{Input matrix of dimensions \eqn{n \times p}{n*p}. Can be a sparse matrix (using class \code{"sparseMatrix"} from the \code{Matrix} package).}

\item{y}{Output vector of dimension \eqn{n}. For \code{type="linear"} should be continuous and for \code{type="logistic"} should be a binary variable.}

\item{groups}{A grouping structure for the input data. Should take the form of a vector of group indices.}

\item{type}{The type of regression to perform. Supported values are: \code{"linear"} and \code{"logistic"}.}

\item{lambda}{The regularisation parameter. Defines the level of sparsity in the model. A higher value leads to sparser models:
\itemize{
\item \code{"path"} computes a path of regularisation parameters of length \code{"path_length"}. The path will begin just above the value at which the first predictor enters the model and will terminate at the value determined by \code{"min_frac"}.
\item User-specified single value or sequence. Internal scaling is applied based on the type of standardisation. The returned \code{"lambda"} value will be the original unscaled value(s).
}}

\item{path_length}{The number of \eqn{\lambda} values to fit the model for. If \code{"lambda"} is user-specified, this is ignored.}

\item{min_frac}{Smallest value of \eqn{\lambda} as a fraction of the maximum value. That is, the final \eqn{\lambda} will be \code{"min_frac"} of the first \eqn{\lambda} value.}

\item{alpha}{The value of \eqn{\alpha}, which defines the convex balance between OSCAR and gOSCAR. Must be between 0 and 1. Recommended value is 0.95.}

\item{max_iter}{Maximum number of ATOS iterations to perform.}

\item{backtracking}{The backtracking parameter, \eqn{\tau}, as defined in Pedregosa and Gidel (2018).}

\item{max_iter_backtracking}{Maximum number of backtracking line search iterations to perform per global iteration.}

\item{tol}{Convergence tolerance for the stopping criteria.}

\item{standardise}{Type of standardisation to perform on \code{X}:
\itemize{
\item \code{"l2"} standardises the input data to have \eqn{\ell_2} norms of one. When using this \code{"lambda"} is scaled internally by \eqn{1/\sqrt{n}}.
\item \code{"l1"} standardises the input data to have \eqn{\ell_1} norms of one. When using this \code{"lambda"} is scaled internally by \eqn{1/n}.
\item \code{"sd"} standardises the input data to have standard deviation of one.
\item \code{"noBaone"} no standardisation applied.
}}

\item{intercept}{Logical flag for whether to fit an intercept.}

\item{screen}{Logical flag for whether to apply screening rules (see Feser and Evangelou (2024)). Screening discards irrelevant groups before fitting, greatly improving speed.}

\item{verbose}{Logical flag for whether to print fitting information.}

\item{w_weights}{Optional vector for the group penalty weights. Overrides the OSCAR penalties when specified. When entering custom weights, these are multiplied internally by \eqn{\lambda} and \eqn{1-\alpha}. To void this behaviour, set \eqn{\lambda = 2} and \eqn{\alpha = 0.5}.}

\item{v_weights}{Optional vector for the variable penalty weights. Overrides the OSCAR penalties when specified. When entering custom weights, these are multiplied internally by \eqn{\lambda} and \eqn{\alpha}. To void this behaviour, set \eqn{\lambda = 2} and \eqn{\alpha = 0.5}.}

\item{warm_start}{Optional list for implementing warm starts. These values are used as initial values in the fitting algorithm. Need to supply \code{"x"} and \code{"u"} in the form \code{"list(warm_x, warm_u)"}. Not recommended for use with a path or CV fit as start from the null model by design.}
}
\value{
A list containing:
\item{beta}{The fitted values from the regression. Taken to be the more stable fit between \code{x} and \code{z}, which is usually the former. A filter is applied to remove very small values, where ATOS has not been able to shrink exactly to zero. Check this against \code{x} and \code{z}.}
\item{group_effects}{The group values from the regression. Taken by applying the \eqn{\ell_2} norm within each group on \code{beta}.}
\item{selected_var}{A list containing the indicies of the active/selected variables for each \code{"lambda"} value. Index 1 corresponds to the first column in X.}
\item{selected_grp}{A list containing the indicies of the active/selected groups for each \code{"lambda"} value. Index 1 corresponds to the first group in the \code{groups} vector. You can see the group order by running \code{unique(groups)}.}
\item{num_it}{Number of iterations performed. If convergence is not reached, this will be \code{max_iter}.}
\item{success}{Logical flag indicating whether ATOS converged, according to \code{tol}.}
\item{certificate}{Final value of convergence criteria.}
\item{x}{The solution to the original problem (see Pedregosa and Gidel (2018)).}
\item{z}{The updated values from applying the first proximal operator (see Pedregosa and Gidel (2018)).}
\item{u}{The solution to the dual problem (see Pedregosa and Gidel (2018)).}
\item{screen_set_var}{List of variables that were kept after screening step for each \code{"lambda"} value. (corresponds to \eqn{\mathcal{S}_v} in Feser and Evangelou (2024)).}
\item{screen_set_grp}{List of groups that were kept after screening step for each \code{"lambda"} value. (corresponds to \eqn{\mathcal{S}_g} in Feser and Evangelou (2024)).}
\item{epsilon_set_var}{List of variables that were used for fitting after screening for each \code{"lambda"} value. (corresponds to \eqn{\mathcal{E}_v} in Feser and Evangelou (2024)).}
\item{epsilon_set_grp}{List of groups that were used for fitting after screening for each \code{"lambda"} value. (corresponds to \eqn{\mathcal{E}_g} in Feser and Evangelou (2024)).}
\item{kkt_violations_var}{List of variables that violated the KKT conditions each \code{"lambda"} value. (corresponds to \eqn{\mathcal{K}_v} in Feser and Evangelou (2024)).}
\item{kkt_violations_grp}{List of groups that violated the KKT conditions each \code{"lambda"} value. (corresponds to \eqn{\mathcal{K}_g} in Feser and Evangelou (2024)).}
\item{pen_slope}{Vector of the variable penalty sequence.}
\item{pen_gslope}{Vector of the group penalty sequence.}
\item{screen}{Logical flag indicating whether screening was performed.}
\item{type}{Indicates which type of regression was performed.}
\item{intercept}{Logical flag indicating whether an intercept was fit.}
\item{lambda}{Value(s) of \eqn{\lambda} used to fit the model.}
}
\description{
Sparse-group OSCAR (SGO) main fitting function. Supports both linear and logistic regression, both with dense and sparse matrix implementations.
}
\details{
\code{fit_sgo()} fits an SGO model (Feser and Evangelou (2024)) using adaptive three operator splitting (ATOS). SGO uses the same model set-up as for SGS, but with different weights (see Bao et al. (2020) and Feser and Evangelou (2024)).
The penalties are given by (for a group \eqn{g} and variable \eqn{i}, with \eqn{p} variables and \eqn{m} groups):
\deqn{
  v_i = \sigma_1 + \sigma_2(p-i), \; w_g = \sigma_1 + \sigma_3(m-g),
}
where
\deqn{
     \sigma_1 = d_i\|X^\intercal y\|_\infty, \; \sigma_2 = \sigma_1/p, \; \sigma_3 = \sigma_1/m, \; d_i = i \times \exp{(-2)}. 
}
}
\examples{
# specify a grouping structure
groups = c(1,1,1,2,2,3,3,3,4,4)
# generate data
data =  gen_toy_data(p=10, n=5, groups = groups, seed_id=3,group_sparsity=1)
# run SGO
model = fit_sgo(X = data$X, y = data$y, groups = groups, type="linear", path_length = 5, 
alpha=0.95, standardise = "l2", intercept = TRUE, verbose=FALSE)
}
\references{
Bao, R., Gu B., Huang, H. (2020). \emph{Fast OSCAR and OWL Regression via Safe Screening Rules}, \url{https://proceedings.mlr.press/v119/bao20b}

Feser, F., Evangelou, M. (2023). \emph{Sparse-group SLOPE: adaptive bi-level selection with FDR-control}, \url{https://arxiv.org/abs/2305.09467}

Feser, F., Evangelou, M. (2024). \emph{Strong screening rules for group-based SLOPE models}, \url{https://arxiv.org/abs/2405.15357}

Pedregosa, F., Gidel, G. (2018). \emph{Adaptive Three Operator Splitting}, \url{https://proceedings.mlr.press/v80/pedregosa18a.html}
}
\seealso{
Other SGS-methods: 
\code{\link{as_sgs}()},
\code{\link{coef.sgs}()},
\code{\link{fit_sgo_cv}()},
\code{\link{fit_sgs}()},
\code{\link{fit_sgs_cv}()},
\code{\link{plot.sgs}()},
\code{\link{predict.sgs}()},
\code{\link{print.sgs}()},
\code{\link{scaled_sgs}()}
}
\concept{SGS-methods}
