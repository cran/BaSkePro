\name{BaSkePro}
\alias{BaSkePro}

\title{

Bayesian Model to Archaeological Faunal Skeletal Profiles
}
\description{

Tool to perform Bayesian inference of carcass processing/transport strategy and bone attrition from archaeological faunal skeletal profiles characterized by percentages of MAU (Minimum Anatomical Units). The approach is based on a generative model for skeletal profiles that replicates the two phases of formation of any archaeofaunal assemblage: initial accumulation as a function of human transport strategies and subsequent attrition.
}
\usage{
BaSkePro(x)
}

\arguments{
  \item{x}{A data frame object with the MAU (Minimum Anatomical Unit) percentage of one species in this specific order: Mandible, Atlas, Axis, Cervical vertebra, Thoracic vertebra, Lumbar vertebra, Rib, Sacrum, Scapula, Humerus, Radius (or radio-ulna), Metacarpal, Pelvis, Femur, Tibia, Metatarsal.

}
}
\details{    Two parameters define this model: 1) the transport preference (alpha), which can take any value between - 1 (mostly axial contribution) and 1 (mostly appendicular contribution) following strategies constructed as a function of butchering efficiency of different anatomical elements and the results of ethnographic studies, and 2) degree of attrition (beta), which can vary between 0 (no attrition) and 10 (maximum attrition) and relates the survivorship of bone elements to their maximum bone density.
    Starting from uniform prior probability distribution functions of alpha and beta, a Monte Carlo Markov Chain sampling based on a random walk Metropolis-Hasting algorithm is adopted to derive the posterior probability distribution functions, which are then available for interpretation. During this process, the likelihood of obtaining the observed percentages of MAU given a pair of parameter values is estimated by the inverse of the Chi2 statistic, multiplied by the proportion of elements within a 1 percent of the observed value.

}
\value{ A table and histograms with  posterior probability distribution functions of alpha and beta parameters (the first and fourth quartiles).

}
\references{Marín-Arroyo, A. B., & Ocio, D. (2018). Disentangling faunal skeletal profiles. A new probabilistic framework. Historical Biology, 30(6), 720-729. DOI:https://doi.org/10.1080/08912963.2017.1336620
}
\author{Ana B. Marin-Arroyo [aut], David Ocio [aut], Marco Vidal-Cordasco [cbt], Delphine Vettese [cbt]

}

\examples{

\donttest{PMAU <- as.data.frame(c(7,0,0,0,0,0,0,0,0,36,14,0,71,71,57,100))

BaSkePro (PMAU)}


}


