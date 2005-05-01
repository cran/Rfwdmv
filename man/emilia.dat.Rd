\name{emilia.dat}
\alias{emilia.dat}
\docType{data}

\title{Municipalities in Emilia-Romagna}

\description{
A large data set with 341 observations on 29 variables.  The data are taken from the 1991 Italian census and cover all municipalities.  Nearly all the variables are indices in which counts have been divided by municipal population.
}

\usage{data(emilia.dat)}
\format{
  A data frame with 341 observations on the following 29 variables.

\describe{
    \item{Province}{a factor with levels \code{BO}, \code{FE}, \code{FO}, \code{MO}, \code{PC}, \code{PR}, \code{RA}, \code{RE} and \code{RN} identifying the province.}
    \item{Pop.Inf}{a numeric vector, the percentage of the population aged less than ten (y1 in ARC).}
    \item{Pop.sen}{a numeric vector, the percentage of the population aged 75 or more (y2 in ARC).}
    \item{Unipers}{a numeric vector, the percentage of the single-member families (y3 in ARC).}
    \item{Divorzi}{a numeric vector, the percentage of the residents who are divorced (y4 in ARC).}
    \item{Vedovi}{a numeric vector, the percentage of widows and widowers (y5 in ARC).}
    \item{Laurea}{a numeric vector, the percentage of the population aged over 25 (y6 in ARC).}
    \item{Notitolo}{a numeric vector, the percentage of those aged over six having no education (y7 in ARC).}
    \item{T.attiv}{a numeric vector, the percentage of those of working age in full-time employment (y8 in ARC).}
    \item{T.disoc.}{a numeric vector, the unemployment rate (y9 in ARC).}
    \item{Mov.nat}{a numeric vector, the standardized natural increase in population (y10 in ARC).}
    \item{Mov.mig}{a numeric vector, the standardized change in population due to migration (y11 in ARC).}
    \item{Natalit.}{a numeric vector the average birth rate from 1992 to 1994 (y12 in ARC).}
    \item{Fecondi1}{a numeric vector, fecundity: the three year average birth rate amongst women of child-bearing age (y13 in ARC).}
    \item{Dopo.82}{a numeric vector, the percentage of occupied houses built since 1982 (y14 in ARC).}
    \item{xbagni}{a numeric vector, the percentage of occupied houses with two or more WCs (y15 in ARC).}
    \item{Impianto}{a numeric vector, the percentage of occupied houses with fixed heating system (y16 in ARC).}
    \item{TV..}{a numeric vector, the percentage of TV licence holders (y17 in ARC).}
    \item{Parco}{a numeric vector, the number of cars per 100 inhabitants (y18 in ARC).}
    \item{Lusso.Imm}{a numeric vector, the percentage of luxury cars (y19 in ARC).}
    \item{Add.H}{a numeric vector, the percentage of those working in hotels and restaurants (y20 in ARC).}
    \item{Add.J}{a numeric vector, the percentage of those working in banking and finance (y21 in ARC).}
    \item{Imp.cont}{a numeric vector, the average decleared income amongst those filing income tax returns (y22 in ARC).}
    \item{Cont.res}{a numeric vector, the percentage of inhabitants filing income tax returns (y23 in ARC).}
    \item{Add.Ul}{a numeric vector, the percentage of the residents employed in factories and public services (y24 in ARC).}
    \item{x0add.Add}{a numeric vector, the percentage of employees working in factories with more than ten employees (y25 in ARC).}
    \item{V68}{a numeric vector, the percentage of employees working in factories with more than 50 employees (y26 in ARC).}
    \item{Artigiane}{a numeric vector, the percentage of those working in artisanal enterprises (y27 in ARC).}
    \item{Imprendi}{a numeric vector, the percentage enterpreneurs and skilled self-emplyed among those of working age (y28 in ARC).}
  }
}
\details{
These 29 variables were selected from 50 available.  The first 13 (\code{y1}-\code{y13}) are demographic variables, the next three (\code{y14}-\code{y16}) measure housing quality, the succeeding seven (\code{y17}-\code{y23}) are measures of individual income and wealth and the last five (\code{y24}-\code{y28}) relate to industrial production. 
}

\source{
Zani (1996).
}

\references{
A selection of these data is given by Atkinson, Riani and Cerioli (2004), p. 560-561; the full data are available on the website http://www.riani.it/arc.
}

\examples{
data(emilia.dat)
}

\keyword{datasets}
