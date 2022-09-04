#' Chemotherapy for Stage B/C colon cancer
#'
#' These are data from one of the first successful trials of adjuvant
#' chemotherapy for colon cancer. Levamisole is a low-toxicity compound
#' previously used to treat worm infestations in animals; 5-FU is a moderately
#' toxic (as these things go) chemotherapy agent. There are only one record per patient for the death outcome (or censoring). This is redistributed from
#' the survival package, with a small modification to include only the death outcome.
#'
#' @format A data frame with 929 rows and 17 variables:
#' \describe{
#'
#' \item{id}{id}
#' \item{study}{1 for all patients}
#' \item{rx}{Treatment - Obs(ervation), Lev(amisole), Lev(amisole)+5-FU}
#' \item{sex}{1=male}
#' \item{age}{in years}
#' \item{obstruct}{obstruction of colon by tumour}
#' \item{perfor}{perforation of colon}
#' \item{adhere}{adherence to nearby organs}
#' \item{nodes}{number of lymph nodes with detectable cancer}
#' \item{time}{days until death or censoring}
#' \item{status}{censoring status}
#' \item{differ}{differentiation of tumour (1=well, 2=moderate, 3=poor)}
#' \item{extent}{Extent of local spread (1=submucosa, 2=muscle, 3=serosa, 4=contiguous structures)}
#' \item{surg}{time from surgery to registration (0=short, 1=long)}
#' \item{node4}{more than 4 positive lymph nodes}
#' \item{etype}{event type: 1=recurrence,2=death}
#' \item{event}{event indicator: censored, death}
#' }
#' @seealso \link[survival]{colon}
"colon"


#' Monoclonal gammopathy data
#'
#' Natural history of 1341 sequential patients with monoclonal gammopathy of undetermined significance (MGUS). This is a superset of the mgus data, at a later point in the accrual process. This dataset is redistributed from the survival package with an added competing risks event indicator.
#'
#' @format A data frame with 1384 observations on the following 10 variables.
#' \describe{
#'     \item{\code{id}}{subject identifier}
#'     \item{\code{age}}{age at diagnosis, in years}
#'     \item{\code{sex}}{a factor with levels \code{F} \code{M}}
#'     \item{\code{dxyr}}{year of diagnosis}
#'     \item{\code{hgb}}{hemoglobin}
#'     \item{\code{creat}}{creatinine}
#'     \item{\code{mspike}}{size of the monoclonal serum spike}
#'     \item{\code{ptime}}{time until progression to a plasma cell
#'         malignancy (PCM) or last contact, in months}
#'     \item{\code{pstat}}{occurrence of PCM: 0=no, 1=yes }
#'     \item{\code{futime}}{time until death or last contact, in months}
#'     \item{\code{death}}{occurrence of death: 0=no, 1=yes}
#'     \item{\code{etime}}{time until either death, pcm, or last contact}
#'     \item{\code{event}}{factor indicating which event occurred first}
#' }
#' @seealso \link[survival]{mgus2}
"mgus2"

