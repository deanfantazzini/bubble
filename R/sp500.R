#' Monthly SP500 index from 1791 till 2015 adjusted for dividends and splits.
#'
#' A dataset containing the monthly SP500 index from 1791 till 2015 adjusted for the dividends and splits.
#'
#' @docType data
#'
#' @usage data(sp500)
#'
#' @format A xts object with 2690 rows and 1 variable:
#' \describe{
#'   \item{price.adj}{SP500 index, ajusted for dividends and splits from 1791M8 till 2015M9}
#' }
#'
#' @importFrom zoo as.yearmon
#'
#' @references Aviral K. Tiwari, Arif B. Dar, Niyati Bhanja, and Rangan Gupta (2016). A Historical Analysis of the US Stock Price Index Using Empirical Mode Decomposition over 1791-2015. Economics Discussion Papers, No 2016-9, Kiel Institute for the World Economy
#'
#' @source \url{https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FZUQDM}
"sp500"
