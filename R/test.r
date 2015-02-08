library(roxygen2)

#' Test function
#'
#' This function is to test the roxygen2 markup language
#' @param  inputList
#' @export
#' @author Alexander Griffith <griffitaj@@gmail.com>
#' @examples
#' testFun(c(1,2,3,4))
testFun<-function(inputList)
    {
        for( i in inputList)
            print(i)
    }
