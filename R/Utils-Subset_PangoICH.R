
# R/Utils-Subset_PangoICH.R

#' @export

Subset_PangoICH <- function(x,
                            i,
                            first,
                            last,
                            except)
{

  UseMethod("Subset_PangoICH")

}

#' @export

Subset_PangoICH.list <- function(x,
                                 i = 1,
                                 first = FALSE,
                                 last = FALSE,
                                 except = FALSE)
  {

  if (first) {

    if (except) {

      return(x[[-1]])

    } else {

      return(x[[1]])

    }

  }

  if (last) {

    if (except) {

      return(x[[-length(x)]])

    } else {

      return(x[[length(x)]])

    }

  }

  if (except) {

    return(x[[-i]])

  } else {

    return(x[[i]])

  }

}

#' @export

Subset_PangoICH.character <- function(x,
                                      i = 1,
                                      first = FALSE,
                                      last = FALSE,
                                      except = FALSE)
  {

  if (first) {

    if (except) {

      return(x[-1])

    } else {

      return(x[1])

    }
  }

  if (last) {

    if (except) {

      return(x[-length(x)])

    } else {

      return(x[length(x)])

    }
  }

  if (except) {

    return(x[-i])

  } else {

    return(x[i])

  }
}
