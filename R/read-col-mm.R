#' readColMM
#'
#' Modified the `Matrix::readMM` function for reading matrices stored in the
#' Harwell-Boeing or MatrixMarket formats but only reading selected column.
#'
#' See \code{\link[Matrix]{readMM}}
#'
#' @param file the name of the file to be read from as a character scalar. Those
#'  storing matrices in the MatrixMarket format usually end in ".mtx".
#' @param which.col An integer scalar, the column index
#' @param chunk, An integer scalar indicating the chunk size to use,
#' i.e., number of rows to read at any one time.
#' @author Ruqian Lyu
#' @export
#' @examples
#' demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
#' readColMM(file = paste0(demo_path,"s1_chr1_vi.mtx"), which.col=2,chunk=2)
#' @return
#' A sparse matrix object that inherits from the "Matrix" class which the original
#' dimensions.To get the vector of the specified column, one need to subset
#' the matrix to select the column with the same index.
#'
readColMM <- function(file,which.col,chunk=1000L)
{
  if (is.character(file))
    file <- if(file == "") stdin() else file(file)
  if (!inherits(file, "connection"))
    stop("'file' must be a character string or connection")
  if (!isOpen(file)) {
    open(file)
    on.exit(close(file))
  }
  scan1 <- function(what, ...)
    scan(file, nmax = 1, what = what, quiet = TRUE, ...)

  if (scan1(character()) != "%%MatrixMarket")# hdr
    stop("file is not a MatrixMarket file")
  if (!(typ <- tolower(scan1(character()))) %in% "matrix")
    stop(gettextf("type '%s' not recognized", typ), domain = NA)
  if (!(repr <- tolower(scan1(character()))) %in% c("coordinate", "array"))
    stop(gettextf("representation '%s' not recognized", repr), domain = NA)
  elt <- tolower(scan1(character()))
  if (!elt %in% c("real", "complex", "integer", "pattern"))
    stop(gettextf("element type '%s' not recognized", elt), domain = NA)

  sym <- tolower(scan1(character()))
  if (!sym %in% c("general", "symmetric", "skew-symmetric", "hermitian"))
    stop(gettextf("symmetry form '%s' not recognized", sym), domain = NA)
  nr <- scan1(integer(), comment.char = "%")
  nc <- scan1(integer())
  nz <- scan1(integer())
  checkIJ <- function(els) {

    if(any(els$i < 1 | els$i > nr))
      stop("readMM(): row	 values 'i' are not in 1:nr", call.=FALSE)
    if(any(els$j < 1 | els$j > nc))
      stop("readMM(): column values 'j' are not in 1:nc", call.=FALSE)
  }
  if (repr == "coordinate") {
    switch(elt,
           "real" = ,
           "integer" = {
             ## TODO: the "integer" element type should be returned as
             ##       an object of an "iMatrix" subclass--once there are


             # Reading it in, chunk by chunk (see behavior of nmax= when what=
             # is a list).
             els <- list(i=NULL,j=NULL,x=NULL)
             repeat {
               current <- scan(file, what=list(i= integer(),
                                               j= integer(),
                                               x= numeric()),
                               nmax=chunk, quiet=TRUE)
               checkIJ(current)
               els <- list(i = c(els$i,current$i[current$j==which.col]),
                           j = c(els$j,current$j[current$j==which.col]),
                           x = c(els$x,current$x[current$j==which.col]))
               if ( length(current$i) < chunk) {
                 break
               }
             }

             switch(sym,
                    "general" = {
                      new("dgTMatrix", Dim = c(nr, nc), i = els$i - 1L,
                          j = els$j - 1L, x = els$x)
                    },
                    "symmetric" = {
                      new("dsTMatrix", uplo = "L", Dim = c(nr, nc),
                          i = els$i - 1L, j = els$j - 1L, x = els$x)
                    },
                    "skew-symmetric" = {
                      stop("symmetry form 'skew-symmetric' not yet implemented
                           for reading")
                      ## FIXME: use dgT... but must expand the (i,j,x) slots!
                      new("dgTMatrix", uplo = "L", Dim = c(nr, nc),
                          i = els$i - 1L, j = els$j - 1L, x = els$x)

                    },
                    "hermitian" = {
                      stop("symmetry form 'hermitian' not yet implemented for
                           reading")
                    },
                    ## otherwise (not possible; just defensive programming):
                    stop(gettextf("symmetry form '%s' is not yet implemented",
                                  sym), domain = NA)
             )
           },
           "pattern" = {
             els <- scan(file, nmax = nz, quiet = TRUE,
                         what = list(i = integer(), j = integer()))
             checkIJ(els)
             switch(sym,
                    "general" = {
                      new("ngTMatrix", Dim = c(nr, nc),
                          i = els$i - 1L, j = els$j - 1L)
                    },
                    "symmetric" = {
                      new("nsTMatrix", uplo = "L", Dim = c(nr, nc),
                          i = els$i - 1L, j = els$j - 1L)
                    },
                    "skew-symmetric" = {
                      stop("symmetry form 'skew-symmetric' not yet implemented
                           for reading")
                      ## FIXME: use dgT... but must expand the (i,j,x) slots!
                      new("ngTMatrix", uplo = "L", Dim = c(nr, nc),
                          i = els$i - 1L, j = els$j - 1L)

                    },
                    "hermitian" = {
                      stop("symmetry form 'hermitian' not yet implemented for
                           reading")
                    },
                    ## otherwise (not possible; just defensive programming):
                    stop(gettextf("symmetry form '%s' is not yet implemented",
                                  sym), domain = NA)
             )
           },
           "complex" = {
             stop("element type 'complex' not yet implemented")
           },
           ## otherwise (not possible currently):
           stop(gettextf("'%s()' is not yet implemented for element type '%s'",
                         "readMM", elt), domain = NA))
  }
  else
    stop(gettextf("'%s()' is not yet implemented for  representation '%s'",
                  "readMM", repr), domain = NA)
}
