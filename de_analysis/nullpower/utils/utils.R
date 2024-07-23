require(Rcpp)
require(Matrix)
require(data.table)
require(SingleCellExperiment)
require(dplyr)
require(tidyr)
require(magrittr)
require(Rcpp)

selectCol <- function(mat, j.col){
  x.col.dense <- rep(0,nrow(mat))
  p.begin <- mat@p[j.col]+1
  p.end <- mat@p[j.col+1]
  i.col <- mat@i[p.begin:p.end]+1 # i counts from 0
  x.col <- mat@x[p.begin:p.end]
  x.col.dense[i.col] <- x.col
  return(x.col.dense)
}

selectCols <- function(mat, j.cols){
  return(sapply(j.cols, selectCol, mat=mat))
}


src <-
  "
#include <Rcpp.h>

// [[Rcpp::export]]
void vec_down_sample(
    Rcpp::NumericVector data,
    const Rcpp::LogicalVector which,
    int begin,
    int end,
    double prob
    ){
    for(int i=begin; i<end; i++){
        if(which[i]){
            data[i] = data[i] + R::rbinom(data[i], prob);
            }
        }
    }
"
sourceCpp(code = src)

downCells <- function(spmat, i.rows, j.cols, p){
  data <- spmat@x
  i.bool <- spmat@i %in% (i.rows-1) # spmat@i begins from 0, i.rows begins from 1
  for (j.col in j.cols){
    begin <- spmat@p[j.col]
    end <- spmat@p[j.col+1]
    vec_down_sample(data, i.bool, begin, end, p)
  }
}

