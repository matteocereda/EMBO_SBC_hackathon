## ***************************
##
## Script name: environment.R
## Purpose of script: config file to set environment variables used in tbe oproject 
## Author: mc
## Date Created: 2023-11-13
##
## ***************************
##
## Notes:
##   
## ***************************




if( Sys.info()['nodename']=='psychobook-2.local'){
  
}else if( Sys.info()['nodename']=='lab07.local'){
  DATA_DIR='/Volumes/mc/repository/EMBO_hackathon'
  CGB_DIR='/Volumes/mc/'
  GIT_LOCAL_DIR='/Volumes/mc/repository/'
  CGB_SHARED='~/Dropbox (HuGeF)/'
  BEDTOOLS='/Applications/miniconda3/bin/bedtools'

}else if( Sys.info()['user'] =='mariachiara.grieco' ){
  DATA_DIR='/hpcnfs/data/cgb/EMBO_hackathon/'
  CGB_DIR='/hpcnfs/data/cgb/'
  GIT_LOCAL_DIR='/hpcnfs/home/mariachiara.grieco/'
  
}else if( Sys.info()['user'] == 'mgrieco'){
  DATA_DIR='/Users/mgrieco/repo/EMBO_hackathon'
  CGB_DIR='/Users/mgrieco/'
  GIT_LOCAL_DIR='/Users/mgrieco/repo/'
  CGB_SHARED='~/Dropbox (HuGeF)/'
  BEDTOOLS='/Users/mgrieco/miniconda/bin/bedtools'
  
}else if( Sys.info()['user'] == 'mgrieco'){
  DATA_DIR='/Users/mgrieco/repo/EMBO_hackathon'
  CGB_DIR='/Users/mgrieco/'
  GIT_LOCAL_DIR='/Users/mgrieco/repo/'
  CGB_SHARED='~/Dropbox (HuGeF)/'
  BEDTOOLS='/Users/mgrieco/miniconda/bin/bedtools'
  
}else if( Sys.info()['nodename'] == 'cgb-tower'){
  DATA_DIR='/home/matteo/Lavoro/EMBO_hackathon'
  CGB_DIR='/home/matteo/Lavoro/'
  GIT_LOCAL_DIR='/home/matteo/Lavoro/'
  CGB_SHARED='/home/matteo/Dropbox (HuGeF)/'
  BEDTOOLS='/usr/bin/bedtools'
  
  
  
  
}else{
}


