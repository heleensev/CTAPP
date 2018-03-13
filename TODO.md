
#### ¿BUGS?
_____________________________________________
#### MongoParsenStore

###### liftover_creator
  * ~~remove globals, client flow like reference_creator~~


###### reference_creator
  * ~~find solution for db (bulk) insert in MultiAlleleException~~

#### MetaReader
  * ~~incorporate mon_vmem~~

##### GWASParse

###### meta_analysis
 * db for bp_max should be ref_db not GWASdb

###### column identifier
 * ~~what to do with duplicate columns~~
 * ~~duplicate header handling (in new_headers/headers)~~

###### val checker
 * ~~sanity check, what happens with effect allele column~~

###### log_writer
 * ~~fix loop in error_logger~~

#### Runner

###### map_chunker
  * double quote list in subprocess
  * double quote \t or escape characters

###### configparser
  * ConfigError does not print message (fixed?)

###### run.py
  * ~~subparallel passes pandas Series to function liftover~~

##### Preprocessor
  * 00000001 prefixed in one of GWAS chunks

##### JobHandler

###### mongo_handler


###### job_computer
  * ~~skip does not go to default 0~~

____________________________________________

#### GOALS
_____________________________________________

##### MetaReader
 * ~~incorporate validator into MetaReader~~
 * ~~new Study object builder for new config~~
 * job config writer module
 * validator: n_col and n_studies
 * parameter for preparse or postparse config validator schema
 * revalidate module

##### GWASParse
 * perform check of header indices if user updated config
 * write errs of study object in log file
 * method for just head printing
 * check study build with hard coded piece
 * include more in regex and accept missing chr and bp column:
 https://huwenboshi.github.io/data%20management/2017/11/23/tips-for-formatting-gwas-summary-stats.html

##### ErrorLogger
 * write head with newly classified headers for successful studies
 * write head with original columns for failed studies?

##### JobHandler

###### mongo_handler
  * ~~use db.currentop() command before db kill~~

#### Runner

##### mapchunker
 * general method for setting attrs

##### drugbankDB
 * pharmalogically active: genbank protein ID matches with 1000g reference file
 * fulldb xml: <atc- codes> keywords
 * https://docs.drugbankplus.com/v1/#introduction

