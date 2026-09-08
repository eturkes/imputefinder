script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE))
setwd(dirname(dirname(normalizePath(script))))
rmarkdown::render("publication/manuscript.Rmd",knit_root_dir=getwd(),quiet=TRUE,
    envir=new.env(parent=globalenv()))
rmarkdown::render("publication/manuscript.Rmd",output_format=rmarkdown::word_document(),
    output_file="manuscript.docx",knit_root_dir=getwd(),quiet=TRUE,envir=new.env(parent=globalenv()))
