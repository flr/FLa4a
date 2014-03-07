PD = FLa4a
SD = $(PD)/inst/admb
BD = $(PD)/inst/bin/linux
LD = /home/millaco/R/x86_64-pc-linux-gnu-library/3.0/FLa4a/R
sourcefiles := $(wildcard $(PD)/R/*.R)

.PHONY = roxygen install compile

all: install

compile: $(BD)/a4a
$(BD)/a4a: $(SD)/a4a.tpl $(SD)/nLogNormal.h 
	rm -rf _tmp
	mkdir _tmp
	cp $(SD)/a4a.tpl $(SD)/nLogNormal.h _tmp
	cd _tmp; admb -s a4a 
	cp _tmp/a4a $(BD)
	rm -rf _tmp

roxygen: $(PD)/man/FLa4a-package.Rd
$(PD)/man/FLa4a-package.Rd: $(sourcefiles)
	rm -rf $(PD)/man/*
	echo 'library(roxygen2); roxygenize("FLa4a", roclets = c("namespace","rd"))' | R --vanilla --slave

install: $(LD)/FLa4a 
$(LD)/FLa4a: $(PD)/man/FLa4a-package.Rd $(BD)/a4a
	R CMD INSTALL FLa4a

