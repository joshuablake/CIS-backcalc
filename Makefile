SHARED_DEPS = utils.R data/STATS18800/predict.rds data/STATS18800/groups.csv data/STATS18596/poststrat.csv

.PHONY: clean

all: outputs/01_plot_prev.html outputs/deconv.rds outputs/groups.csv outputs/poststrat.csv outputs/predict.rds outputs/README.txt

outputs/01_plot_prev.html: 01_plot_prev.Rmd $(SHARED_DEPS)
	Rscript -e 'rmarkdown::render("01_plot_prev.Rmd", output_dir = "outputs", quiet = TRUE)'

outputs/deconv.rds: 02_do_backcalc.R $(SHARED_DEPS)
	Rscript 02_do_backcalc.R

outputs/groups.csv: data/STATS18800/groups.csv
	cp $< $@

outputs/poststrat.csv: data/STATS18596/poststrat.csv
	cp $< $@

outputs/predict.rds: data/STATS18800/predict.rds
	cp $< $@

outputs/README.txt: data/STATS18800/README.txt
	cp $< $@

clean:
	rm -f outputs/*