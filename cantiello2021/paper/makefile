.PHONY: main.pdf all clean

all: rednoise.pdf

rednoise.pdf: rednoise.tex refs.bib
	latexmk -pdf -pdflatex="pdflatex -interaction=nonstopmode" -use-make rednoise.tex

clean:
	latexmk -CA
