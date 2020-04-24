.PHONY: main.pdf all clean

all: paper.pdf

figures/schematic-crop.pdf: figures/schematic.pdf
	cd figures; pdfcrop schematic.pdf; cd ..

paper.pdf: paper.tex refs.bib figures/schematic-crop.pdf
	latexmk -pdf -pdflatex="pdflatex -interaction=nonstopmode" -use-make paper.tex

clean:
	latexmk -CA
