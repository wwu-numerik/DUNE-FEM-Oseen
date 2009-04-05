TEX = latex
TEX_SRC = $(shell grep '%%% TeX-master' *.tex  | cut -d':' -f1)
TARGET = femhowto.dvi
TEXFLAGS = --src-specials
# TEXFLAGS += \\nonstopmode\\input

all: $(TARGET)

%.dvi: %.tex
	$(TEX) $(TEXFLAGS) $<;
	if grep $*.bbl $*.log > /dev/null; then \
	  bibtex $*; \
	fi;
	while grep Rerun $*.log > /dev/null ; do \
          $(TEX) $(TEXFLAGS) $< ; \
        done

%.ps: %.dvi
	dvips $< -o

clean:
	rm -rf auto
	rm -f *~
	rm -f *.aux *.bbl *.blg *.log *.toc
	rm -f  \#* *.cb *.idx *.not *.rel TAGS

distclean: clean
	rm -f *.dvi *.ps *.pdf
