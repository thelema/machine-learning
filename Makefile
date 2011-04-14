SOURCES = predict.ml gen_bigarray.ml test.ml train.ml best.ml

all: predict.native

$(SOURCES:.ml=.native): $(SOURCES)
	ocamlbuild $(SOURCES:.ml=.native)

prof: predict.p.native

$(SOURCES:.ml=.p.native): $(SOURCES)
	ocamlbuild $(SOURCES:.ml=.p.native)


clean:
	ocamlbuild -clean