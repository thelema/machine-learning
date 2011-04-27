SOURCES = predict.ml gen_bigarray.ml test.ml train.ml best.ml hedge.ml cov_test.ml upgradeOVO.ml pred_test.ml best_many.ml

all: predict.native

$(SOURCES:.ml=.native): $(SOURCES)
	ocamlbuild $(SOURCES:.ml=.native)

prof: predict.p.native

$(SOURCES:.ml=.p.native): $(SOURCES)
	ocamlbuild $(SOURCES:.ml=.p.native)

byte: predict.byte

$(SOURCES:.ml=.byte): $(SOURCES)
	ocamlbuild $(SOURCES:.ml=.byte)

clean:
	ocamlbuild -clean