all: predict.native

predict.native: predict.ml
	ocamlbuild predict.native

clean:
	ocamlbuild -clean