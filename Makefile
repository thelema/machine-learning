all: predict.native

predict.native gen_bigarray.native: predict.ml gen_bigarray.ml
	ocamlbuild predict.native gen_bigarray.native

clean:
	ocamlbuild -clean