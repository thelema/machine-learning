open Batteries_uni
open Printf

open Lacaml.Impl.D

let perceptron n =
  let w = Vec.make0 n in
  let wsum = Vec.make0 n in
  let count = ref 0 in
  let wpred x = dot w x in
  let learn xi yi =
    let pred = wpred xi in
    if yi > 0. && pred < 0. || yi > 0. && pred < 0. then 
      add ~z:w ~x:x ~y:w;
    incr count;
    add ~z:wsum ~x:wsum ~y:w;
  in
  (learn, fun () -> let n = 1. /. float !count in map (fun i -> i *. n) wsum)

(*
let to_mat e = 
  let rows = Enum.count e in
  let cols = Enum.peek e |> Option.get |> List.length in
  let m = Mat.create rows cols in
  Enum.iteri (fun i l -> List.iteri (fun j x -> m.{i,j} <- float_of_string x) l) e;
  m
*)
let read_data_line line = String.nsplit line " "
let read_data_file fn = File.lines_of fn /@ read_data_line |> to_mat
let read_label_file fn = File.lines_of fn /@ int_of_string 

let () = printf "Reading..%!"
let t0 = Sys.time()
let alldata = read_data_file "development.txt"
let alllabels = read_label_file "development_label.txt"
let () = printf "Done(%.2f)\n%!" (Sys.time () -. t0)

(*
let vectors_train = 
let labels_train = 
let vectors_test = read_data_file "test.txt"
let labels_test = read_label_file "test_label.txt"
*)



let () = printf "Reading..%!"
let t0 = Sys.time()
let v1 = svm_scale ~lower:0. alldata
let m = svm_train ~kernel_type:LINEAR v1
let () = printf "Done(%.2f)\n%!" (Sys.time () -. t0)
