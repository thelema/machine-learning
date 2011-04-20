open Batteries_uni
open Printf
open Bigarray
open Lacaml.Impl.S (* Single-precision reals *)

open Predict

let n = 100

let () =
  let test fn = 
    let p = pred_read fn |> predict_cat in
    let cs = group_by_cat n in

    let correct = ref 0 in
    let total = ref 0 in
    let test_xi i xi = 
      incr total; if !total land 0xff = 0 then printf ".%!"; 
      if snd (p (get_i train_data xi)) = i then incr correct
    in
    Array.iteri (fun i csi -> Array.iter (test_xi i) csi) cs;
    let total = Array.fold_left (fun acc csi -> acc + Array.length csi) 0 cs in
    printf "%s: %d correct of %d\n" fn !correct total;
  in
  args () |> iter test
