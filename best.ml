open Batteries_uni
open Printf
open Predict
open Bigarray


let test_data = Datafile.get_matrix "testing.txt.ba"
let batch_size = 10000

let () = 
  let make_best infile = 
    let outall = infile ^ ".all" in
    let predictor = pred_read infile in
    printf "Predicting%!";
    let t0 = Sys.time() in
    let n = Array2.dim2 test_data in
    let batches = n / batch_size in
    let oc = File.open_out outall in
    for i = 0 to batches do (* extra loop at the end for any leftover *)
      let i0 = i * batch_size in
      let slice = Array2.sub_right test_data (i0+1) (i0+batch_size) in
      batch_predict_cat predictor slice |> Array.iter (print_pred_conf oc);
      printf ",";
    done;
    IO.close_out oc;
    printf "Done (%.2fs)\n%!" (Sys.time () -. t0);
  in
  args () |> iter make_best
