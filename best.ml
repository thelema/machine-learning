open Batteries_uni
open Printf
open Predict
open Bigarray

let test_data = Datafile.get_matrix "testing.txt.ba"

let () = 
  let infile = Sys.argv.(1) in
  let outall = infile ^ ".all" in
  let outbest = infile ^ ".best" in
  let predictor = pred_read infile |> predict_cat in
  printf "Predicting%!";
  let t0 = Sys.time() in
  let n = if Array.length Sys.argv > 2 then int_of_string Sys.argv.(2) else -1 in
  let test_data = if n = -1 then test_data else Array2.sub_right test_data 1 n in
  let preds = predict_many predictor test_data in
  printf "Done (%.2fs)\n%!" (Sys.time () -. t0);
  File.with_file_out outall (fun oc -> Array.iter (print_pred_conf oc) preds);
  let best = Array.enum preds |> best_predictions  in
  File.with_file_out outbest (fun oc -> output_preds oc best)
