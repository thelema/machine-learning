open Batteries_uni
open Printf
open Predict

let test_data = Datafile.get_matrix "testing.txt.ba"

let () = 
  let infile = Sys.argv.(1) in
  let outall = infile ^ ".all" in
  let outbest = infile ^ ".best" in
  let predictor = pred_read Sys.argv.(1) |> predict_cat in
  printf "Predicting...%!";
  let t0 = Sys.time() in
  let preds = predictor test_data in
  printf "Done (%.2fs)\n%!" (Sys.time () -. t0);
  File.with_file_out outall (fun oc -> Array.iter (print_pred_conf oc) preds);
  let best = Array.enum preds |> best_predictions  in
  File.with_file_out outbest (fun oc -> output_preds oc best)
