open Batteries_uni
open Printf
open Predict
open Bigarray

let test_data = Datafile.get_matrix "testing.txt.ba"

let () = 
  let make_best infile = 
    let outall = infile ^ ".all" in
    let outbest = infile ^ ".best" in
    let predictor = pred_read infile in
    printf "Predicting%!";
    let t0 = Sys.time() in
    let preds = batch_predict_cat predictor test_data in
    printf "Done (%.2fs)\n%!" (Sys.time () -. t0);
    File.with_file_out outall (fun oc -> Array.iter (print_pred_conf oc) preds);
    let best = Array.enum preds |> best_predictions  in
    File.with_file_out outbest (fun oc -> output_preds oc best)
  in
  args () |> iter make_best
