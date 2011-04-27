open Batteries_uni
open Printf
open Predict
open Bigarray

let test_data = Datafile.get_matrix "testing.txt.ba"


let () = 
  let make_best infile = 
    let outall = infile ^ ".all" in
    let predictor = pred_read infile in
    printf "Predicting%!";
    let t0 = Sys.time() in
    let preds = predict_many (predict_cat predictor) test_data |> List.enum in
    File.with_file_out outall (fun oc -> Enum.iter (print_pred_conf oc) preds);
    printf "Done (%.2fs)\n%!" (Sys.time () -. t0);
  in
  args () |> iter make_best
