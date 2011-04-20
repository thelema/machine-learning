open Batteries_uni
open Printf
open Bigarray
open Lacaml.Impl.S (* Single-precision reals *)

open Predict

let n = 200

let () =
  printf "Train_rows: %d\n" train_rows;
  match pred_read Sys.argv.(1) with
    | One_one (pairs, preds) as cp -> 
      let cs = group_by_cat n in
      let correct = ref 0 in
      let results = Array.create_matrix (category_count+1) (category_count+1) None in
      for idx = 0 to Array.length pairs - 1do
	let (i,j) = pairs.(idx) in
	let p = predict_b preds.(idx) in
	let iright = Array.fold_left (fun acc x -> 
	  if p (get_i train_data x) > 0. then acc+1 else acc)
	  0 cs.(i) in
	let jright = Array.fold_left (fun acc x -> if p (get_i train_data x) < 0. then acc+1 else acc) 0 cs.(j) in
	(*	printf "(%d,%d:%d,%d) " i j iright jright; *)
	if iright + jright > n + n/2 then (
	  printf "+%!"; 
	  incr correct;
	) else (
	  printf "x%!"; 
	  let data = Array.append cs.(i) cs.(j) in
	  let labels = Array.init (Array.length data) (fun x -> if x < Array.length cs.(i) then 1. else -1.) in
	  preds.(idx) <- klsq ~cutoff:0.3 ~lambda:1.0 ~sigma:1.0 data labels;
	);
	results.(i).(j) <- Some ((iright + jright) / (n/5));
      done;
      Array.print ~first:"" ~last:"" ~sep:"\n"
	(Array.print ~first:"" ~last:"" ~sep:"" (fun oc -> function None -> fprintf oc "." | Some d -> fprintf oc "%1d" d)) stdout results;
      printf "\n%d correct of %d\n" !correct (Array.length pairs);
      marshal_file (Filename.basename Sys.argv.(1)) cp (* write the predictor back to disk in its modified form *)
    | _ -> failwith "Not a one-one predictor"
