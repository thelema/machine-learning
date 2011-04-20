open Batteries_uni
open Printf
open Bigarray
open Lacaml.Impl.S (* Single-precision reals *)

open Predict

let n = 100

let () =
  printf "Train_rows: %d\n" train_rows;
  match pred_read Sys.argv.(1) with
    | One_one (pairs, preds) -> 
      let cs = group_by_cat n in
      let correct = ref 0 in
      let results = Array.create_matrix (category_count+1) (category_count+1) None in
      let test_pred (i,j) p = 
	let p = predict_b p in
	let iright = Array.fold_left (fun acc x -> 
	  if p (get_i train_data x) > 0. then acc+1 else acc)
	  0 cs.(i) in
	let jright = Array.fold_left (fun acc x -> if p (get_i train_data x) < 0. then acc+1 else acc) 0 cs.(j) in
(*	printf "(%d,%d:%d,%d) " i j iright jright; *)
	if iright + jright > n + n/2 then (printf "+%!"; incr correct;) else printf "x%!";
	results.(i).(j) <- Some ((iright + jright) / (n/5))
      in
      Array.iter2 test_pred pairs preds;
      Array.print ~first:"" ~last:"" ~sep:"\n"
	(Array.print ~first:"" ~last:"" ~sep:"" (fun oc -> function None -> fprintf oc "." | Some d -> fprintf oc "%1d" d)) stdout results;
      printf "\n%d correct of %d\n" !correct (Array.length pairs)
    | _ -> failwith "Not a one-one predictor"
