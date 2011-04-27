open Batteries_uni
open Printf

(* Read from many foo.all and produce a global.best with top guesses across all files *)

let top_n = 100
let keep = 250

let push (k:float) v (map,min_k,count as acc) = 
  if count < keep then 
    if Map.mem k map then acc 
    else (Map.add k v map, min k min_k, count+1)
  else if k < min_k then acc
  else 
    if Map.mem k map then acc else
    let m = Map.remove min_k map |> Map.add k v in 
    (m, Map.min_binding m |> fst, count)

let push k v a = 
  let (m,_,c as b) = push k v a in
   if (Map.enum m |> Enum.count <> c) then (
     printf "Count: %d Map keys: (%d) %a" c (Map.enum m |> Enum.count) (Enum.print Float.print) (Map.keys m);
     assert false;
   ) else b


let push_pred acc (id, i, pred, str) = 
  Map.modify_def (Map.create Float.compare, 1.,0) pred (push str (id, i, str)) acc

let best_predictions preds = 
  (* highest to lowest prediction strength *)
  let to_list (m,k,c) = 
    printf "threshold:%f  count: %d  map_entries:%d\n" k c (Map.enum m |> Enum.count);
    Map.enum m |> Enum.map snd |> List.of_backwards 
    |> List.unique ~cmp:(fun (_,(i1:int),_) (_,i2,_) -> i1 = i2) 
    |> List.take top_n 
  in
  Enum.fold push_pred Map.empty preds |> Map.map to_list

let print_pred oc (label,ps) = 
(*  List.iter (fun (id, i, str) -> fprintf oc "%d\t%d\t#%.3f from %s\n" label i str id) ps *)
  List.iter (fun (_id, i, _str) -> fprintf oc "%d\t%d\n" label i) ps 

let output_preds oc ps =
  Enum.print ~first:"" ~last:"" ~sep:"" print_pred oc (Map.enum ps)
    
let parse_pred fn i str = 
  let label, strength = String.split str "\t" in
  (fn, i+1, int_of_string label, float_of_string strength)
    
let get_preds infile = File.lines_of infile |> Enum.mapi (parse_pred infile)

let outbest = "global.best"

let () = 
  let best = args () |> map get_preds |> Enum.flatten |> best_predictions in
  File.with_file_out outbest (fun oc -> output_preds oc best)
