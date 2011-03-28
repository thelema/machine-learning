open Batteries_uni
open Printf

open Bigarray
open Lacaml.Impl.D

let debug = false
let rows = 125000
let cols = 900

let cats = 2
let top_n = 100

let pred w x = let p = dot ~x w in abs_float p, p>0.
let right = ref 0
let count = ref 0

let perceptron () =
  let w = Vec.make0 cols in
  let wsum = Vec.make0 cols in
  let learn (xi, yi) =
    let _,pred = pred w xi in
    if pred <> yi then 
      if yi 
      then (if debug then printf "+%!"; ignore(Vec.add ~n:cols ~z:w xi w))
      else (if debug then printf "-%!"; ignore(Vec.sub ~n:cols ~z:w xi w))
    else incr right;
    incr count;
    ignore(Vec.add ~n:cols ~z:wsum wsum w);
  in
  (learn, fun () -> scal (1. /. float !count) wsum; wsum)

let to_vec l = List.map float_of_string l |> Vec.of_list

(*
let read_data_line line = String.nsplit line " " |> to_vec
let read_data_file fn = File.lines_of fn /@ read_data_line *)


let read_data_file fn = 
  let fh = Unix.openfile fn [Unix.O_RDONLY] 0o755 in
  let mmap = Array2.map_file fh int8_signed fortran_layout false cols rows in
  let get_row i = Array2.slice_right mmap i |> Array1.map float float64 in
  (1--rows) |> Enum.map get_row
    

let read_label_file fn = File.lines_of fn /@ int_of_string /@ (fun y -> y land 1 = 0)

let read_data data label = 
  Enum.combine (read_data_file data, read_label_file label)

let rec insert x = function
  | [] -> [x]
  | y::ys when y < x -> y::(insert x ys)
  | ys -> x::ys;;

let top n xs =
  let rec loop1 acc xs = function
    | 0 -> (acc, xs)
    | n -> loop1 (insert (List.hd xs) acc) (List.tl xs) (n - 1)
  in
  let rec loop2 acc = function
    | [] -> List.rev acc
    | x::xs -> loop2 (List.tl (insert x acc)) xs
  in
  let (acc, xs) = loop1 [] xs n in
  loop2 acc xs;;

let push (k,v) (map,min_k,count as _acc) = 
  if count < top_n then 
    (Map.add k v map, min k min_k, count+1) 
(*  else if k < min_k then acc *)
  else 
    let m = Map.add k v map |> Map.remove min_k in 
    (m, Map.min_binding m |> fst, count)

let test w acc (x, y) = 
  let str, pred = pred w x in 
  Map.modify_def (Map.create Float.compare, 0.,0) y (push (str, pred)) acc

let score_map y items = 
  let add_points pts pred acc = if pred = y then acc + pts else acc in
  let points = Enum.fold2 add_points 0 (top_n --- 1) (List.enum items) in
  float points /. float top_n

let train_file = "development.ba"
let train_label_file = "development_label.txt"
let test_file = "development.ba"
let test_label_file = "development_label.txt"

let () = 
  let data = read_data train_file train_label_file in
  let lrn, ret = perceptron () in
  printf "Training..%!";
  let t0 = Sys.time() in
  Enum.iter lrn data;
  printf "Done training(%.2fs)\n%!" (Sys.time () -. t0);
  printf "Correct during training: %d of %d\n" !right !count;
  let w = ret() in
(*  Vec.to_array w |> Array.print Float.print stdout; print_newline ();*)
  let data = read_data train_file train_label_file in
  let t0 = Sys.time() in
  let predict_maps = Enum.fold (test w) Map.empty data in
  let top_n_map = Map.map (fun (m,_,_) -> Map.enum m |> map snd |> List.of_backwards) predict_maps in
  let scores = Map.mapi score_map top_n_map in
  printf "Scores: %a (%.2fs)\n" (Map.print Bool.print Float.print) scores (Sys.time () -. t0);
  let overall = Map.enum scores |> map snd |> Enum.reduce (+.) in
  printf "Overall: %.2f\n" (overall /. float cats) ;
  ()


(*
let () = printf "Reading..%!"
let t0 = Sys.time()
let v1 = svm_scale ~lower:0. alldata
let m = svm_train ~kernel_type:LINEAR v1
let () = printf "Done(%.2f)\n%!" (Sys.time () -. t0)
 *)
 
