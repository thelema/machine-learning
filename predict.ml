open Batteries_uni
open Printf

open Libosvm
open Lacaml.Impl.S (* Single-precision reals *)
open Ocamlviz
open Bigarray

let print_float oc x = fprintf oc "%.2f" x
let print_float5 oc x = fprintf oc "%.5f" x

let real_run = Array.length Sys.argv > 1
let debug = false
let cat_bits = ref 16
let category_count = 146
let top_n = 10
let loops = 10

let train_file = "training.ba"
let train_label_file = "training_label.txt"
let test_file = if real_run then Sys.argv.(1) else "development.ba"
let test_label_file = "development_label.txt"


module type Private_int = sig type t val of_int : int -> t val to_int : t -> int val print : 'a IO.output -> t -> unit val compare : t -> t -> int end
module Bit_cat : Private_int = Int
module Cat : Private_int = Int

let pred_split p = abs_float p, (if p=0. then Random.bool () else p>0.)

let norm ais = Array.fold_left (fun acc a -> a *. a +. acc) 0. ais |> sqrt
let zeros ais = Array.fold_left (fun acc a -> if a = 0. then acc+1 else acc) 0 ais


(********************************************)
(** Two-category learners *******************)
(********************************************)
let right = ref 0
let count = ref 0

(*
let perceptron acc =
  let w = Vec.make0 Datafile.cols in
  fun xi yi -> 
    if debug then printf "XI(%d):\n%a\n%!" !count (Array.print print_float) (Vec.to_array xi);
    if dot xi w *. yi <= 0. then axpy ~alpha:yi ~x:xi w else incr right;
    incr count;
    ignore(Vec.add ~z:acc acc w);
    let n2 = Vec.sqr_nrm2 w in
(*    if !count land 0xff = 0 then printf "i:%d D:%.2f, dec:%B really:%B wnorm:%.2f\n" !count d y yi n2; *)
    if n2 > 1000. then (printf "s"; scal (500. /. n2) w;)
 *)

let per2 data labels =
  count := 0; right := 0;
  let accum w (xi,yi) =
    incr count;
    if dot xi w *. yi > 0. then (incr right; w)
    else if yi > 0. then Vec.add w xi else Vec.sub w xi
  in
  let n = Array2.dim1 data in
  let data_label = (0--^n) |> map (fun i -> Array2.slice_right data i, labels.{i}) in
  let ws = Enum.scanl accum (Vec.make0 Datafile.cols) data_label in
  let wsum = Vec.make0 Datafile.cols in
  Enum.iter (fun wi -> ignore(Vec.add ~z:wsum wsum wi)) ws;
  scal (1. /. float n) wsum;
  (fun x -> dot x wsum)

let per3 data labels =
  let w = Vec.make0 Datafile.cols in
  let acc = Vec.make0 Datafile.cols in
  for i = 1 to Array2.dim1 data do
    incr count;
    let xi = Array2.slice_right data i in
    if dot xi w *. yi > 0. then incr right
    else axpy ~alpha:yi ~x:xi w;
    ignore(Vec.add ~z:acc acc w);
    let n2 = Vec.sqr_nrm2 w in
    if n2 > 1000. then (scal (800. /. n2) w;)
  done;
  scal (1. /. float n) acc;
  acc

let kp_core n kij labels ais =
  let rec sum i acc j = 
    if j > n then acc 
    else sum i (ais.(j) *. kij.{i,j} +. acc) (j+1) 
  in
  for i = 1 to n do
    if labels.{i} *. (sum i 0. 1) <= 0. then 
      ais.{i} <- ais.{i} +. yi
  done
 
let kperceptron k data labels =
  let n = Array1.dim labels in
  let ais = Vec.make0 n in
  let kij = Array2.create float32 fortran_layout n n in
  printf "Generating kij...%!";
  for i = 1 to n do
    let xi = (Array2.slice_right data i) in
    for j = 1 to n do
      kij.{i,j} <- k xi (Array2.slice_right data j)
    done;
  done;
  printf "\nLooping through predictions...%!";
  for l = 1 to 10 * loops do kp_core n kij labels ais done;
  printf "Done. (ai_norm,max,zeros = %.2f,%.0f,%d)\n%!" (norm ais) (Array.reduce max ais) (zeros ais);
  (fun x -> let t = ref 0. in for i = 0 to n-1 do t := !t +. ais.(i) *. k x (Array2.slice_right data i); done; !t +. 0.)

(*
let kp_core_nocache n k data labels ais =
  let rec sum i acc j = 
    if j >= n then acc 
    else sum i (ais.(j) *. k data.(i) data.(j) +. acc) (j+1) 
  in
  Array.iteri (fun i yi -> if yi *. (sum i 0. 0) <= 0. then ais.(i) <- ais.(i) +. yi) labels

let kperceptron_nocache k data labels =
  let n = Array.length data in
  let ais = Array.create n 0. in
  for l = 1 to loops do kp_core_nocache n k data labels ais done;
  (fun x -> let t = ref 0. in for i = 0 to n-1 do t := !t +. ais.(i) *. k x data.(i); done; !t +. 0.)
  *)

let k_2 x1 x2 = let a = (1. +. dot x1 x2) in a *. a
let k_3 x1 x2 = let a = (1. +. dot x1 x2) in a *. a *. a
let k_rbf two_s_sqr x1 x2 = let a = Vec.sqr_nrm2 (Vec.sub x1 x2) in exp (~-. a /. two_s_sqr)

(*
let () = printf "Reading..%!"
let t0 = Sys.time()
let v1 = svm_scale ~lower:0. alldata
let m = svm_train ~kernel_type:LINEAR v1
let () = printf "Done(%.2f)\n%!" (Sys.time () -. t0)
 *)

let svm data labels = 
  let _v1 = svm_train ~kernel_type:LINEAR data in
  assert false

(********************************************)
(** EXTENDING 2-category to multi-category **)
(********************************************)

let mask = (1 lsl !cat_bits) - 1

type ('a,'b) manual_cache = {
  get : 'a -> 'b; 
  del : 'a -> unit; 
  enum: unit -> ('a * 'b) BatEnum.t
}

let make_map ~gen =
  let m = ref BatMap.empty in
  {get = (fun k -> 
    try BatMap.find k !m
    with Not_found -> gen k |> tap (fun v -> m := BatMap.add k v !m));
   del = (fun k -> m := BatMap.remove k !m);
   enum = (fun () -> BatMap.enum !m) }


let {get = map_cat; enum=cat_mapping} = 
  make_map (fun _ -> Bit_cat.of_int (Random.bits () land mask))

let push ((k:float),v) (map,min_k,count as acc) = 
  if count < top_n then 
    (Map.add k v map, min k min_k, count+1) 
  else if k < min_k then acc 
  else 
    let m = Map.remove min_k map |> Map.add k v in 
    (m, Map.min_binding m |> fst, count)

let rec popcount c x = if x = 0 then c else popcount (c+1) (x land (x-1))
let ham_dist x y = popcount 0 (x lxor Bit_cat.to_int y) (*|> tap (fun d -> printf "HD(%x,%x) = %d " x (Bit_cat.to_int y) d) *)

let nearest_in label_map pred_lab = 
    Enum.arg_min (fun (_,bit_lab) -> ham_dist pred_lab bit_lab) (Array.enum label_map) |> fst

let rec to_bits n x = if n = 0 then [] else (if x land 1 = 1 then 1. else -1.) :: (to_bits (n-1) (x asr 1))
let merge_qb (q,p) pred = (q *. abs_float pred, if pred >= 0. then p lsl 1 + 1 else p lsl 1)

type 'a pred_f_t = vec -> float * 'a

let extend_hamm_incr incr_learner data =
  let ws = List.init !cat_bits (fun _ -> Vec.make0 Datafile.cols) in
  let ps = List.map incr_learner ws in
  let lrn_data () = 
    Enum.iter (fun (xi,yi) -> 
      let bits = to_bits !cat_bits (map_cat yi |> Bit_cat.to_int) in
      List.iter2 (fun p b -> p xi b) ps bits) (data ())
  in
  printf "Training hamm_incr..\n%!";
  let t0 = Sys.time() in
  for i = 1 to loops do
    let err_in = !count - !right in
    lrn_data ();
    printf "Loop: %d, err: %d, nrm2s: %a\n%!" i (!count - !right - err_in) (List.print print_float5) (List.map (fun w -> Vec.sqr_nrm2 w /. float !count /. float !count) ws); 
    
  done;
  printf "Done training(%.2fs)\n%!" (Sys.time () -. t0);
  printf "Correct during training: %d of %d\n" !right !count;
  let scal_f = 1. /. float !count in
  List.iter (fun w -> scal scal_f w) ws;
  let cat_map = cat_mapping () |> Array.of_enum in
  let decode_arr = Array.init (1 lsl !cat_bits) (fun i -> nearest_in cat_map i) in
  printf "CatMap: %a\n" (Array.print (Pair.print Cat.print Bit_cat.print)) cat_map;
  let pred_n x = 
    List.rev_map (fun w -> dot w x) ws 
    |> List.fold_left merge_qb (1.,0) 
    |> second (fun i -> decode_arr.(i))
  in
  (pred_n : Cat.t pred_f_t)


let train_t = Time.create "train"
let eh_pred_p = Point.create "eh_pred"

let extend_hamm gen_classifier data labels =
  Time.start train_t;
  let mask = (1 lsl !cat_bits) - 1 in
  let {get = map_cat; enum=cat_mapping} = 
    make_map (fun _ -> Random.bits () land mask) in
  let bit_labels i x = if ((map_cat x) asr i) land 1 = 1 then 1. else -1. in
  printf "Training Classifiers...\n%!";
  let classifiers = List.init !cat_bits (fun i -> gen_classifier data (Array.map (bit_labels i) labels)) in
  printf "Done training\n%!";
  let cat_map = cat_mapping () |> Enum.map (second Bit_cat.of_int) |> Array.of_enum in
  let decode_arr = Array.init (1 lsl !cat_bits) (fun i -> nearest_in cat_map i) in
  Time.stop train_t;
  let pred_n x = 
    Point.observe eh_pred_p;
    List.rev_map (fun classifier -> classifier x) classifiers
    |> List.fold_left merge_qb (1.,0) 
    |> (fun (i,j) -> i, decode_arr.(j))
  in
  (pred_n : Cat.t pred_f_t)


let extend_one_one_incr incr_learner data =
  right := 0; count := 0;
  printf "Training 1-1..\n%!";
  let t0 = Sys.time() in
  let mm = (1 --^ Array2.dim1 data) |> Enum.fold (fun acc i -> Map.modify_def [] (Cat.to_int y) (List.cons x) acc) (Map.create Int.compare) (data_f ()) in
  let pairs = Map.keys mm |> map (fun i -> Map.keys mm // (fun x -> x > i) /@ (fun j -> i,j)) |> Enum.flatten |> List.of_enum in
  let idata i y = try Map.find i mm |> List.enum |> Enum.map (fun x -> (x,y)) with Not_found -> failwith ("No values of category " ^ string_of_int i) in
  let train_p (i,j) = 
    let err_in = !count - !right in
    let ws = Vec.make0 Datafile.cols in
    let pc_learn = incr_learner ws in
    Enum.append (idata i 1.) (idata j (-1.)) 
    |> Random.shuffle |> Array.enum
    |> iter (fun (i,j) -> pc_learn i j);
    printf "1-1: (%d,%d), err: %d, nrm2s: %a\n%!" i j (!count - !right - err_in) print_float5 (Vec.sqr_nrm2 ws /. float (!count * !count)); 
    ws
  in
  let ws = List.map train_p pairs in
  printf "Done training(%.2fs)\n%!" (Sys.time () -. t0);
  printf "Correct during training: %d of %d\n" !right !count;
  let max_cat = Enum.reduce max (Map.keys mm) in
  let cat_count = Map.keys mm |> Enum.count in
  let pred_n x =
    let votes = Array.create (max_cat+1) 0 in
    let pc_vote (i,j) w = 
      let vote = if dot w x > 0. then i else j in
      votes.(vote) <- votes.(vote) + 1
    in
    List.iter2 pc_vote pairs ws;
    let winner = Enum.arg_max (fun i -> votes.(i)) (1--category_count) in
    float votes.(winner) /. float cat_count, Cat.of_int winner
  in
  (pred_n: Cat.t pred_f_t)

(*
let extend_one_one gen_classifier data labels =
  right := 0; count := 0;
  printf "Training 1-1..\n%!";
  let t0 = Sys.time() in
  let mm = Enum.fold (fun acc (x,y) -> Map.modify_def [] (Cat.to_int y) (List.cons x) acc) (Map.create Int.compare) (data_f ()) in
  let pairs = Map.keys mm |> map (fun i -> Map.keys mm // (fun x -> x > i) /@ (fun j -> i,j)) |> Enum.flatten |> List.of_enum in
  let idata i y = try Map.find i mm |> List.enum |> Enum.map (fun x -> (x,y)) with Not_found -> failwith ("No values of category " ^ string_of_int i) in
  let train_p (i,j) = 
    let err_in = !count - !right in
    let predictor = Enum.append (idata i true) (idata j false) 
		   |> Random.shuffle |> Array.enum |> gen_classifier in
    printf "1-1: (%d,%d), err: %d, nrm2s: %a\n%!" i j (!count - !right - err_in) print_float5 (Vec.sqr_nrm2 ws /. float (!count * !count));
    predictor
  in
  let ps = List.map train_p pairs in
  printf "Done training(%.2fs)\n%!" (Sys.time () -. t0);
  printf "Correct during training: %d of %d\n" !right !count;
  let max_cat = Enum.reduce max (Map.keys mm) in
  let cat_count = Map.keys mm |> Enum.count in
  let pred_n x =
    let votes = Array.create (max_cat+1) 0 in
    let pc_vote (i,j) p = 
      let vote = if p x > 0. then i else j in
      votes.(vote) <- votes.(vote) + 1
    in
    List.iter2 pc_vote pairs ps;
    let winner = Enum.arg_max (fun i -> votes.(i)) (1--category_count) in
    float votes.(winner) /. float cat_count, Cat.of_int winner
  in
  (pred_n: Cat.t pred_f_t)
*)

(***************************************)
(** Scoring and framework **************)
(***************************************)


let score_map pred items = 
  let len = List.length items in
  let add_points pts (i,y) acc = if pred = y then acc + pts else acc in
  let points = Enum.fold2 add_points 0 (len --- 1) (List.enum items) in
  float points /. float len

let score_map is_correct guess_cat guess_ids = 
  let correct_enum = Enum.map (is_correct guess_cat) (List.enum guess_ids) in
  let is_correct _ = match Enum.get correct_enum with Some true -> true | _ -> false in
  let points = Enum.filter is_correct (top_n --- 1) |> Enum.fold (+) 0 in
  float points /. float top_n

let push_pred i (str, pred) acc = 
  Map.modify_def (Map.create Float.compare, 0.,0) pred (push (str, i)) acc

let best_predictions preds = 
  (* highest to lowest prediction strength *)
  let to_list (m,_,_) = Map.enum m |> Enum.map snd |> List.of_backwards in
  Enum.foldi push_pred Map.empty preds |> Map.map to_list

let read_label_file fn = File.lines_of fn /@ int_of_string |> Array.of_enum |> Array1.of_array int8_unsigned Datafile.layout

let train_full pred_gen =
  pred_gen (Datafile.get_matrix train_file) (read_label_file train_label_file)
    
let train_slice ?(skip=0) n pred_gen =
  let data = Array2.sub_right (Datafile.get_matrix train_file) (skip+1) n in
  let labels = Array1.sub (read_label_file train_label_file) (skip+1) n in
  pred_gen data labels

let train_slices n pred_gen =
  let fulldata = Datafile.get_matrix train_file in
  let slices = Array2.dim1 fulldata / n in
  let datas = List.init slices (fun i -> Array2.sub_right fulldata (1+n*i) n) in
  let fulllabels = read_label_file train_label_file in
  let labels = List.init slices (fun i -> Array1.sub fulllabels (1+n*i) n) in
  List.combine datas labels

let evaluate best = 
(*  printf "Predictions:\n%a\n" (Map.print Cat.print (List.print Int.print)) best; *)
  let labels = read_label_file test_label_file in
  let is_correct guess i = try labels.{i-1} = guess with Invalid_argument _ -> printf "Item %d out of range\n" i; false in
  let scores = Map.mapi (score_map is_correct) best in
  printf "Scores: %a\n" (Map.print ~sep:", " Int.print print_float) scores;
  let overall = Map.enum scores |> map snd |> Enum.reduce (+.) in
  overall /. float category_count

let output_preds ps =
  let print_pred oc (i,ps) = 
    List.iter (fun p -> fprintf oc "%d\t%d\n" i p) ps
  in
  Enum.print ~first:"" ~last:"" ~sep:"" print_pred stdout (Map.enum ps)

let pred_t = Time.create "pred"
let pred_p = Point.create "pred"

let run_test ?(n = 100_000) ?(print=false) name pred = 
  Time.start pred_t; Point.observe pred_p;
  let data = Datafile.get_matrix test_file in
  (0--^n) |> Enum.map (fun i -> pred (Array2.slice_right data i))
  |> best_predictions 
  |> tap (fun _ -> Time.stop pred_t)
  |> evaluate |> printf "%s Overall: %.2f\n%!" name
 
(***************************************)
(**        MAIN        *****************)
(***************************************)

let () = 
    (*  check_data (); *)
(*
  extend_hamm (kperceptron k_3) |> train_slice 2000 |> 
      run_test ~print:real_run "KP3H";
*)

  extend_hamm (kperceptron k_3) |> train_slices 1000 |>
      Enum.iter (run_test ~print:real_run "KP3*" ~n:10000);

  extend_hamm_incr perceptron |> train_full |>
      run_test ~print:real_run "Hamm";

  extend_one_one_incr perceptron |> train_full |>
      run_test ~print:real_run "1-1 ";
  ()
