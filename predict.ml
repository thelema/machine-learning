open Batteries_uni
open Printf

open Libosvm
open Lacaml.Impl.S (* Single-precision reals *)
open Ocamlviz
open Bigarray

type matrix = (float, float32_elt, fortran_layout) Array2.t
type matrix64 = (float, float64_elt, fortran_layout) Array2.t
type cats = (int, int8_unsigned_elt, fortran_layout) Array1.t
type 'a pred_f_t = matrix -> (float * 'a) Enum.t

let print_float oc x = fprintf oc "%.2f" x
let print_float5 oc x = fprintf oc "%.5f" x

let real_run = Array.length Sys.argv > 1
let debug = false
let cat_bits = ref 16
let category_count = 146
let top_n = 10
let loops = 10

let train_data = Datafile.get_matrix "training.ba"
let train_labels = Datafile.read_label_file "training_label.txt"

let test_data = Datafile.get_matrix (if real_run then Sys.argv.(1) else "development.ba")
let test_labels = Datafile.read_label_file "development_label.txt"

let pred_split p = abs_float p, (if p=0. then Random.bool () else p>0.)

let norm ais = Vec.sqr_nrm2 ais |> sqrt
let zeros ais = (1--Array1.dim ais) |> Enum.fold (fun acc i -> if ais.{i} = 0. then acc+1 else acc) 0
let get_i m i = Array2.slice_right m i

type bpredictor = 
  | Dot of float array (* weights of different features *)
  | Kern_3 of int * int * float array (* offset/len of data used, ais*)
  | Kern_rbf of float * int * int * float array (* sigma, offset/len of data used, ais*)
  | Svm of string (* filename with serialized predictor *)

type cpredictor = 
  | Hamm of bpredictor array * int array (* n predictors, n-bit codewords *)

(********************************************)
(** Two-category learners *******************)
(********************************************)
let right = ref 0
let count = ref 0

let perceptron (data: matrix) (labels: vec) =
  let w = Vec.make0 Datafile.cols in
  let acc = Vec.make0 Datafile.cols in
  let n = Array1.dim labels in
  for i = 1 to n do
    let xi = Array2.slice_right data i in
    let yi = labels.{i} in
    if dot xi w *. yi > 0. then incr right
    else axpy ~alpha:yi ~x:xi w;
    ignore(Vec.add ~z:acc acc w);
    let n2 = Vec.sqr_nrm2 w in
    if n2 > 1000. then (scal (800. /. n2) w;)
  done;
  count := !count + n;
  scal (1. /. float n) acc;
  Dot (Array1.to_array acc)

let kp_core n (kij: matrix) (labels: vec) (ais: vec) =
  let rec sum i acc j = 
    if j > n then acc 
    else sum i (ais.{j} *. kij.{i,j} +. acc) (j+1) 
  in
  for i = 1 to n do
    let yi = labels.{i} in
    if yi *. (sum i 0. 1) <= 0. then 
      ais.{i} <- ais.{i} +. yi
  done
 
let kperceptron kij (labels:vec) =
  let n = Array1.dim labels in
  let ais = Vec.make0 n in
  printf "\nLooping through %d training data items %d times...%!" n (10*loops);
  for l = 1 to 10 * loops do kp_core n kij labels ais done;
  printf "Done. (ai_norm,zeros = %.2f,%d)\n%!" (norm ais) (zeros ais);
  ais

let k_rbf two_s_sqr x1 x2 = let a = Vec.sqr_nrm2 (Vec.sub x1 x2) in exp (~-. a /. two_s_sqr)
let gen_kij_rbf two_sig_sq data =
  let n = Array2.dim2 data in
  let kij = Array2.create float32 fortran_layout n n in
  printf "Generating kij...%!";
  for i = 1 to n do
    let xi = (Array2.slice_right data i) in
    for j = 1 to n do
      kij.{i,j} <- k_rbf two_sig_sq xi (Array2.slice_right data j)
    done;
  done;
  kij

let k_3 x1 x2 = let a = (1. +. dot x1 x2) in a *. a *. a
let gen_kij_3 data =
  let n = Array2.dim2 data in
  let kij = Array2.create float32 fortran_layout n n in
  printf "Generating kij...%!";
  for i = 1 to n do
    let xi = (Array2.slice_right data i) in
    for j = 1 to n do
      kij.{i,j} <- k_3 xi (Array2.slice_right data j)
    done;
  done;
  kij


let kperceptron3_slice off len labels_arr =
  let kij = gen_kij_3 (Array2.sub_right train_data off len) in
  let tag a = Kern_3 (off, len, Array1.to_array a) in
  Array.map (kperceptron kij |- tag) labels_arr

let kperceptron_rbf_slice sigma off len labels_arr =
  let kij = gen_kij_rbf sigma (Array2.sub_right train_data off len) in
  let tag a = Kern_rbf (sigma, off, len, Array1.to_array a) in
  Array.map (kperceptron kij |- tag) labels_arr

(*
let () = printf "Reading..%!"
let t0 = Sys.time()
let v1 = svm_scale ~lower:0. alldata
let m = svm_train ~kernel_type:LINEAR v1
let () = printf "Done(%.2f)\n%!" (Sys.time () -. t0)
 *)

let svm (data: matrix64) (labels: vec) = 
  let _model = svm_train ~svm_type:C_SVC ~kernel_type:RBF data (* labels *) in
  assert false

(** PREDICTION USING A BPREDICTOR **)

let predict_many d f = 
  (1--Array2.dim2 d) |> Enum.map (fun i -> f (get_i d i))

let predict_b = 
  function
  | Dot warr -> let w = Array1.of_array Datafile.kind Datafile.layout warr in (fun x -> dot w x)
  | Kern_3 (off, len, ais) ->
    let data = Array2.sub_right train_data off len in
    let ais = Array1.of_array Datafile.kind Datafile.layout ais in
    (fun x -> 
      let t = ref 0. in 
      for i = 1 to len do 
	t := !t +. ais.{i} *. k_3 x (Array2.slice_right data i); 
      done; 
      !t +. 0.)
  | Kern_rbf (sigma, off, len, ais) -> 
    let data = Array2.sub_right train_data off len in
    let ais = Array1.of_array Datafile.kind Datafile.layout ais in
    (fun x -> 
      let t = ref 0. in 
      for i = 1 to len do 
	t := !t +. ais.{i} *. k_rbf sigma x (Array2.slice_right data i); 
      done; 
      !t +. 0.)
  | Svm _fn -> assert false

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

let push ((k:float),v) (map,min_k,count as acc) = 
  if count < top_n then 
    (Map.add k v map, min k min_k, count+1) 
  else if k < min_k then acc 
  else 
    let m = Map.remove min_k map |> Map.add k v in 
    (m, Map.min_binding m |> fst, count)

let rec popcount c x = if x = 0 then c else popcount (c+1) (x land (x-1))
let ham_dist x y = popcount 0 (x lxor y) (*|> tap (fun d -> printf "HD(%x,%x) = %d " x (Bit_cat.to_int y) d) *)

let nearest_in label_map pred_lab = 
    (Array.range label_map) |> Enum.arg_min (fun i -> ham_dist pred_lab label_map.(i))

let rec to_bits n x = if n = 0 then [] else (if x land 1 = 1 then 1. else -1.) :: (to_bits (n-1) (x asr 1))
let merge_qb (q,p) pred = (q *. abs_float pred, if pred >= 0. then p lsl 1 + 1 else p lsl 1)

let train_t = Time.create "train"
let eh_pred_p = Point.create "eh_pred"

let extend_hamm gen_classifier =
  Time.start train_t;
  let n = Array1.dim train_labels in
  let mask = (1 lsl !cat_bits) - 1 in
  let {get = map_cat; enum=cat_mapping} = 
    make_map (fun _ -> Random.bits () land mask) in
  let bit_labels i x = if ((map_cat x) asr i) land 1 = 1 then 1. else -1. in
  let gen_labels i = 
    let ls = Vec.make0 n in 
    for j = 1 to n do ls.{j} <- bit_labels i train_labels.{j} done; 
    ls 
  in
  let all_labels = Array.init !cat_bits gen_labels in
  printf "Training Classifiers...\n%!";
  let classifiers = gen_classifier all_labels in
  printf "Done training\n%!";
  let cat_map = Array.create category_count 0 in
  Enum.iter (fun (i,b) -> cat_map.(i) <- b) (cat_mapping ());
  Time.stop train_t;
  Hamm (classifiers, cat_map)

let get_heads enum_l = 
  Array.fold_left (fun acc x -> match Enum.get x with 
      None -> raise Enum.No_more_elements 
    | Some v -> v :: acc) 
    [] enum_l

let predict_cat = function
  | Hamm (classifiers, cat_map) ->
    let decode_arr = Array.init (1 lsl !cat_bits) (fun i -> nearest_in cat_map i) in
    let pred_n xs = 
      Point.observe eh_pred_p;
      let decisions = Array.map (fun classifier -> predict_many xs (predict_b classifier)) classifiers in
      let decode l = 
	let (str, bits) = List.fold_left merge_qb (1.,0) l in
	str, decode_arr.(bits)
      in
      Enum.from (fun () -> get_heads decisions) |> Enum.map decode
    in
    pred_n

(* TODO
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

let train_full pred_gen = pred_gen train_data train_labels
    
let train_slice ?(skip=0) n pred_gen = pred_gen skip n

let train_slices n pred_gen =
  let slices = Array1.dim train_labels / n in
  List.init slices (fun i -> pred_gen (1+n*i) n)

let evaluate best = 
(*  printf "Predictions:\n%a\n" (Map.print Cat.print (List.print Int.print)) best; *)
  let is_correct guess i = try train_labels.{i+1} = guess with Invalid_argument _ -> printf "Item %d out of range\n" i; false in
  let scores = Map.mapi (score_map is_correct) best in
  printf "Scores: %a\n" (Map.print ~sep:", " Int.print print_float) scores;
  let overall = Map.enum scores |> map snd |> Enum.reduce (+.) in
  overall /. float category_count

let output_preds ps =
  let print_pred oc (i,ps) = 
    List.iter (fun p -> fprintf oc "%d\t%d\n" i p) ps
  in
  Enum.print ~first:"" ~last:"" ~sep:"" print_pred stdout (Map.enum ps)

let run_test ?(n = 100_000) name (pred: int pred_f_t) = 
  let d = if real_run then test_data else Array2.sub_right test_data 0 n in
  pred d |> best_predictions |> evaluate |> printf "%s Overall: %.2f\n%!" name
 
let marshal_file fn_base x =
  let tm = Unix.time() |> Unix.localtime in
  let fn = sprintf "%s.%d.%d.%d" fn_base tm.Unix.tm_mday tm.Unix.tm_hour tm.Unix.tm_min in
  let oc = Pervasives.open_out fn in
  Legacy.Marshal.to_channel oc x [];
  Pervasives.close_out oc

(***************************************)
(**        MAIN        *****************)
(***************************************)

let train () = 
    (*  check_data (); *)

  extend_hamm (kperceptron3_slice 0 2000) |> marshal_file "hamm_kp3_2k";
(*
      run_test ~print:real_run "KP3H";


  extend_hamm (kperceptron k_3) |> train_slices 1000 |>
      List.iter (run_test ~print:real_run "KP3*" ~n:10000);
*)
(*
  extend_hamm_incr perceptron |> train_full |>
      run_test ~print:real_run "Hamm";

  extend_one_one_incr perceptron |> train_full |>
      run_test ~print:real_run "1-1 ";
 *)

  ()

let predict () = ()


let () = if exe = "predict" then predict () else train ()
