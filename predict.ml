open Batteries_uni
open Printf

open Libosvm
open Lacaml.Impl.S (* Single-precision reals *)
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
let category_count = 164
let top_n = 100
let loops = 10

let train_data = Datafile.get_matrix "training.ba"
let train_labels = Datafile.read_label_file "training_label.txt"

let test_data = Datafile.get_matrix (if real_run then Sys.argv.(1) else "development.ba")
let test_labels = Datafile.read_label_file "development_label.txt"

let pred_split p = abs_float p, (if p=0. then Random.bool () else p>0.)

let norm ais = Vec.sqr_nrm2 ais |> sqrt
let zeros ais = (1--Array1.dim ais) |> Enum.fold (fun acc i -> if ais.{i} = 0. then acc+1 else acc) 0
let get_i (m:matrix) i = Array2.slice_right m i
let predict_many d f = (1--Array2.dim2 d) |> Enum.map (fun i -> f (get_i d i))


type bpredictor = 
  | Dot of float array (* weights of different features *)
  | Kern_3 of int array * float array (* offsets of data used, ais*)
  | Kern_rbf of float * int array * float array (* sigma, offsets of data used, ais*)
  | Svm of string (* filename with serialized predictor *)

type cpredictor = 
  | Hamm of bpredictor array * int array (* n predictors, n-bit codewords *)
  | One_one of (int * int) list * bpredictor array

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

let rec pred (kij: matrix) (ais: vec) i acc j = 
  if j < 1 then acc else 
    pred kij ais i (ais.{j} *. kij.{j,i} +. acc) (j-1) 

let kp_core n kij (labels: vec) (ais: vec) =
  for i = 1 to n do
    let yi = labels.{i} in
    if yi *. (pred kij ais i 0. n) <= 0. then 
      ais.{i} <- ais.{i} +. yi
  done
 
let kperceptron n kij labels =
  let ais = Vec.make0 n in
(*  printf "Training: %d data items %d times...%!" n (10*loops); 
  let t0 = Sys.time () in*)
  for l = 1 to 10 * loops do kp_core n kij labels ais done;
(*  printf "Done(%.2fs). (ai_norm,zeros = %.2f,%d)\n%!" (Sys.time () -. t0) (norm ais) (zeros ais); *)
  ais

let k_rbf two_s_sqr x1 x2 = let a = Vec.sqr_nrm2 (Vec.sub x1 x2) in exp (~-. a /. two_s_sqr)
let gen_kij_rbf two_sig_sq n offsets =
  let kij = Array2.create float32 fortran_layout n n in
(*  printf "Generating kij(%dx%d) using rbf(x,y)...%!" n n; *)
  for i = 1 to n do
    let xi = (Array2.slice_right train_data offsets.(i)) in
    for j = 1 to n do
      kij.{i,j} <- k_rbf two_sig_sq xi (Array2.slice_right train_data offsets.(j))
    done;
  done;
(*  printf "Done\n%!"; *)
  kij

let k_3 x1 x2 = let a = (1. +. dot x1 x2) in a *. a *. a
let gen_kij_3 n offsets =
  let kij = Array2.create float32 fortran_layout n n in
(*  printf "Generating kij(%dx%d) using (1 + x dot y)^3...%!" n n;*)
  for i = 1 to n do
    let xi = (Array2.slice_right train_data offsets.(i)) in
    for j = 1 to n do
      kij.{j,i} <- k_3 xi (Array2.slice_right train_data offsets.(j))
    done;
  done;
(*  printf "Done\n%!"; *)
  kij

let k_pow pow x1 x2 = let a = (1. +. dot x1 x2) in a ** pow
let gen_kij_pow pow n offsets =
  let kij = Array2.create float32 fortran_layout n n in
(*  printf "Generating kij(%dx%d) using (1 + x dot y)^3...%!" n n;*)
  for i = 1 to n do
    let xi = (Array2.slice_right train_data offsets.(i)) in
    for j = 1 to n do
      kij.{j,i} <- k_3 xi (Array2.slice_right train_data offsets.(j))
    done;
  done;
(*  printf "Done\n%!"; *)
  kij


let clean a offs =
  let offs = Array.filteri (fun i _ -> a.{i} <> 0.) offs in
  let a = Array1.to_array a |> Array.filter ((<>) 0.) in
  (a, offs)
  

let kperceptron_offs genk len offs labels_arr =
  let kij = genk len offs in
  let tag a = let a, offs = clean a offs in Kern_3 (offs, a) in
  Array.map (kperceptron len kij |- tag) labels_arr

let kperceptron3_offs = kperceptron_offs gen_kij_3
let kperceptron3_offs sigma = kperceptron_offs (gen_kij_rbf sigma)

let kperceptron_slice genk off len labels_arr =
  let offs = (off --^ (off + len) |> Array.of_enum) in 
  kperceptron_offs genk len offs labels_arr

let kperceptron3_slice = kperceptron_slice gen_kij_3
let kperceptron_rbf_slice sigma = kperceptron_slice (gen_kij_rbf sigma)
let kperceptron_pow_slice pow = kperceptron_slice (gen_kij_pow pow)


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

let extend_hamm gen_classifier =
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
  let cat_map = Array.create (category_count+1) 0 in
  Enum.iter (fun (i,b) -> cat_map.(i) <- b) (cat_mapping ());
  printf "Training Hamming Classifiers...\n%!";
  let t0 = Sys.time () in
  let classifiers = gen_classifier all_labels in
  printf "Done training (%.2fs)\n%!" (Sys.time () -. t0);
  Hamm (classifiers, cat_map)

(* binary search, returning +1. on found and -1. on not found *)
let rec search a ?(l=0) ?(u=Array.length a - 1) k =
    if l > u then -1. else
    let m = l + (u - l) / 2 in
    match compare k a.(m) with
    | -1 -> search a ~l ~u:(m - 1) k
    | 1 -> search a ~l:(m + 1) ~u k
    | _ -> 1. ;;

let extend_one_one gen_classifier =
  printf "Training 1-1..\n%!";
  let t0 = Sys.time() in
  let grouped_data_by_category = Enum.foldi (fun i l acc -> Map.modify_def [] l (List.cons i) acc) (Map.create Int.compare) (Array1.enum train_labels) in
  let gdc = Map.map (Array.of_list |- tap (Array.sort Int.compare)) grouped_data_by_category in
  let category_pairs = 
    Map.keys gdc |> map (fun i -> Map.keys gdc // (fun x -> x > i) /@ (fun j -> i,j)) 
    |> Enum.flatten |> List.of_enum 
  in
  let datapoints i j =
    let get c = Map.find c gdc |> Array.enum in
    Enum.append (get i) (get j) |> Random.shuffle
  in
  let classes i data = let pos = Map.find i gdc in Array.map (search pos) data in
  let make_classifier (i,j) = let data = datapoints i j in gen_classifier data (classes i data) in
  let ps = List.map make_classifier category_pairs |> Array.of_list in
  printf "Done training(%.2fs)\n%!" (Sys.time () -. t0);
  One_one (category_pairs, ps)

(***************************************)
(** Prediction functions ***************)
(***************************************)

let predict_b = 
  function
  | Dot warr -> let w = Array1.of_array Datafile.kind Datafile.layout warr in (fun x -> dot w x)
  | Kern_3 (offs, ais) ->
    let len = Array.length offs in
(*    printf "Kern3 decoder (off %d; len %d, %d ais)\n" off len (Array.length ais); *)
    let ais = Array1.of_array Datafile.kind Datafile.layout ais in
    (fun x -> 
      let t = ref 0. in 
      for i = 1 to len do 
	t := !t +. ais.{i} *. k_3 x (Array2.slice_right train_data offs.(i)); 
      done; 
      !t +. 0.)
  | Kern_rbf (sigma, offs, ais) -> 
    let len = Array.length offs in
    let ais = Array1.of_array Datafile.kind Datafile.layout ais in
    (fun x -> 
      let t = ref 0. in 
      for i = 1 to len do 
	t := !t +. ais.{i} *. k_rbf sigma x (Array2.slice_right train_data offs.(i)); 
      done; 
      !t +. 0.)
  | Svm _fn -> assert false

let get_heads enum_l = 
  Array.fold_left (fun acc x -> match Enum.get x with 
      None -> raise Enum.No_more_elements 
    | Some v -> v :: acc) 
    [] enum_l |> List.rev

let rec popcount c x = if x = 0 then c else popcount (c+1) (x land (x-1))
let ham_dist x y = popcount 0 (x lxor y) (*|> tap (fun d -> printf "HD(%x,%x) = %d " x (Bit_cat.to_int y) d) *)

let nearest_in label_map pred_lab = 
    (1--category_count) |> Enum.arg_min (fun i -> ham_dist pred_lab label_map.(i))

let rec to_bits n x = if n = 0 then [] else (if x land 1 = 1 then 1. else -1.) :: (to_bits (n-1) (x asr 1))
let merge_qb (q,p) pred = (q *. abs_float pred, if pred >= 0. then p lsl 1 + 1 else p lsl 1)

let predict_cat = function
  | Hamm (bpreds, cat_map) ->
    let decode_arr = Array.init (1 lsl !cat_bits) (fun i -> nearest_in cat_map i) in
    let classifiers = Array.map (fun bp -> predict_b bp) bpreds in
    let decode_bits l = 
      let (str, bits) = List.fold_left merge_qb (1.,0) l in
      str, decode_arr.(bits)
    in
    (fun xs -> 
      let decisions = Array.map (fun cl -> predict_many xs cl) classifiers in
      Enum.from (fun () -> get_heads decisions) |> Enum.map decode_bits
    )
  | One_one (pairs,bpreds) ->
    let classifiers = Array.map (fun bp -> predict_b bp) bpreds in
    let decode decisions = 
      let votes = Array.create (category_count+1) 0 in
      List.iter2 (fun (i,j) d -> if d > 0. then votes.(i) <- votes.(i) + 1 else votes.(j) <- votes.(j) + 1) pairs decisions;
      let winner = Enum.arg_max (fun i -> votes.(i)) (1--category_count) in
      float votes.(winner) /. float category_count, winner
    in
    (fun xs -> 
      let decisions = Array.map (fun cl -> predict_many xs cl) classifiers in
      Enum.from (fun () -> get_heads decisions) |> Enum.map decode
    )

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

let push ((k:float),v) (map,min_k,count as acc) = 
  if count < top_n then 
    (Map.add k v map, min k min_k, count+1) 
  else if k < min_k then acc 
  else 
    let m = Map.remove min_k map |> Map.add k v in 
    (m, Map.min_binding m |> fst, count)

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
  printf "Scores: %a\n" (Map.print ~first:"" ~last:"" ~sep:", " Int.print print_float) scores;
  let overall = Map.enum scores |> map snd |> Enum.reduce (+.) in
  overall /. float category_count

let print_pred oc (i,ps) = 
  List.iter (fun p -> fprintf oc "%d\t%d\n" i p) ps

let output_preds oc ps =
  Enum.print ~first:"" ~last:"" ~sep:"" print_pred oc (Map.enum ps)

let run_test ?(n = 10_000) name (cpred: cpredictor) = 
  let rand_offset = Random.int (Array2.dim2 test_data - n) in
  Array2.sub_right test_data rand_offset n |> predict_cat cpred 
  |> best_predictions |> evaluate |> printf "%s Overall: %.2f\n%!" name
 
let marshal_file fn_base x =
  let tm = Unix.time() |> Unix.localtime in
  let fn = sprintf "preds/%s.%d.%d.%d" fn_base tm.Unix.tm_mday tm.Unix.tm_hour tm.Unix.tm_min in
  let ctr = ref 0 in
  let fn = 
    if Sys.file_exists (fn ^ ".pred") then 
      ( while Sys.file_exists (fn ^ "." ^ string_of_int !ctr ^ ".pred") do incr ctr; done; 
	(fn ^ "." ^ string_of_int !ctr ^ ".pred")) else (fn ^ ".pred")
  in
  let oc = Pervasives.open_out fn in
  Legacy.Marshal.to_channel oc x [];
  Pervasives.close_out oc

let crossval partial_pred param_values =
    List.iter (fun i -> partial_pred i |> run_test ("crossvalidate:" ^ string_of_float i)) param_values

let crossval_int partial_pred param_values =
    List.iter (fun i -> partial_pred i |> run_test ("crossvalidate:" ^ string_of_int i)) param_values

(***************************************)
(**        MAIN        *****************)
(***************************************)

let train () = 
(*  kperceptron3_slice 1 2000 |> extend_hamm |> marshal_file "hamm_kp3_0_2k"; *)
(*  train_slices 1000 kperceptron3_slice |> List.iter (extend_hamm |- marshal_file "hamm_kp3_slc_1k"); *)
  kperceptron_rbf_slice 50. 1 2000 |> extend_hamm |> marshal_file "hamm_kpr_1_2k";
  train_slices 1000 (kperceptron_rbf_slice 50.) |> List.iter (extend_hamm |- marshal_file "hamm_kpr_slc_1k")

(*
  extend_hamm_incr perceptron |> train_full |>
      run_test ~print:real_run "Hamm";
  extend_one_one_incr perceptron |> train_full |>
      run_test ~print:real_run "1-1 ";
 *)

let print_pred_conf oc i (conf, pred) =
  fprintf oc "%d\t%f\t%d\n" pred conf i
let predict p oc =
  predict_cat p test_data |> Enum.iteri (print_pred_conf oc)

let pred_read fn =
  let ic = Pervasives.open_in_bin fn in
  let (p: cpredictor) = Legacy.Marshal.from_channel ic in
  Pervasives.close_in ic;
  p

let predict_scan () = ()
  (* scan for *.pred files, move them to *.pred_act, create a *.output file and write predictions there *)
  (* rename to used.* when done *)

let test () = 
(* best: 2K?
  crossval_int (fun i ->   
    let rand_offset = Random.int (Array2.dim2 train_data - i) in
    kperceptron3_slice rand_offset i |> extend_hamm) [500; 1000; 2000; 4000; 8000; 10000; 15000; 20000]
*)
  crossval_int (fun i ->   
    let rand_offset = Random.int (Array2.dim2 train_data - 2000) in
    kperceptron_pow_slice i rand_offset 2000 |> extend_hamm) [1; 2; 3; 4; 5; 6; 7; 8; 9; 10]



(*  pred_read "kp3_0_2k.pred" |> run_test "kp3_0_2k" *)
    
