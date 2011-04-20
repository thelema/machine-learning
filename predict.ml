open Batteries_uni
open Printf

open Libosvm 
open Lacaml.Impl.S (* Single-precision reals *)
open Bigarray

type matrix = (float, float32_elt, fortran_layout) Array2.t
type mat64 = (float, float64_elt, fortran_layout) Array2.t
type cats = (int, int8_unsigned_elt, fortran_layout) Array1.t
type 'a pred_f_t = matrix -> (float * 'a) Enum.t

let print_float oc x = fprintf oc "%.2f" x
let print_float5 oc x = fprintf oc "%.5f" x

let debug = false
let category_count = 164
let top_n = 100

let set_mhs n = Gc.set { (Gc.get()) with Gc.minor_heap_size = n; } 
let () = set_mhs 1_000_000; Random.self_init ();;

let train_data = Datafile.get_matrix "training.ba"
let train_labels = Datafile.read_label_file "training_label.txt"
let train_rows = Array1.dim train_labels

let train_data64 = lazy(Datafile.get_matrix64 "training.txt.ba64")
let d64_slice = Datafile.get_matrix64 "training.txt.ba64.slice"
let d64_labels = Datafile.read_label_file "training_label.txt.slice"

let pred_split p = abs_float p, (if p=0. then Random.bool () else p>0.)

let norm ais = Vec.sqr_nrm2 ais |> sqrt
let zeros ais = (1--Array1.dim ais) |> Enum.fold (fun acc i -> if ais.{i} = 0. then acc+1 else acc) 0
let get_i (m:matrix) i = Array2.slice_right m i
let vec_of_arr a = Array1.of_array Datafile.kind Datafile.layout a
let vec64_of_arr a = Array1.of_array Datafile.kind64 Datafile.layout a
let mat_of_arr a = Array2.of_array Datafile.kind64 Datafile.layout a

let matrix_print oc (m:mat64) = 
  for i = 1 to Array2.dim1 m do
    for j = 1 to Array2.dim2 m do
      fprintf oc "%.2f " m.{i,j};
    done;
    fprintf oc "\n";
  done

type bpredictor = 
  | Dot of float array (* weights of different features *)
  | Kern_3 of int array * float array (* offsets of data used, ais*)
  | Kern_rbf of float * int array * float array (* sigma, offsets of data used, ais*)
  | Kern_pow of float * int array * float array (* exponent, offsets of data used, ais*)
  | Dot_plus of float array * float
  | Notest
  | Svm_b of string (* filename of saved model *)

type cpredictor = 
  | Hamm of bpredictor array * int array (* n predictors, n-bit codewords *)
  | One_one of (int * int) array * bpredictor array (* i vs j pairs, with a predictor for each *)
  | Svm of string (* filename of saved model *)
  | Hedge of cpredictor array * float array
  | Pre_predicted of string (* filename with strength\tprediction lines *)

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
    let xi = get_i data i in
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

let perceptron_b offsets labels =
  let w = Vec.make0 Datafile.cols in
  let b = ref 0. in
  let acc = Vec.make0 Datafile.cols in
  let bacc = ref 0. in
  let n = Array.length offsets in
  printf "pb(%d)%!" n; 
  for ti = 1 to n do
    if ti land 0xfff = 0 then printf ".%!";
    let xi = get_i train_data offsets.(ti-1) in
    let yi = labels.(ti-1) in
    if !b +. dot xi w *. yi > 0. then incr right
    else (axpy ~alpha:yi ~x:xi w; b := !b +. 0.2 *. yi;);
    let n2 = Vec.sqr_nrm2 w +. !b *. !b in
    if n2 > 1000. then (scal (800. /. n2) w; b := (800. /. n2) *. !b);
    ignore(Vec.add ~z:acc acc w); bacc := !bacc +. !b;
  done;
  count := !count + n;
  scal (1. /. float n) acc;
  bacc := !bacc /. float n;
  Dot_plus (Array1.to_array acc, !bacc)

let batch_perb os ls = 
  Array.map (perceptron_b os) ls

let rec pred (kij: matrix) (ais: vec) i acc j = 
  if j < 1 then acc else 
    pred kij ais i (ais.{j} *. kij.{j,i} +. acc) (j-1) 

let kp_core kij (ais: vec) (labels: vec) =
  let n = Array1.dim ais in
  let errs = ref 0 in
  for i = 1 to n do
    let yi = labels.{i} in
    if yi *. (pred kij ais i 0. n) <= 0. then 
      (incr errs; ais.{i} <- ais.{i} +. yi)
  done;
  !errs

let cap_psi sr phi mu = (sr *. phi) *. (sr *. phi) +. 2. *. sr *. phi *. (1. -. phi *. mu)

(*
let rec sparse_pred_aux (kij: matrix) (ais: vec) i acc j = 
  if ais.(j) = 0. then sparse_pred_aux kij ais i acc (j+1) 
  else sparse_pred_aux kij ais i (ais.{j} *. kij.{j,i} +. acc) t
 *)

(* magical formula for phi *)
let solve_phi sr mu q m = 
  let a = sr *. sr -. 2. *. sr *. mu in
  let b = 2. *. sr in
  let c = q -. (15. /. 32.) *. float m in
  if a = 0. then -.c /. b else 
    (-. b +. sqrt (b *. b -. 4. *. a *. c) ) /. (2. *. a) 

let rec nonzeros_at_least a lim i = 
  if i >= Array.length a then false else
  nonzeros_at_least a (if a.(i) = 0. then lim else lim-1) (i+1)

let shuffle a = 
  for n = Array.length a - 1 downto 1 do
    let k    = Random.int ( n + 1 ) in
    if k <> n then
      let buf  = Array.get a n in
      Array.set a n (Array.get a k);
      Array.set a k buf;
  done


  (* Implementation of forgetron, a bounded memory perceptron *)
let rec ft_core ~b ?(b0=10) ?(bincr=10) (kij: matrix) (ais: vec) =
  let n = Array1.dim ais in
  let forget_queue = ref Deque.empty in
  let q = ref 0. in (* sum of all cap_psi so far *)
  let m = ref 0 in (* total mistakes *)
  let bnow = ref (b0 - bincr) in
  printf "ft(%d)%!" n; 
  let process_order = Array.init n (fun i -> i+1) in
  fun (labels: float array) ->
    let m0 = !m in
    shuffle process_order;
    bnow := min b (!bnow + bincr);
    for i = 0 to n-1 do
      let i = process_order.(i) in
      (*      if i land 0xfff = 0 then printf ".%!"; *)
      let yi = labels.(i-1) in
      let mu_i = dot (get_i kij i) ais in
      if yi *. mu_i <= 0. then ( (* guessed wrong or no guess *)
	if ais.{i} = 0. then ( (* i is not in queue *)
	  (* put i in the queue *)
	  forget_queue := Deque.cons i !forget_queue;
	  if Deque.size !forget_queue > !bnow then (* overflow *)
	    (* r is oldest in queue, will be removed *)
	    let fq, r = Option.get (Deque.rear !forget_queue) in
	    forget_queue := fq;
	    (* the weight for the oldest *)
	    let sr = abs_float ais.{r} in 
	    ais.{i} <- yi; (* add i to ais *)
	    (* the current prediction for x_r *)
	    let mu = labels.(r-1) *. dot (get_i kij r) ais in 
	    ais.{r} <- 0.; (* remove r from ais *)
	    let phi = solve_phi sr mu !q !m in  (* optimal phi *)
	    if phi <= 0. then failwith "negative phi";
	    if phi < 1. then (
	      scal phi ais; (* scale ais down by phi *)
	      q := !q +. cap_psi sr phi mu; (* accumulate q *)
	    ) else (* don't scale, cap phi at 1.0 *)
	      q := !q +. cap_psi sr 1.0 mu; (* accumulate q *)
	);
	incr m;
	ais.{i} <- yi (* no matter what, include i *)
      )
    done;
    (*    printf "(q:%.2f, m:%d)" !q !m; *)
    !m - m0 (* the number of mistakes this loop *)


(* Implementation of forgetron, a bounded memory perceptron *)
(* this version doesn't use a gram matrix, but computes the kernel on the fly *)
let rec ft_core_onfly ~b ?(b0=10) (k: vec -> vec -> float) (offsets: int array) (ais: vec) =
  let n = Array1.dim ais in
  let forget_queue = ref Deque.empty in
  let q = ref 0. in (* sum of all cap_psi so far *)
  let m = ref 0 in (* total mistakes *)
  let bincr = 10 in
  let bnow = ref (b0 - bincr) in
  printf "ft(%d)%!" n; 
  let f x = 
    let acc = ref 0. in 
    for i = 1 to n do 
      if ais.{i} <> 0. then 
	acc := !acc +. ais.{i} *. k (get_i train_data offsets.(i)) x
    done;
    !acc +. 0.
  in
  let process_order = Array.init n (fun i -> i) in
  fun (labels: float array) ->
    let m0 = !m in
    shuffle process_order;
    bnow := min b (!bnow + bincr);
    for i = 0 to n-1 do
      let i = process_order.(i) in (* 0..n-1 -> 0..n-1 *)
      (*      if i land 0xfff = 0 then printf ".%!"; *)
      let yi = labels.(i) in
      let mu_i = f (get_i train_data offsets.(i)) in
      if yi *. mu_i <= 0. then ( (* guessed wrong or no guess *)
	if ais.{i+1} = 0. then ( (* i is not in queue *)
	  (* put i in the queue *)
	  forget_queue := Deque.cons i !forget_queue;
	  if Deque.size !forget_queue > !bnow then (* overflow *)
	    (* r is oldest in queue, will be removed *)
	    let fq, r = Option.get (Deque.rear !forget_queue) in
	    forget_queue := fq;
	    (* the weight for the oldest *)
	    let sr = abs_float ais.{r} in 
	    ais.{i+1} <- yi; (* add i to ais *)
	    (* the current prediction for x_r *)
	    let mu = labels.(r) *. f (get_i train_data r) in 
	    ais.{r+1} <- 0.; (* remove r from ais *)
	    let phi = solve_phi sr mu !q !m in  (* optimal phi *)
	    if phi <= 0. then failwith "negative phi";
	    if phi < 1. then (
	      scal phi ais; (* scale ais down by phi *)
	      q := !q +. cap_psi sr phi mu; (* accumulate q *)
	    ) else (* don't scale, cap phi at 1.0 *)
	      q := !q +. cap_psi sr 1.0 mu; (* accumulate q *)
	);
	incr m;
	ais.{i+1} <- yi (* no matter what, include i *)
      )
    done;
    (*    printf "(q:%.2f, m:%d)" !q !m; *)
    !m - m0 (* the number of mistakes this loop *)

(* RBF kernel function and Gram Matrix generator *)	
let k_rbf two_s_sqr x1 x2 = exp (~-. (Vec.ssqr_diff x1 x2) /. two_s_sqr)
let gen_kij_rbf ?(debug=false) two_sig_sq offsets =
  let n = Array.length offsets in
  let kij = Array2.create float32 fortran_layout n n in
  if debug then printf "rbf(%d^2)%!" n;
  for i = 1 to n do
    if debug && i land 0xff = 0 then printf ".%!";
    let xi = (get_i train_data offsets.(i-1)) in
    for j = i to n do
      let k = k_rbf two_sig_sq xi (get_i train_data offsets.(j-1)) in
      kij.{i,j} <- k; kij.{j,i} <- k
    done;
  done;
  if debug then printf "Done\n%!";
  kij

(* x^3 kernel function and Gram Matrix generator *)
let k_3 x1 x2 = let a = (1. +. dot x1 x2) in a *. a *. a
let gen_kij_3 offsets =
  let n = Array.length offsets in
  let kij = Array2.create float32 fortran_layout n n in
  printf "k3(%d^2)%!" n;
  for i = 1 to n do
    if i land 0xff = 0 then printf ".%!";
    let xi = (get_i train_data offsets.(i-1)) in
    for j = i to n do
      let k = k_3 xi (get_i train_data offsets.(j-1)) in
      kij.{i,j} <- k; kij.{j,i} <- k
    done;
  done;
  kij

let k_pow pow x1 x2 = let a = (1. +. dot x1 x2) in a ** pow
let gen_kij_pow ?(debug=false) pow offsets =
  let n = Array.length offsets in
  let kij = Array2.create float32 fortran_layout n n in
  if debug then printf "x^%f(%d)%!" pow n;
  for i = 1 to n do
    if i land 0xff = 0 then printf ".%!";
    let xi = (get_i train_data offsets.(i-1)) in
    for j = 1 to n do
      let k = k_3 xi (get_i train_data offsets.(j-1)) in
      kij.{i,j} <- k; kij.{j,i} <- k
    done;
  done;
(*  printf "Done\n%!"; *)
  kij


let clean2 pred_a a b =
  let n = Array.length a in
  (* Use a bitset to store which elements will be in the final array. *)
  let bs = BatBitSet.create n in
  for i = 0 to n-1 do
    if pred_a a.(i) then BatBitSet.set bs i
  done;
  (* Allocate the final array and copy elements into it. *)
  let new_len = BatBitSet.count bs in
  let j = ref 0 in
  let filtered_a = Array.init new_len
    (fun _ ->
       (* Find the next set bit in the BitSet. *)
       while not (BatBitSet.is_set bs !j) do incr j done;
       let r = a.(!j) in
       incr j;
       r) in
  j := 0;
  let filtered_b = Array.init new_len
    (fun _ ->
       (* Find the next set bit in the BitSet. *)
       while not (BatBitSet.is_set bs !j) do incr j done;
       let r = b.(!j) in
       incr j;
       r) in
  printf "(%d) %!" new_len;
  (filtered_a, filtered_b)

let clean (offs: int array) a =
  let a = Array1.to_array (a: vec) in
  clean2 (fun a -> a <> 0.) a offs

let kperceptron_offs (genk,tag_v) (core: matrix -> vec -> float array -> int) ~loops offs labels_arr =
  let n = Array.length offs in
  let kij = genk offs in
  let run_perc labels =
    let ais = Vec.make0 (Array.length offs) in
    let core = core kij ais in
    let l = ref 0 in
    let cutoff = n / 10 in
    let wrongs = ref n in
    while !l < loops && !wrongs > cutoff do 
      incr l; wrongs := core labels 
    done;
    if !l >= loops then printf "F%d" !wrongs else printf "%d" !l; 
    clean offs ais |> tag_v
  in  
  Array.map run_perc labels_arr

let kperceptron_slice gk core ~loops (off, len) labels_arr =
  let offs = (off --^ (off + len) |> Array.of_enum) in 
  kperceptron_offs gk core ~loops offs labels_arr

let kperceptron_elems (genk, tag_v) core ~loops offs labels =
  let n = Array.length offs in
  if Array.length labels <> n then invalid_arg "Labels must have the same length as offs";
  let (kij:matrix) = genk offs in
  let ais = Vec.make0 n in
  let core = core kij ais in
  let l = ref 0 in
  let cutoff = n / 10 in
  let wrongs = ref n in
  while !l < loops && !wrongs > cutoff do 
    incr l; wrongs := core labels 
  done;
  if !l >= loops then printf "F%d" !wrongs else printf "%d" !l; 
  clean offs ais |> tag_v

let gt_k3 = (gen_kij_3, fun (a,o) -> Kern_3 (o,a))
let gt_rbf sigma = (gen_kij_rbf sigma, fun (a,o) -> Kern_rbf(sigma,o,a))
let gt_pow p = (gen_kij_pow p, fun (a,o) -> Kern_pow (p, o, a))

let kt_rbf sigma = (k_rbf sigma, fun (a,o) -> Kern_rbf(sigma, o, a))

let kperceptron_onfly (k, tag_v) core_fly ~loops offs labels =
  let n = Array.length offs in
  if Array.length labels <> n then invalid_arg "Labels must have the same length as offs";
  let ais = Vec.make0 n in
  let core = core_fly k ais in
  let l = ref 0 in
  let cutoff = n / 40 in
  let wrongs = ref n in
  while !l < loops && !wrongs > cutoff do 
    incr l; wrongs := core labels 
  done;
  if !l >= loops then printf "F%d" !wrongs else printf "%d" !l; 
  clean offs ais |> tag_v

  
let print_statistics ?(oc=stdout) e =
  match Enum.get e with
      None -> fprintf oc "Empty enumeration - no statistics\n"
    | Some x0 ->
	let m = ref x0
	and k = ref 1 
	and s = ref 0. in
	let t = ref 0. in
	Enum.iter (fun x -> 
		     t := !t +. x;
		     incr k; 
		     let mk = !m +. (x -. !m)/.(float !k) in 
		     s := !s +. (x -. !m) *. (x -. mk);
		     m := mk) e;
	let stdev = sqrt(!s /. float (!k-1)) in
	fprintf oc "N: %d Sum: %.2f Mean: %.2f Stdev: %.2f\n" !k !t !m stdev


(* Kernel least squares  ~sigma is 2*sigma^2 *)
let klsq ~lambda ~sigma ~cutoff d l = 
  let k = gen_kij_rbf sigma d in
  for i = 1 to Array.length d do
    k.{i,i} <- k.{i,i} +. lambda
  done;
  let b = Mat.of_col_vecs [|vec_of_arr l|] in
  getrs k b;
  let a = Array2.slice_right b 1 |> Array1.to_array in
  print_statistics (Array.enum a);
  let a,o = clean2 (fun ai -> abs_float ai > cutoff) a d in
  Kern_rbf (sigma, o, a)


let svm (l: vec) = 
  let fn = Filename.temp_file ~temp_dir:"svm/" "svmc" "" in 
  let m = svm_train ~kernel_type:RBF d64_slice (Obj.magic l) in
  svm_save_model ~file:fn ~m;
  Svm fn

let genmat n = 
  Array2.create Datafile.kind64 Datafile.layout Datafile.cols n
let dbuf = ref (genmat 1)

let project_data_rows64 offs =
  let n = Array.length offs in
  if Array2.dim2 !dbuf <> n then dbuf := genmat n;
  let dout = !dbuf in
  let d = Lazy.force train_data64 in
  let out_pos = ref 1 in
  Array.iter (fun i -> Array1.blit (Array2.slice_right d i) (Array2.slice_right dout !out_pos); incr out_pos) offs;
  dout
  

let svm offs labels = 
  let d = project_data_rows64 offs in
  let fn = Filename.temp_file ~temp_dir:"svm/" "svmb" "" in 
  let m = svm_train ~kernel_type:RBF d (Obj.magic (vec64_of_arr labels)) in
  svm_save_model ~file:fn ~m;
  Svm_b fn

    

(*
let () = printf "Reading..%!"
let t0 = Sys.time()
let v1 = svm_scale ~lower:0. alldata
let m = svm_train ~kernel_type:LINEAR v1
let () = printf "Done(%.2f)\n%!" (Sys.time () -. t0)
 *)

(*

let svm (data: matrix64) (labels: vec) = 
  let _model = svm_train ~svm_type:C_SVC ~kernel_type:RBF data (* labels *) in
  assert false
*)
(********************************************)
(** EXTENDING 2-category to multi-category **)
(********************************************)

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

let extend_hamm cat_bits gen_classifier =
  printf "Training Hamming Classifiers%!";
  let n = train_rows in
  let mask = (1 lsl cat_bits) - 1 in
  let {get = map_cat; enum=cat_mapping} = 
    make_map (fun _ -> Random.full_range () land mask) in
  let bit_labels i x = if ((map_cat x) asr i) land 1 = 1 then 1. else -1. in
  let gen_labels i = Array.init n (fun j -> bit_labels i train_labels.{j+1}) in
  printf ".%!";
  let all_labels = Array.init cat_bits gen_labels in
  printf ".%!";
  let cat_map = Array.create (category_count+1) 0 in
  printf ".%!";
  Enum.iter (fun (i,b) -> cat_map.(i) <- b) (cat_mapping ());
  printf ".\n%!";
  let t0 = Sys.time () in
  let classifiers = gen_classifier all_labels in
  printf "Done training (%.2fs)\n%!" (Sys.time () -. t0);
  Hamm (classifiers, cat_map)

(* group all datapoints in training set into categories, with max [cap] values in each group *)
let group_by_cat ?(category_count=category_count) cap =
  let cats = Array.create (category_count+1) Vect.empty in
  let insert_random l i = 
    let len = Vect.length cats.(l) in
    let pos = Random.int (1+len) - 1 in
    if pos = -1 then 
      cats.(l) <- Vect.prepend i cats.(l) 
    else if pos = len-1 then 
      cats.(l) <- Vect.append i cats.(l)
    else
      cats.(l) <- Vect.insert pos (Vect.singleton i) cats.(l)
  in
  for i = 1 to train_rows do
    let l = train_labels.{i} in 
    if l < category_count then insert_random l i;
  done;
  let cap_slice v = 
    Vect.to_array (if Vect.length v > cap then Vect.sub 0 cap v else v)
  in
  Array.map cap_slice cats

let extend_one_one cats_a rejector gen_classifier =
  printf "Training 1-1..\n%!";
  let cap = Array.fold_left (fun acc a -> let l = Array.length a in if l > acc then l else acc) 0 cats_a in
  let t0 = Sys.time() in
  let cat_pairs = ref Vect.empty in
  let ps = ref Vect.empty in
  let common_buf = Array.create (cap * 2) 0 in
  let common_label = Array.create (cap * 2) 0. in
  for i = 1 to category_count do
    for j = i+1 to category_count do
      let ilen = Array.length cats_a.(i) in
      let jlen = Array.length cats_a.(j) in
      let data, labels = 
	if ilen + jlen = cap * 2 then ( (* use common buf *)
	  Array.blit cats_a.(i) 0 common_buf 0 ilen;
	  Array.blit cats_a.(j) 0 common_buf ilen jlen;
	  Array.fill common_label 0 ilen (1.);
	  Array.fill common_label ilen jlen (-1.);
	  common_buf, common_label
	) else (
	  Array.append cats_a.(i) cats_a.(j), 
	  Array.init (ilen + jlen) (fun i -> if i < ilen then 1. else -1.)
	) in
      let t1 = Unix.gettimeofday () in
      let cij = gen_classifier data labels in
      let t2 = Unix.gettimeofday () in
      printf "gen:%.2f\t" (t2-.t1);
      if not (rejector cij) then (
	ps := Vect.append cij !ps;
	cat_pairs := Vect.append (i,j) !cat_pairs;
      )
    done
  done;
  let category_pairs = Vect.to_array !cat_pairs in
  let ps = Vect.to_array !ps in
  printf "Done training(%.2fs)\n%!" (Sys.time () -. t0);
  One_one (category_pairs, ps)

(***************************************)
(***********     SVM      **************)
(***************************************)


(** TODO **)



(***************************************)
(** Prediction functions ***************)
(***************************************)

let predict_b = 
  function
  | Dot warr -> let w = vec_of_arr warr in (fun x -> dot w x)
  | Dot_plus (warr, b) -> let w = vec_of_arr warr in (fun x -> dot w x +. b)
  | Kern_3 (offs, ais) ->
(*    printf "Kern3 decoder (off %d; len %d, %d ais)\n" off len (Array.length ais); *)
    (fun x -> 
      let len = Array.length offs in
      let t = ref 0. in 
      for i = 0 to len-1 do 
	t := !t +. ais.(i) *. k_3 x (get_i train_data offs.(i)); 
      done; 
      !t +. 0.)
  | Kern_rbf (sigma, offs, ais) -> 
    (fun x -> 
      let len = Array.length offs in
      let t = ref 0. in 
      for i = 0 to len-1 do 
	t := !t +. ais.(i) *. k_rbf sigma x (get_i train_data offs.(i)); 
      done; 
      !t +. 0.)
  | Kern_pow (p, offs, ais) -> 
    (fun x -> 
      let len = Array.length offs in
      let t = ref 0. in 
      for i = 1 to len do 
	t := !t +. ais.(i-1) *. k_pow p x (get_i train_data offs.(i-1)); 
      done; 
      !t +. 0.)
  | Notest -> (fun x -> 0.)
  | Svm_b fn -> 
    let m = svm_load_model fn in
    (fun x -> svm_predict_probability m (mat_of_arr [|Array1.to_array x|]) |> Pair.print2 matrix_print stdout; 0.)

let print_bpred oc = function
  | Dot w -> Array.length w |> fprintf oc "Dot(%d)"
  | Dot_plus (w,b) -> fprintf oc "Dot(%d+%.2fb)" (Array.length w) b
  | Kern_3 (o,a) -> fprintf oc "K3(%d)" (Array.length o)
  | Kern_rbf (s,o,a) -> fprintf oc "RBF%.2f(%d)" s (Array.length o)
  | Kern_pow (p,o,a) -> fprintf oc "K%.1f(%d)" p (Array.length o)
  | Notest -> fprintf oc "NT"
  | Svm_b fn -> fprintf oc "Svmb(%s)" fn

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
let predict_many f d = 
  let t0 = Unix.gettimeofday () in 
  let r = Array.init (Array2.dim2 d) (fun i -> if i land 0xfff = 0 then printf ".%!"; f (get_i d (i+1))) in
  printf "%.0f/s" (float (Array2.dim2 d) /. (Unix.gettimeofday () -. t0));
  r

let rec predict_cat = function
  | Hamm (bpreds, cat_map) ->
    let classifiers = Array.map (fun bp -> predict_b bp) bpreds in
    let bits = Array.length bpreds in
    let decode_bits = 
      if bits < 20 then (* precompute decoding map when it's not too large *)
	let decode_arr = Array.init (1 lsl bits) (fun i -> nearest_in cat_map i) in
	(fun l -> let (str, bits) = Array.fold_left merge_qb (1.,0) l in str, decode_arr.(bits))
      else
	(fun l -> let (str, bits) = Array.fold_left merge_qb (1.,0) l in str, nearest_in cat_map bits)
    in
    (fun (d:vec) -> Array.map (fun cl -> cl d) classifiers |> decode_bits)
  | One_one (pairs,bpreds) ->
    let classifiers = Array.map (fun bp -> predict_b bp) bpreds in
    let decode decisions = 
      let votes = Array.create (category_count+1) 0 in
      Array.iter2 (fun (i,j) d -> if d > 0. then votes.(i) <- votes.(i) + 1 else votes.(j) <- votes.(j) + 1) pairs decisions;
      let winner = Enum.arg_max (fun i -> votes.(i)) (1--category_count) in
      float votes.(winner) /. float category_count, winner
    in
    (fun d -> Array.map (fun cl -> cl d) classifiers |> decode)
  | Svm fn -> 
    let m = svm_load_model fn in
    (fun x -> svm_predict_probability m (mat_of_arr [|Array1.to_array x|]) |> Pair.print2 matrix_print stdout; (0., 1))
  | Hedge (es, ws) ->
    let classifiers = Array.map (fun p -> predict_cat p) es in
    let weight_sum = Array.reduce (+.) ws in
    let decode ds =
      let votes = Array.create (category_count + 1) 0. in
      Array.iteri (fun i (str,l) -> votes.(l) <- votes.(l) +. ws.(i) *. str) ds;
      let winner = Enum.arg_max (fun i -> votes.(i)) (1--category_count) in
      votes.(winner) /. weight_sum, winner
    in
    (fun d -> Array.map (fun cl -> cl d) classifiers |> decode)
  | Pre_predicted fn ->
    let parse _ = assert false in
    let ps = File.lines_of fn |> map parse |> Array.of_enum in
    let i = ref 0 in
    (fun _ -> let r = ps.(!i) in incr i; r)

let rec batch_predict_cat = function
  | Hamm (bpreds, cat_map) ->
    let classifiers = Array.map (fun bp -> predict_b bp) bpreds in
    let bits = Array.length bpreds in
    let decode_bits = 
      if bits < 20 then (* precompute decoding map when it's not too large *)
	let decode_arr = Array.init (1 lsl bits) (fun i -> nearest_in cat_map i) in
	(fun bits -> decode_arr.(bits))
      else
	(fun bits -> nearest_in cat_map bits)
    in
    (fun (m: matrix) -> 
      let n = Array2.dim2 m in
      let strs = Array.create n 1. in
      let results = Array.create n 0 in 
      for i = 0 to bits-1 do 
	for j = 0 to n-1 do
	  let pred = classifiers.(i) (get_i m (j+1)) in
	  strs.(j) <- strs.(j) *. abs_float pred;
	  results.(j) <- if pred >= 0. then results.(j) lsl 1 + 1 else results.(j) lsl 1;
	done;
      done;
      Array.map2 (fun s r -> s, decode_bits r) strs results)
  | One_one (pairs,bpreds) ->
    let classifiers = Array.map (fun bp -> predict_b bp) bpreds in
    let decode votes = 
      let winner = Enum.arg_max (fun i -> votes.(i)) (1--category_count) in
      float votes.(winner) /. float category_count, winner
    in
    (fun (m: matrix) -> 
      let n = Array2.dim2 m in
      let results = Array.create_matrix n (category_count+1) 0 in 
      for i = 0 to Array.length pairs - 1 do 
	for j = 1 to n do
	  let pred = classifiers.(i) (get_i m j) in
	  let (a,b) = pairs.(i) in
	  if pred >= 0. 
	  then results.(j).(a) <- results.(j).(a) + 1
	  else results.(j).(b) <- results.(j).(b) + 1
	done;
      done;
      Array.map decode results)
  | Svm _file -> assert false
(*    let m = svm_load_model ~file in
    (fun xs -> svm_predict_32 ~m ~x:xs) *)
  | Hedge (es, ws) -> assert false
  | Pre_predicted fn -> assert false


let rec print_cpred oc = function
  | Hamm (bps,cmap) -> fprintf oc "Hamm(%ax%d)" print_bpred bps.(0) (Array.length bps)
  | One_one (pairs, bpreds) -> fprintf oc "OVO(%ax%d)" print_bpred bpreds.(0) (Array.length bpreds)
  | Svm _ -> fprintf oc "Svm()"
  | Hedge (es, ws) -> 
    fprintf oc "Hedge("; 
    for i = 0 to Array.length es - 1 do 
      fprintf oc "%a,%.3f" print_cpred es.(i) ws.(i) 
    done;
    fprintf oc ")"
  | Pre_predicted fn -> fprintf oc "PP(%s)" fn

(***************************************)
(** COMBINING MULTIPLE CPREDICTORS *****)
(***************************************)

let hedge experts (off, len) =
  let n = Array.length experts in
  let w = Array.create n 0. in
  let preds = Array.map predict_cat experts in
  let votes = Array.create (category_count + 1) 0. in
  let strs = Array.create n 0. in
  let guesses = Array.create n 0 in
  for i = off to off+len-1 do
    let xi = get_i train_data i in
    Array.fill votes 0 n 0.;
    for i = 0 to n-1 do
      let str,l = preds.(i) xi in
      strs.(i) <- str; guesses.(i) <- l;
      votes.(l) <- votes.(l) +. exp w.(i) *. str;
    done;
    let joint_prediction = Enum.arg_max (fun i -> votes.(i)) (1--category_count) in
    let yi = train_labels.{i} in
    if joint_prediction <> yi then
      let test_expert j = if strs.(j) > 0.01 && guesses.(j) <> yi then w.(j) <- w.(j) -. strs.(j) in
      for i = 0 to n-1 do test_expert i; done
  done;
  printf "Hedge weights: e^%a\n" (Array.print print_float5) w;
  Hedge (experts, w)


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
  let slices = train_rows / n in
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

let print_pred_conf oc (conf, pred) =
  fprintf oc "%d\t%f\n" pred conf

let print_cloned_head n ro ps = 
  let preds = Enum.clone ps |> Enum.take n |> Array.of_enum  in
  let data = Array2.sub_right train_data ro n in
  for i = 0 to n-1 do
    printf "Data:%a\npred:%a\n" 
      (Array.print print_float) (get_i data (i+1)|> Array1.to_array) 
      print_pred_conf preds.(i)
  done

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


let rand_slice ?(range=train_rows) len = (1+Random.int (range-len), len)

let grouped_slice len = 
  let per_group = len / category_count + 1 in
  group_by_cat per_group |> Array.to_list |> Array.concat

let test_accuracy ?(n=10_000) name cpred =
  marshal_file name cpred;
  let t0 = Sys.time () in
  let off,len = rand_slice ~range:(Array2.dim2 train_data) n in
  let right = 
    Array2.sub_right train_data off len 
  |> predict_many (predict_cat cpred)
  |> Array.fold_lefti (fun acc i (_,l) -> if train_labels.{i+1} = l then acc+1 else acc) 0     
  in
  let accuracy = 100. *. float right /. float n in
  printf "%s Accuracy: %d of %d (%.2f%%) (%.2fs)\n%!" name right n accuracy (Sys.time () -. t0);
  accuracy

let run_test ?(n = 10_000) name (cpred: cpredictor) = 
  marshal_file name cpred;
  let t0 = Sys.time () in
  let off,len = rand_slice ~range:(Array2.dim2 train_data) n in
  let eval = Array2.sub_right train_data off len |> predict_many (predict_cat cpred) |> Array.enum
(*  |> tap (print_cloned_head 10 rand_offset) *)
  |> best_predictions 
  |> tap (fun ps -> File.with_file_out ("bests/" ^ name ^ ".bests") (fun oc -> output_preds oc ps))
  |> evaluate 
  in
  printf "%s Overall: %.2f (%.2fs)\n%!" name eval (Sys.time () -. t0);
  eval

let avg ~n f x = (1--n |> Enum.map (fun _ -> f x) |> Enum.reduce (+.)) /. float n

let crossval ~n ~f ~xs =
  let test x = f x |> run_test ("crossvalidate:" ^ string_of_float x) in
  List.iter (avg ~n test |- printf "Avg: %f") xs

let crossval_int ~n ~f ~xs =
  let test x = f x |> run_test ("crossvalidate:" ^ string_of_int x) in
    List.iter (avg ~n test |- printf "Avg: %f") xs

let cv_one_one cats stringify ~loops ~f ~xs =
  let ca = 1 + Random.int category_count in
  let cb = 1 + Random.int category_count in

  let alen = Array.length cats.(ca) in
  let blen = Array.length cats.(cb) in
  printf "CV 1-1 cats %d(%d) and %d(%d)..\n%!" ca alen cb blen;
  let splita = alen * 4 / 5 in
  let cata, (testa: vec array) = 
    Array.sub cats.(ca) 0 splita,
    Array.sub cats.(ca) splita (alen - splita) |> Array.map (get_i train_data) in
  let splitb = blen * 4 / 5 in
  let catb, testb = 
    Array.sub cats.(cb) 0 splitb, 
    Array.sub cats.(cb) splitb (blen - splitb) |> Array.map (get_i train_data) in
  let test_count = Array.length testa + Array.length testb in
  let data = Array.append cata catb in
  let labels = Array.init (Array.length data) (fun i -> if i < Array.length cata then 1. else -1.) in
  List.iter (fun x -> 
    printf "x=%s " (stringify x);
    let run_test () = 
      let p = f x data labels |> predict_b in
      let wrong_a = ref 0 in
      Array.iter (fun x -> if p x <= 0. then incr wrong_a) testa;
      let wrong_b = ref 0 in
      Array.iter (fun x -> if p x >= 0. then incr wrong_b) testb;
      float (!wrong_a + !wrong_b) /. float test_count
    in
    let t1 = Sys.time () in
    let avg_wrong = avg ~n:loops run_test () in
    let t2 = Sys.time () in
    printf "Acc: %.1f%% in %.2fs\n%!" (100. *. (1. -. avg_wrong)) ((t2 -. t1) /. float loops);
  ) xs

let slice_shuffle n slice_len = 
  let slices = n / slice_len in
  1--slices |> Random.shuffle |> Array.enum |> Enum.map (fun i -> i --^ (i+slice_len)) |> Enum.flatten |> Array.of_enum

let predict p oc =
  predict_many (predict_cat p) train_data |> Array.iter (print_pred_conf oc)

let pred_read fn =
  let ic = Pervasives.open_in_bin fn in
  let (p: cpredictor) = Legacy.Marshal.from_channel ic in
  Pervasives.close_in ic;
  printf "Loaded predictor: %a\n%!" print_cpred p;
  p
