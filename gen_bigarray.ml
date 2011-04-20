open Batteries_uni
open Printf
open Bigarray

let () = if Array.length Sys.argv < 2 then printf "Usage: %s [input_file]" Sys.argv.(0)

let use_std_trans = true

let cols = Datafile.cols
let category_count = 164

let gen_output input_file = 
  let output_file = input_file ^ ".ba64" in
  let rows = Sys.argv.(2) |> int_of_string in
  let mmap = Datafile.write64 rows output_file in
  let maxs = Array.create cols min_float in
  let push_x row col x =
    if row land 0xfff = 0 then print_string ".";
    mmap.{col+1,row+1} <- x;
    if x > maxs.(col) then maxs.(col) <- x;
  in  
  let data = Datafile.read_text input_file in
  Enum.iteri (fun row d -> List.iteri (push_x row) d; mmap.{cols, row+1} <- 1.) data;
  let t0 = Sys.time () in
  printf "Reading took: %.2f s\n" t0;
  let maxs = if use_std_trans then (printf "Using standard scaling\n"; Datafile.scaling) else 
      (Array.print ~last:"|]\n" (fun oc f -> fprintf oc "%2.0f" f) stdout maxs; maxs) in
  IO.flush_all ();
  for j = 1 to rows do
    for i = 1 to cols do
      mmap.{i,j} <- mmap.{i,j} /. maxs.(i-1)
    done;
  done;
  printf "Scaling took: %.2f s\n" (Sys.time () -. t0);
  ()
    
let group_by_cat train_labels ?(category_count=category_count) cap =
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
  for i = 1 to Array1.dim train_labels do
    let l = train_labels.{i} in 
    if l < category_count then insert_random l i;
  done;
  let cap_slice v = 
    Vect.to_array (if Vect.length v > cap then Vect.sub 0 cap v else v)
  in
  Array.map cap_slice cats

let gen_slice32 n input_file input_labels_file =
  let output_file = input_file ^ ".slice" in
  let output_labels = File.open_out (input_labels_file ^ ".slice") in
  let indata = Datafile.get_matrix input_file in
  let cs = group_by_cat (Datafile.read_label_file input_labels_file) n in
  let rows = Array.fold_left (fun acc ci -> acc + Array.length ci) 0 cs in
  let mmap = Datafile.write rows output_file in
  let row = ref 1 in
  let push l i =
    if !row land 0xfff = 0 then print_string ".";
    let src = Array2.slice_right indata i in
    let dst = Array2.slice_right mmap !row in
    Array1.blit src dst;
    incr row;
    fprintf output_labels "%d\n" l;
  in
  Array.iteri (fun l -> Array.iter (fun i -> push l i)) cs;
  let t0 = Sys.time () in
  printf "Slicing took: %.2f s\n" t0;
  ()

let gen_slice64 n input_file input_labels_file =
  let output_file = input_file ^ ".slice" in
  let output_labels = File.open_out (input_labels_file ^ ".slice") in
  let indata = Datafile.get_matrix64 input_file in
  let cs = group_by_cat (Datafile.read_label_file input_labels_file) n in
  let rows = Array.fold_left (fun acc ci -> acc + Array.length ci) 0 cs in
  let mmap = Datafile.write64 rows output_file in
  let row = ref 1 in
  let push l i =
    if !row land 0xfff = 0 then printf ".%!";
    let src = Array2.slice_right indata i in
    let dst = Array2.slice_right mmap !row in
    Array1.blit src dst;
    incr row;
    fprintf output_labels "%d\n" l;
  in
  Array.iteri (fun l -> Array.iter (fun i -> push l i)) cs;
  let t0 = Sys.time () in
  printf "Slicing took: %.2f s\n" t0;
  ()

(* let () = gen_output Sys.argv.(1) *)
let () = gen_slice64 100 Sys.argv.(1) Sys.argv.(2)

