open Batteries_uni
open Printf

let () = if Array.length Sys.argv < 2 then printf "Usage: %s [input_file]" Sys.argv.(0)

let use_std_trans = true
let sliced = true

let input_file = Sys.argv.(1)
let cols = Datafile.cols

let gen_output () = 
  let output_file = input_file ^ ".ba" in
  let rows = Sys.argv.(2) |> int_of_string in
  let mmap = Datafile.write rows output_file in
  let maxs = Array.create cols min_float in
  let push_x row col x =
    if row land 0xff = 0 then print_string ".";
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
    
(*
let gen_sliced () = 
  let label_data = read_label_file Sys.argv.(2) in

  let num_cats = 164 in
  let label_files = Array.init num_cats (fun i -> input_file ^ ".ba." ^ (string_of_int i)) in
  let ocs = Array.map open_out_bin label_files in
*)
    
  
  
  
  
let () = if sliced then gen_sliced () else gen_output ()
