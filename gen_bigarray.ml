open Batteries_uni
open Printf

let () = if Array.length Sys.argv < 2 then printf "Usage: %s [input_file]" Sys.argv.(0)

let input_file = Sys.argv.(1)
let output_file = input_file ^ ".ba"

let use_std_trans = false

let cols = Datafile.cols
let rows = Datafile.rows

let () = 
  let mmap = Datafile.write output_file in
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
  if use_std_trans then (
    printf "Using standard scaling\n";
    List.iteri (fun i x -> maxs.(i) <- x) Datafile.scaling;
  ) else
    Array.print ~last:"|]\n" (fun oc f -> fprintf oc "%2.0f" f) stdout maxs;
  IO.flush_all ();
  for j = 1 to rows do
    for i = 1 to cols do
      mmap.{i,j} <- mmap.{i,j} /. maxs.(i-1)
    done;
  done;
  printf "Scaling took: %.2f s\n" (Sys.time () -. t0);
  ()