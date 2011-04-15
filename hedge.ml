open Batteries_uni
open Printf
open Predict

let n = 1000

let () = 
  let experts = args () |> map pred_read |> Array.of_enum in
  hedge experts (rand_slice ~range:train_rows n) |> test_accuracy ~n "Hedge" |> ignore
