open Batteries_uni
open Predict

let () = 
(*  kperceptron3_slice 1 2000 |> extend_hamm |> marshal_file "hamm_kp3_0_2k"; *)
(*  train_slices 1000 kperceptron3_slice |> List.iter (extend_hamm |- marshal_file "hamm_kp3_slc_1k"); *)
(* DONE  kperceptron_slice (gt_rbf 50.) kp_core 100 1 2000 |> extend_hamm 32 |> tap (marshal_file "hamm_kpr_1_2k") |> run_test "hamm_kpr_1_2k"; *)
(* DONE  train_slices 1000 (kperceptron_slice (gt_rbf 50.) kp_core 100) |> List.iter (extend_hamm 32 |- tap (marshal_file "hamm_kpr_slc_1k") |- run_test "hamm_kpr_slc_1k") *)
(* DONE  Array.init train_rows (fun i -> i+1) |> batch_perb |> extend_hamm 32 |> run_test "hamm32_percb_full" |> ignore; *)

  svm |> extend_one_one (group_by_cat 500) (function _ -> true) |> run_test "svm_OVO" |> ignore;
(*
  Bigarray.Array1.to_array d64_labels |> Array.map float |> vec_of_arr |> svm |> run_test "svmslice" |> ignore *)
(*
 kperceptron_elems (gt_rbf 1.0) (ft_core ~b:200 ~b0:20 ~bincr:1) ~loops:1000
  |> extend_one_one (group_by_cat 500) (function (*| Kern_rbf (_,os,_) -> Array.length os >= 200*) | _ -> true)
  |> run_test "kprbf.4_ft200_100loops_1-1cap500" |> ignore
*)

(*
  kperceptron_offs (gt_rbf 0.4) (ft_core ~b:4000 ~b0:4000) ~loops:1000 (grouped_slice 10000)
  |> extend_hamm 60
  |> run_test "kprbf.4_ft400_50loops_hamm60" |> ignore
*)

(*  slice_shuffle train_rows 2048 |> batch_perb |> extend_hamm 32 |> run_test "hamm32_percb_full" ~n:100_000 |> ignore *)

