open Batteries_uni
open Predict

let () = 
  (* best: 2K?
     crossval_int (fun i ->   
     let rand_offset = Random.int (Array2.dim2 train_data - i) in
     kperceptron3_slice rand_offset i |> extend_hamm) [500; 1000; 2000; 4000; 8000; 10000; 15000; 20000]
  *)

 
  (*kperceptron_elems gt_k3 kp_core 100 |> extend_one_one |> tap (marshal_file "oneone_kp3") |> run_test "kp3_5k_1-1"; *)

(*  crossval ~n:3 ~xs:[0.01; 0.05; 0.1; 0.5; 1.; 2.; 4.; 7.; 9.; 10.; 100.]
    ~f:(fun i -> kperceptron_slice (gt_rbf i) kp_core 100 (rand_slice 2000) |> extend_hamm 50);


  crossval ~n:3 ~xs:[0.01; 0.05; 0.1; 0.5; 1.; 2.; 4.; 7.; 9.; 10.; 100.]
    ~f:(fun i -> kperceptron_slice (gt_rbf i) kp_core 100 (rand_slice 4000) |> extend_hamm 50);
*)

(*  extend_one_one (1,train_rows) perceptron_b ~cap:500 |> run_test "perb_full_1-1" |> ignore; *)


(*  cv_one_one (group_by_cat 625) string_of_float ~loops:5
    ~xs:[0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1.0; 1.2; 1.5; 2.0; 5.0] 
    ~f:(fun sigma -> klsq ~lambda:1.0 ~sigma);*)

  cv_one_one (group_by_cat 625) string_of_float ~loops:1
    ~xs:[0.01; 0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1.0; 1.2; 1.5; 2.0; 5.0] 
    ~f:(fun s -> klsq ~lambda:0.1 ~sigma:0.4 ~cutoff:s);

(* SEARCH FOR PARAMETERS *)
(* find the best sigma (0.3-1.0 is good) *)
(*  cv_one_one (group_by_cat 625) string_of_float ~loops:5
    ~xs:[0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1.0; 1.2; 1.5; 2.0; 5.0] 
    ~f:(fun s -> kperceptron_elems (gt_rbf s) (ft_core ~b:500) ~loops:250); *)

(* find a good bound for forgetron memory (10 is enough)  *)
(*
  cv_one_one (group_by_cat 625) string_of_int ~loops:5
    ~xs:[10; 15; 20; 25; 30; 35; 40; 50; 60; 80; 100; 150] 
    ~f:(fun n -> kperceptron_elems (gt_rbf 0.4) (ft_core ~b:max_int ~b0:n) ~loops:250); 
*)
(* Find a good number of perceptron iterations (40 passes costs no more time than one pass, it seems)
  cv_one_one ~cap:500 string_of_int ~loops:5
    ~xs:[1; 4; 10; 40; 100; 250; 500; 1000; 2500; 4000; 10000]
    ~f:(fun n -> kperceptron_elems (gt_rbf 0.4) (ft_core ~b:50) ~loops:n);
 *)
(*
  kperceptron_elems (gt_rbf 0.1) (ft_core ~b:300) ~loops:100 |> extend_one_one ~cap:500 (1, train_rows) |> run_test "kprbf0.1_ft300_100loops_1-1cap500" |> ignore;
*)

(*  crossval ~n:3 ~xs:[0.01; 0.05; 0.1; 0.5; 1.; 2.; 4.; 7.; 9.; 10.]
    ~f:(fun i -> kperceptron_slice (gt_rbf i) (ft_core ~b:1000) ~loops:1000 (rand_slice 4000) |> extend_hamm 50);*)


  crossval_int ~n:3 ~xs:(8--31 |> List.of_enum)
    ~f:(fun i -> kperceptron_offs (gt_rbf 0.3) (ft_core ~b:150 ~b0:10) ~loops:100 (grouped_slice 10000) |> extend_hamm i);
