RStudio not running
Traceback:
7: withCallingHandlers(for (j in seq_along(x)) {
       x[[j]] <- vectbl_assign(x[[j]], i, recycled_value[[j]])
   }, vctrs_error = function(cnd) {
       vectbl_recycle_rhs_rows(value, length(i), i_arg, value_arg, 
           call)
       abort_assign_incompatible_type(x, recycled_value, j, value_arg, 
           cnd, call = call)
   })
6: tbl_subassign_row(xo, i, value, i_arg, value_arg, call)
5: tbl_subassign(x, i, j, value, i_arg, j_arg, substitute(value))
4: `[<-.tbl_df`(`*tmp*`, counter, , value = list(parent = 1029146L, 
       node = 450L, branch.length = 1L, label = "H-450", host = "H-450", 
       state = "Denmark", state.x = NA, state.y = NA, time.parent = 73L, 
       time = 74L))
3: `[<-`(`*tmp*`, counter, , value = list(parent = 1029146L, node = 450L, 
       branch.length = 1L, label = "H-450", host = "H-450", state = "Denmark", 
       state.x = NA, state.y = NA, time.parent = 73L, time = 74L))
2: getTransmissionTree(trans_simulation)
1: run_simulation(max_infections = max_infections, sim_length = sim_length, 
       run_no = run_no, out_dir = out_dir)

