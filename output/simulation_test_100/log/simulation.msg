RStudio not running
Traceback:
14: as.data.table.list(x, keep.rownames = keep.rownames, check.names = check.names, 
        .named = nd$.named)
13: data.table(V1 = i)
12: `[.data.table`(res$table.hosts, Move.ID[i], `:=`(current.in = going.to))
11: res$table.hosts[Move.ID[i], `:=`(current.in = going.to)]
10: moveFunction(res, pres.time, moving.full, structure.matrix = structure.matrix, 
        sdMoveParsed = sdMoveParsed, structure.raster = structure.raster, 
        attracted.by.raster = attracted.by.raster, max.raster = max.raster)
9: makeMoves(res$host.info.A, pres.time, moving.full, structure.matrix = structure.matrix)
8: singleDiscrete(...)
7: nosoiSim(type = "single", popStructure = "discrete", length.sim = sim_length, 
       max.infected = max_infections, init.individuals = 1, init.structure = "Denmark", 
       structure.matrix = transition_matrix, pExit = p_exit_func, 
       param.pExit = list(t_incubation = t_incub_func), timeDep.pExit = FALSE, 
       diff.pExit = FALSE, pMove = p_move_func, param.pMove = NA, 
       timeDep.pMove = FALSE, diff.pMove = FALSE, nContact = n_contact_func, 
       param.nContact = NA, timeDep.nContact = FALSE, diff.nContact = FALSE, 
       pTrans = p_trans_func_diff, param.pTrans = list(p_max = p_max_func, 
           t_incubation = t_incub_func), timeDep.pTrans = FALSE, 
       diff.pTrans = TRUE, hostCount.pTrans = FALSE, prefix.host = "H", 
       print.progress = TRUE, print.step = 10)
6: run_nosoi(transition_matrix, max_infections = max_infections, 
       sim_length = sim_length, out_dir = simulation_dir)
5: doTryCatch(return(expr), name, parentenv, handler)
4: tryCatchOne(expr, names, parentenv, handlers[[1L]])
3: tryCatchList(expr, classes, parentenv, handlers)
2: tryCatch({
       run_nosoi(transition_matrix, max_infections = max_infections, 
           sim_length = sim_length, out_dir = simulation_dir)
   }, error = function(e) {
       msg <- paste0("Attempt ", attempts, " failed with error: ", 
           conditionMessage(e))
       log_print(msg)
       return(NULL)
   })
1: run_simulation(max_infections = max_infections, sim_length = sim_length, 
       run_no = run_no, out_dir = out_dir)

