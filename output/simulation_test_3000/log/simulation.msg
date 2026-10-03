RStudio not running
Traceback:
19: lapply(args, as.name)
18: lapply(args, as.name)
17: standardGeneric("pmax")
16: pmax(0, rbeta(x, shape1 = 3, shape2 = 5))
15: x(1)
14: FUN(X[[i]], ...)
13: lapply(X = X, FUN = FUN, ...)
12: sapply(ParamHost, function(x) x(1))
11: as.list(sapply(ParamHost, function(x) x(1)))
10: newLine(hosts.ID, as.character(df.meetTransmit[[i, "active.hosts"]]), 
        infected.in = getPositionInfected(res, df.meetTransmit, i), 
        pres.time, ParamHost, current.environmental.value = df.meetTransmit[[i, 
            "current.env.value"]], current.cell.number.raster = df.meetTransmit[[i, 
            "current.cell.raster"]], current.count.A = df.meetTransmit[[i, 
            "host.count.A"]], current.count.B = df.meetTransmit[[i, 
            "host.count.B"]])
9: writeInfected(df.meetTransmit, res$host.info.A, pres.time, ParamHost)
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

