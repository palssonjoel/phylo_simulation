RStudio not running
Traceback:
17: (function (t, current.in, p_max, t_incubation) 
    {
        R <- NULL
        if (current.in == "Albania") 
            R <- R_local[current.in]
        if (current.in == "Austria") 
            R <- R_local[current.in]
        if (current.in == "Belgium") 
            R <- R_local[current.in]
        if (current.in == "Bulgaria") 
            R <- R_local[current.in]
        if (current.in == "Croatia") 
            R <- R_local[current.in]
        if (current.in == "Cyprus") 
            R <- R_local[current.in]
        if (current.in == "Czechia") 
            R <- R_local[current.in]
        if (current.in == "Denmark") 
            R <- R_local[current.in]
        if (current.in == "Estonia") 
            R <- R_local[current.in]
        if (current.in == "Finland") 
            R <- R_local[current.in]
        if (current.in == "France") 
            R <- R_local[current.in]
        if (current.in == "Germany") 
            R <- R_local[current.in]
        if (current.in == "Greece") 
            R <- R_local[current.in]
        if (current.in == "Hungary") 
            R <- R_local[current.in]
        if (current.in == "Ireland") 
            R <- R_local[current.in]
        if (current.in == "Italy") 
            R <- R_local[current.in]
        if (current.in == "Latvia") 
            R <- R_local[current.in]
        if (current.in == "Lithuania") 
            R <- R_local[current.in]
        if (current.in == "Luxembourg") 
            R <- R_local[current.in]
        if (current.in == "Malta") 
            R <- R_local[current.in]
        if (current.in == "Netherlands") 
            R <- R_local[current.in]
        if (current.in == "Poland") 
            R <- R_local[current.in]
        if (current.in == "Portugal") 
            R <- R_local[current.in]
        if (current.in == "Romania") 
            R <- R_local[current.in]
        if (current.in == "Serbia") 
            R <- R_local[current.in]
        if (current.in == "Slovakia") 
            R <- R_local[current.in]
        if (current.in == "Slovenia") 
            R <- R_local[current.in]
        if (current.in == "Spain") 
            R <- R_local[current.in]
        if (current.in == "Sweden") 
            R <- R_local[current.in]
        if (current.in == "Switzerland") 
            R <- R_local[current.in]
        if (is.null(R)) {
            stop("Unknown state: ", current.in)
        }
        if (t < t_incubation) {
            return(0)
        }
        return(p_max * R)
    })(t = 9L, current.in = "Poland", p_max = 0.275051111523576, 
        t_incubation = 1.18955932714837)
16: do.call(pFunc, c(list(t = t), x))
15: (function (prestime, inf.time, ...) 
    {
        t = prestime - inf.time
        x <- list(...)
        do.call(pFunc, c(list(t = t), x))
    })(prestime = 98L, inf.time = 89L, current.in = "Poland", p_max = 0.275051111523576, 
        t_incubation = 1.18955932714837)
14: do.call(pFunc_eval, c(list(prestime = prestime), parameters))
13: pasedFunction$vect(prestime = pres.time, .SD)
12: `[.data.table`(res$table.hosts, active.hosts, pasedFunction$vect(prestime = pres.time, 
        .SD), by = "hosts.ID", .SDcols = pasedFunction$vectArgs)
11: res$table.hosts[active.hosts, pasedFunction$vect(prestime = pres.time, 
        .SD), by = "hosts.ID", .SDcols = pasedFunction$vectArgs]
10: applyFunctionToHosts(res, pres.time, pTransParsed, active.hosts)
9: meetTransmit(res$host.info.A, pres.time, positions = c("current.in", 
       "host.count.A"), nContactParsed, pTransParsed)
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
       print.progress = TRUE, print.step = 100)
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

