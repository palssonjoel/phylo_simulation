RStudio not running
Traceback:
8: stopifnot(length(n) == 1L)
7: sample.int(length(x), size, replace, prob)
6: sample(bases, size = 1, prob = P_t[b, ])
5: FUN(X[[i]], ...)
4: vapply(infector_seq, function(b) {
       sample(bases, size = 1, prob = P_t[b, ])
   }, character(1))
3: evolve(host_table, E, E_1, eigvals, bases)
2: hky_nosoi(ref_genome = ref_genome, host_data = host_data)
1: run_simulation(max_infections = max_infections, sim_length = sim_length, 
       seed = seed, out_dir = out_dir)

