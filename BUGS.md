**BUG1**

`./bin/model_monostar 10 5 -d half_pi_with -j 1`
`./bin/model_monostar 20 7 -m af -d half_pi_with  --hamiltonian_J_classical 1.7 --hamiltonian_J_quantum 0.7`

Nan in results as LA fails

**BUG2**

`./bin/model_monostar 10 2 -i -a -r eg -m fo`
assert fails (debug build)
