#bin/bash

model_type="af"
for n_sites in `seq 20 20`; do
    for n_pt in `seq 1 7`; do
        model_monostar -m "${model_type}" -n "${n_sites}" -p "${n_pt}" > "model_monostar_${model_type}_${n_sites}_${n_pt}.log"
    done
done
