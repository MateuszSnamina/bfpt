#bin/bash

model_type="af"
for n_sites in `seq 20 20`; do
    for n_pt in `seq 22 24`; do
        model_monostar -m "${model_type}" -n "${n_sites}" -p "${n_pt}" > "model_monostar_${model_type}_${n_sites}_${n_pt}.log"
    done
done

#model_type="af"
#for n_sites in `seq 24 24`; do
#    for n_pt in `seq 0 15`; do
#        model_monostar -m "${model_type}" -n "${n_sites}" -p "${n_pt}" > "model_monostar_${model_type}_${n_sites}_${n_pt}.log"
#    done
#done
