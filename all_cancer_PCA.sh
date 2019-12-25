# CANCER=("acc" "blca" "cesc" "chol" "coadread" "dlbc" "esca" "gbm" "kich" "kirc" "kirp" "lgg" "lihc" "luad" "lusc" "meso" "ov" "paad" "pcpg" "prad" "sarc" "tgct" "thca" "thym" "ucs" "uvm")
CANCER=("acc" "blca" "cesc" "chol" "dlbc" "esca" "gbm" "kich" "kirc" "kirp" "lgg" "lihc" "luad" "lusc" "meso" "ov" "paad" "pcpg" "prad" "sarc" "tgct" "thca" "thym" "ucs" "uvm")
CANCER2=("brca" "coadread" "hnsc" "skcm" "stad" "ucec")

for i in `seq 1 $((${#CANCER[@]}-1))`
do
	# for j in `seq $(($i+1)) ${#CANCER[@]}`
	for j in `seq 1 ${#CANCER2[@]}`

	do
	# echo "${CANCER[i]}" "${CANCER[j]}"
	# echo ${CANCER[i]}  ${CANCER[j]}
	Rscript tSNE_2data.R ${CANCER[i]}  ${CANCER2[j]}
	done
done

for i in `seq 1 $((${#CANCER2[@]}-1))`
do
	for j in `seq 1 ${#CANCER2[@]}`
	do
	Rscript tSNE_2data.R ${CANCER2[i]}  ${CANCER2[j]}
	done
done

# "acc" "blca" 
# #"brca" 
# "cesc" "chol" 
# #"coadread" 
# "dlbc" "esca" "gbm"
# #"hnsc"
# "kich" "kirc" "kirp" "lgg" "lihc" "luad" "lusc" "meso" "ov" "pcpg" "paad" "prad" "sarc"
# #"skcm"
# #"stad"
# "tgct" "thca" "thym"
# #"ucec"
# "ucs" "uvm"