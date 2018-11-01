
for i in `seq 0 10 101`;
do
	tmp=$(($i + 10))


	run_name="#PBS -N dens_field_${i}_512"
	
	loop_name="for((i = ${i}; i < ${tmp}; i++))"

	sed "7s|.*|${run_name}|" base_create_grids.pbs > runscripts/tmp
	sed "16s|.*|${loop_name}|" runscripts/tmp > runscripts/create_grids_${i}_512actual.pbs 

	qsub runscripts/create_grids_${i}_512actual.pbs
done
